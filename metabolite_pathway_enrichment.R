#!/usr/bin/env Rscript
# Metabolite Pathway Enrichment Analysis
#
# Input:  a CSV/TSV/Excel file with columns:
#           name  – compound names (strings or numeric feature IDs)
#           mz    – observed m/z
#           rt    – retention time
#           sample columns named (plasma_)?(NC|NN|TC|TN)<replicate_number>
#
# Conditions:
#   NC and NN are treated as the same baseline → pooled as "Control"
#   TC and TN are the two treatment arms
#   The "plasma_" prefix on column names is stripped and ignored
#
# Output: pathway_enrichment_dotplot.pdf
#         pathway_enrichment_results.csv

# ── 0. Null-coalescing helper ──────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── 1. Package bootstrap ───────────────────────────────────────────────────────
bioc_pkgs <- c("KEGGREST")
cran_pkgs <- c("tidyverse", "readxl", "httr", "jsonlite", "ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

for (p in bioc_pkgs)
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)

for (p in cran_pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(KEGGREST)
  library(httr)
  library(jsonlite)
  library(ggplot2)
})

# ── 2. Data loading ────────────────────────────────────────────────────────────
load_data <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))
  data <- switch(ext,
    csv  = read_csv(file_path, show_col_types = FALSE),
    tsv  = ,
    txt  = read_tsv(file_path, show_col_types = FALSE),
    xlsx = ,
    xls  = read_excel(file_path),
    stop("Unsupported format '", ext, "'. Use csv, tsv, txt, xlsx, or xls.")
  )
  message(sprintf("Loaded: %d compounds × %d columns", nrow(data), ncol(data)))
  data
}

# ── 3. Sample parsing & condition grouping ─────────────────────────────────────
# NC and NN are pooled → "Control"; TC and TN remain separate treatment arms.
parse_samples <- function(data) {
  cols <- colnames(data)
  # Pattern handles:
  #   plain conditions:      NC1, NN2, TC3, TN4, C1
  #   plasma_ prefix:        Plasma_NC1, Plasma_C2
  #   technical-rep suffix:  Plasma_NC2_1, Plasma_TC2_2
  # NC, NN, and bare C are all pooled as "Control"
  pat  <- "^(?:plasma_)?(NC|NN|TC|TN|C)(\\d+)(_\\d+)?$"
  sample_cols <- cols[grepl(pat, cols, perl = TRUE, ignore.case = TRUE)]

  if (!length(sample_cols))
    stop(
      "No sample columns found.\n",
      "Expected names like NC1, NN2, TC1, C3, plasma_TN3, Plasma_NC2_1, etc.\n",
      "Columns seen: ", paste(head(cols, 20), collapse = ", ")
    )

  meta <- tibble(raw = sample_cols) %>%
    mutate(
      stripped  = sub("(?i)^plasma_", "", raw, perl = TRUE),
      grp_raw   = sub("(?i)^(NC|NN|TC|TN|C)(\\d+)(_\\d+)?$", "\\1", stripped, perl = TRUE),
      replicate = as.integer(sub("(?i)^(NC|NN|TC|TN|C)(\\d+)(_\\d+)?$", "\\2", stripped, perl = TRUE)),
      condition = case_when(
        toupper(grp_raw) %in% c("NC", "NN", "C") ~ "Control",
        toupper(grp_raw) == "TC"                  ~ "TC",
        toupper(grp_raw) == "TN"                  ~ "TN",
        TRUE                                       ~ grp_raw
      )
    )

  message("Sample → condition mapping:")
  print(count(meta, grp_raw, condition))

  list(sample_cols = sample_cols, meta = meta)
}

# ── 4. Average intensities per condition ───────────────────────────────────────
average_intensities <- function(data, sample_info) {
  find_col <- function(pattern)
    grep(pattern, colnames(data), ignore.case = TRUE, value = TRUE)[1]

  name_col <- find_col("^name$")
  mz_col   <- find_col("^(mz|m\\.z|m/z)")   # no $ — matches "m/z meas." etc.
  rt_col   <- find_col("^rt\\b|^retention")  # matches "RT [min]" and "retention time"

  if (any(is.na(c(name_col, mz_col, rt_col))))
    stop("Could not locate name/mz/rt columns. Found: ", paste(colnames(data), collapse = ", "))

  data %>%
    select(compound = !!name_col, mz = !!mz_col, rt = !!rt_col,
           all_of(sample_info$sample_cols)) %>%
    mutate(compound = as.character(compound)) %>%
    pivot_longer(all_of(sample_info$sample_cols),
                 names_to = "raw", values_to = "intensity") %>%
    left_join(select(sample_info$meta, raw, condition), by = "raw") %>%
    group_by(compound, mz, rt, condition) %>%
    summarise(mean_intensity = mean(intensity, na.rm = TRUE), .groups = "drop")
}

# ── 5a. KEGG pathway mapping (name-based, with on-disk cache) ─────────────────
map_kegg <- function(avg, cache_file = "kegg_mapping_cache.rds") {
  cache <- if (file.exists(cache_file)) readRDS(cache_file) else list(hits = list(), pw_names = c())

  # Only string (non-numeric) names can be searched by name in KEGG
  search_names <- unique(avg$compound) %>%
    .[!grepl("^[0-9.\\s]+$", .)]

  message(sprintf("KEGG: querying %d named compounds (may take several minutes)…", length(search_names)))

  # Fetch the full HSA pathway name table once
  if (!length(cache$pw_names)) {
    message("  Fetching KEGG pathway list…")
    cache$pw_names <- tryCatch({
      pw <- keggList("pathway", "hsa")
      setNames(sub("\\s*-\\s*Homo sapiens.*$", "", pw), names(pw))
    }, error = function(e) c())
  }

  for (i in seq_along(search_names)) {
    nm <- search_names[i]
    if (!is.null(cache$hits[[nm]])) next

    Sys.sleep(0.35)  # respect KEGG rate limit

    kegg_id <- tryCatch({
      res <- keggFind("compound", nm)
      if (length(res)) names(res)[1] else NA_character_
    }, error = function(e) NA_character_)

    pathways <- character(0)
    if (!is.na(kegg_id)) {
      Sys.sleep(0.35)
      pathways <- tryCatch({
        lk <- keggLink("pathway", kegg_id)
        unname(lk[grepl("^path:(hsa|map)", lk)])
      }, error = function(e) character(0))
    }

    cache$hits[[nm]] <- pathways

    if (i %% 20 == 0) {
      message(sprintf("  KEGG: %d / %d", i, length(search_names)))
      saveRDS(cache, cache_file)
    }
  }

  saveRDS(cache, cache_file)

  map_df <- bind_rows(lapply(names(cache$hits), function(nm) {
    pws <- cache$hits[[nm]]
    if (!length(pws)) return(NULL)
    tibble(
      compound     = nm,
      pathway_id   = pws,
      pathway_name = ifelse(pws %in% names(cache$pw_names), cache$pw_names[pws], pws),
      db           = "KEGG"
    )
  }))

  message(sprintf("KEGG: %d compounds → %d pathways",
                  n_distinct(map_df$compound), n_distinct(map_df$pathway_id)))
  map_df
}

# ── 5b. HMDB pathway mapping (name + m/z fallback, with on-disk cache) ────────
.hmdb_get <- function(url) {
  tryCatch({
    r <- GET(url, timeout(15), add_headers(Accept = "application/json"))
    if (status_code(r) != 200) return(NULL)
    fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  }, error = function(e) NULL)
}

hmdb_search_name <- function(nm) {
  .hmdb_get(paste0(
    "https://hmdb.ca/metabolites.json?filter%5Bname_query%5D=",
    URLencode(nm, reserved = TRUE)
  ))
}

hmdb_search_mz <- function(mz, tol_ppm = 10) {
  # Assume [M+H]+ adduct; adjust if your data uses a different ionisation mode
  mass <- mz - 1.007276
  tol  <- mass * tol_ppm / 1e6
  .hmdb_get(sprintf(
    "https://hmdb.ca/metabolites.json?filter%%5Bmass_min%%5D=%.5f&filter%%5Bmass_max%%5D=%.5f",
    mass - tol, mass + tol
  ))
}

hmdb_pathways <- function(hmdb_id) {
  dat <- .hmdb_get(paste0("https://hmdb.ca/metabolites/", hmdb_id, ".json"))
  if (is.null(dat)) return(character(0))
  pws <- dat$pathways
  if (!length(pws)) return(character(0))
  vapply(pws, function(p) p$name %||% "", character(1))
}

first_accession <- function(resp) {
  mets <- resp$metabolites
  if (!length(mets)) return(NA_character_)
  mets[[1]]$accession %||% NA_character_
}

map_hmdb <- function(avg, cache_file = "hmdb_mapping_cache.rds") {
  cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()

  named_cmpds <- unique(avg$compound) %>%
    .[!grepl("^[0-9.\\s]+$", .)]

  unknown_cmpds <- avg %>%
    filter(grepl("^[0-9.\\s]+$", compound)) %>%
    select(compound, mz) %>%
    distinct()

  message(sprintf("HMDB: %d named compounds + %d m/z lookups…",
                  length(named_cmpds), nrow(unknown_cmpds)))

  results <- list()

  # ── name-based ────────────────────────────────────────────────────────────
  for (i in seq_along(named_cmpds)) {
    nm  <- named_cmpds[i]
    key <- paste0("name:", nm)
    if (!is.null(cache[[key]])) { results[[nm]] <- cache[[key]]; next }

    Sys.sleep(0.6)
    resp    <- hmdb_search_name(nm)
    hmdb_id <- if (!is.null(resp)) first_accession(resp) else NA_character_

    pws <- if (!is.na(hmdb_id)) { Sys.sleep(0.6); hmdb_pathways(hmdb_id) } else character(0)

    results[[nm]] <- tibble(
      compound     = nm,
      hmdb_id      = hmdb_id,
      pathway_name = if (length(pws)) pws else character(0),
      db           = "HMDB"
    )
    cache[[key]] <- results[[nm]]

    if (i %% 10 == 0) {
      message(sprintf("  HMDB names: %d / %d", i, length(named_cmpds)))
      saveRDS(cache, cache_file)
    }
  }

  # ── m/z-based for numeric/unknown features ────────────────────────────────
  for (i in seq_len(nrow(unknown_cmpds))) {
    cmpd  <- as.character(unknown_cmpds$compound[i])
    mz_v  <- unknown_cmpds$mz[i]
    key   <- paste0("mz:", round(mz_v, 4))
    if (!is.null(cache[[key]])) { results[[cmpd]] <- cache[[key]]; next }

    Sys.sleep(0.6)
    resp    <- hmdb_search_mz(mz_v)
    hmdb_id <- if (!is.null(resp)) first_accession(resp) else NA_character_

    pws <- if (!is.na(hmdb_id)) { Sys.sleep(0.6); hmdb_pathways(hmdb_id) } else character(0)

    results[[cmpd]] <- tibble(
      compound     = cmpd,
      hmdb_id      = hmdb_id,
      pathway_name = if (length(pws)) pws else character(0),
      db           = "HMDB"
    )
    cache[[key]] <- results[[cmpd]]

    if (i %% 10 == 0) {
      message(sprintf("  HMDB m/z: %d / %d", i, nrow(unknown_cmpds)))
      saveRDS(cache, cache_file)
    }
  }

  saveRDS(cache, cache_file)

  out <- bind_rows(results) %>%
    filter(!is.na(pathway_name), nchar(pathway_name) > 0)

  message(sprintf("HMDB: %d compounds → %d pathways",
                  n_distinct(out$compound), n_distinct(out$pathway_name)))
  out
}

# ── 6. Over-representation analysis (Fisher's exact test) ─────────────────────
compute_fc <- function(avg, treatment) {
  avg %>%
    filter(condition %in% c("Control", treatment)) %>%
    pivot_wider(id_cols     = c(compound, mz, rt),
                names_from  = condition,
                values_from = mean_intensity) %>%
    mutate(
      treat = .data[[treatment]],
      fc     = (treat + 1) / (Control + 1),   # pseudocount avoids Inf/0
      log2fc = log2(fc)
    )
}

ora <- function(fg, bg, pw_map) {
  fg  <- unique(fg)
  bg  <- unique(bg)
  N   <- length(bg)
  n   <- length(fg)
  pws <- unique(pw_map$pathway_name)

  res <- lapply(pws, function(pw) {
    pw_cpds <- unique(pw_map$compound[pw_map$pathway_name == pw])
    K <- sum(bg %in% pw_cpds)   # pathway hits in background
    k <- sum(fg %in% pw_cpds)   # pathway hits in foreground
    if (K == 0 || k == 0) return(NULL)

    # 2×2 contingency table for one-sided Fisher test (enrichment)
    mat <- matrix(c(k, n - k, K - k, N - n - K + k), nrow = 2)
    p   <- fisher.test(mat, alternative = "greater")$p.value

    tibble(
      pathway_name    = pw,
      n_pathway       = K,
      n_hit           = k,
      n_background    = N,
      n_foreground    = n,
      fold_enrichment = (k / n) / (K / N),
      p_value         = p
    )
  })

  bind_rows(res) %>%
    mutate(
      p_adj      = p.adjust(p_value, "BH"),
      neg_log10p = -log10(p_value)
    ) %>%
    arrange(p_value)
}

enrich_conditions <- function(avg, kegg_map, hmdb_map, fc_thresh = 1.5) {
  combined_pw <- bind_rows(
    select(kegg_map, compound, pathway_name),
    select(hmdb_map, compound, pathway_name)
  ) %>% distinct()

  if (!nrow(combined_pw))
    stop("No pathway mappings found. Check KEGG/HMDB connectivity and compound names.")

  background <- intersect(unique(avg$compound), unique(combined_pw$compound))
  message(sprintf("ORA background: %d compounds with pathway annotations", length(background)))

  treatment_arms <- intersect(c("TC", "TN"), unique(avg$condition))

  bind_rows(lapply(treatment_arms, function(cond) {
    fc_tbl <- compute_fc(avg, cond)

    fg <- fc_tbl %>%
      filter(abs(log2fc) >= log2(fc_thresh)) %>%
      pull(compound) %>%
      intersect(background)

    message(sprintf("%s vs Control: %d foreground compounds (|FC| ≥ %.1f, in pathway map)",
                    cond, length(fg), fc_thresh))

    # Relax threshold automatically if too few compounds pass
    if (length(fg) < 5) {
      message("  Fewer than 5 foreground compounds; relaxing |FC| threshold to 1.2")
      fg <- fc_tbl %>%
        filter(abs(log2fc) >= log2(1.2)) %>%
        pull(compound) %>%
        intersect(background)
      message(sprintf("  → %d foreground compounds after relaxation", length(fg)))
    }

    if (length(fg) < 3) {
      message("  Skipping enrichment for ", cond, " (< 3 foreground compounds)")
      return(NULL)
    }

    ora(fg, background, combined_pw) %>% mutate(condition = cond)
  }))
}

# ── 7. Dot plot ────────────────────────────────────────────────────────────────
dotplot_enrichment <- function(
    enr,
    top_n    = 20,
    p_cutoff = 0.05,
    out_file = "pathway_enrichment_dotplot.pdf"
) {
  plt <- enr %>%
    filter(p_value < p_cutoff, fold_enrichment > 1, n_hit >= 2) %>%
    group_by(condition) %>%
    slice_min(p_value, n = top_n, with_ties = FALSE) %>%
    ungroup()

  if (!nrow(plt)) {
    message(
      "No pathways pass p < ", p_cutoff, " with fold enrichment > 1.\n",
      "Showing top ", top_n, " pathways per condition (no significance filter)."
    )
    plt <- enr %>%
      filter(fold_enrichment > 1, n_hit >= 2) %>%
      group_by(condition) %>%
      slice_min(p_value, n = top_n, with_ties = FALSE) %>%
      ungroup()
  }

  if (!nrow(plt)) stop("No enrichment results to plot.")

  # Order pathways by maximum fold enrichment across conditions
  pw_order <- plt %>%
    group_by(pathway_name) %>%
    summarise(fe = max(fold_enrichment, na.rm = TRUE), .groups = "drop") %>%
    arrange(fe) %>%
    pull(pathway_name)

  plt <- plt %>%
    mutate(
      pathway_name = factor(pathway_name, levels = pw_order),
      dot_size     = pmin(-log10(p_value), 8),    # cap for visual clarity
      fdr_label    = if_else(p_adj < 0.05, "FDR<5%", "")
    )

  cond_pal <- c(TC = "#E64B35", TN = "#4DBBD5")

  p <- ggplot(plt, aes(x = fold_enrichment, y = pathway_name,
                       size = dot_size, colour = condition)) +
    geom_point(alpha = 0.80) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey55", linewidth = 0.4) +
    scale_size_continuous(
      name   = expression(-log[10](italic(p))),
      range  = c(2, 11),
      breaks = c(1, 2, 3, 5, 8),
      labels = c("0.1", "0.01", "0.001", "1e-5", "≥1e-8")
    ) +
    scale_colour_manual(values = cond_pal, name = "Condition") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(
      title    = "Metabolite Pathway Enrichment (ORA)",
      subtitle = "Treatment vs Control  |  dot size = −log₁₀(p-value)  |  KEGG + HMDB pathways",
      x        = "Fold Enrichment",
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y        = element_text(size = 8),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 9, colour = "grey40"),
      legend.position    = "right",
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )

  plot_height <- max(6, n_distinct(plt$pathway_name) * 0.35 + 3)
  ggsave(out_file, p, width = 13, height = plot_height, dpi = 300)
  message("Plot saved → ", out_file)
  print(p)
  invisible(p)
}

# ── 8. Main pipeline ───────────────────────────────────────────────────────────
run_pipeline <- function(
  file_path,
  fc_threshold = 1.5,     # minimum |fold change| to include compound in ORA foreground
  p_cutoff     = 0.05,    # p-value threshold for dot plot display
  top_n        = 20,      # max pathways shown per condition in the plot
  use_hmdb     = TRUE,    # set FALSE to skip HMDB (much faster; KEGG only)
  out_plot     = "pathway_enrichment_dotplot.pdf",
  out_table    = "pathway_enrichment_results.csv"
) {
  message("\n══ Loading data ══════════════════════════════════════════")
  dat <- load_data(file_path)

  message("\n══ Parsing sample columns ════════════════════════════════")
  si  <- parse_samples(dat)

  message("\n══ Averaging intensities by condition ════════════════════")
  avg <- average_intensities(dat, si)

  message("\n══ KEGG pathway mapping ══════════════════════════════════")
  kegg <- map_kegg(avg)

  hmdb <- tibble(compound = character(), pathway_name = character(), db = character())
  if (use_hmdb) {
    message("\n══ HMDB pathway mapping ══════════════════════════════════")
    hmdb <- map_hmdb(avg)
  }

  message("\n══ Enrichment analysis ═══════════════════════════════════")
  enr <- enrich_conditions(avg, kegg, hmdb, fc_threshold)

  write_csv(enr, out_table)
  message("Results table → ", out_table)

  message("\n══ Dot plot ══════════════════════════════════════════════")
  dotplot_enrichment(enr, top_n, p_cutoff, out_plot)

  message("\n══ Done ══════════════════════════════════════════════════")
  invisible(list(averaged = avg, kegg = kegg, hmdb = hmdb, enrichment = enr))
}

# ── Run ────────────────────────────────────────────────────────────────────────
# Change file_path to the path of your data file.
# All other arguments have sensible defaults.
results <- run_pipeline(
  file_path    = "your_data_file.csv",  # ← update this
  fc_threshold = 1.5,
  p_cutoff     = 0.05,
  top_n        = 20,
  use_hmdb     = TRUE
)
