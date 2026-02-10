# =======================
# Common-pathways GSEA runner
# Compare two expression datasets via GSEA and find:
#  - common significant pathways
#  - common UP (NES>0 in both)
#  - common DOWN (NES<0 in both)
# Exports: TSV tables + PNG/PDF plots + optional RDS
# =======================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(msigdbr)
  library(clusterProfiler)
  library(ggplot2)
  library(readr)
})

# ---- Helpers ----

collapse_rank_by_log2fc <- function(df, symbol_col = "symbol", lfc_col = "log2FoldChange") {
  df %>%
    filter(!is.na(.data[[lfc_col]]), !is.na(.data[[symbol_col]])) %>%
    mutate(
      !!symbol_col := as.character(.data[[symbol_col]]),
      !!lfc_col := as.numeric(.data[[lfc_col]])
    ) %>%
    group_by(.data[[symbol_col]]) %>%
    slice_max(order_by = abs(.data[[lfc_col]]), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(.data[[lfc_col]])) %>%
    select(all_of(symbol_col), all_of(lfc_col)) %>%
    deframe()
}

make_outdir <- function(outdir = ".", prefix = "Run", stamp = TRUE) {
  prefix_safe <- gsub("[^A-Za-z0-9_\\-]+", "_", prefix)
  if (isTRUE(stamp)) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    outdir <- file.path(outdir, paste0(prefix_safe, "_", ts))
  } else {
    outdir <- file.path(outdir, prefix_safe)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(outdir, winslash = "/", mustWork = FALSE)
}

# ---- Main function ----

run_common_gsea <- function(
    file1 = NULL,
    file2 = NULL,
    label1 = "Dataset1",
    label2 = "Dataset2",
    prefix = "Run",
    outdir = ".",
    timestamp_dir = TRUE,
    
    # Input columns
    symbol_col = "symbol",
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    
    # Filters for "sig genes" tables (not required for GSEA itself)
    padj_cutoff = 0.05,
    lfc_cutoff = 0.5,
    
    # GSEA parameters
    gsea_padj_cutoff = 0.1,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    
    # MSigDB parameters
    species = "Mus musculus",
    category = "H",
    subcategory = NULL,
    
    # Outputs
    plot_set = c("drivers", "common"),  # "drivers" = significant in both; "common" = all overlaps
    export_plots = TRUE,
    export_tables = TRUE,
    save_rds = TRUE
) {
  plot_set <- match.arg(plot_set)
  
  # ---- 0) Resolve inputs ----
  if (is.null(file1)) {
    message("Choose Dataset 1 Excel file...")
    file1 <- file.choose()
  }
  if (is.null(file2)) {
    message("Choose Dataset 2 Excel file...")
    file2 <- file.choose()
  }
  
  prefix_safe <- gsub("[^A-Za-z0-9_\\-]+", "_", prefix)
  label1_safe <- gsub("[^A-Za-z0-9_\\-\\s]+", "_", label1)
  label2_safe <- gsub("[^A-Za-z0-9_\\-\\s]+", "_", label2)
  
  outdir_final <- make_outdir(outdir = outdir, prefix = prefix_safe, stamp = timestamp_dir)
  message("Outputs will be written to: ", outdir_final)
  
  # ---- 1) Read ----
  df1 <- readxl::read_excel(file1)
  df2 <- readxl::read_excel(file2)
  
  required_cols <- c(symbol_col, lfc_col, padj_col)
  missing_1 <- setdiff(required_cols, colnames(df1))
  missing_2 <- setdiff(required_cols, colnames(df2))
  if (length(missing_1)) stop("File 1 missing columns: ", paste(missing_1, collapse = ", "))
  if (length(missing_2)) stop("File 2 missing columns: ", paste(missing_2, collapse = ", "))
  
  # ---- 2) Optional sig gene tables ----
  sig1 <- df1 %>%
    mutate(
      .padj = as.numeric(.data[[padj_col]]),
      .lfc  = as.numeric(.data[[lfc_col]])
    ) %>%
    filter(.padj < padj_cutoff, abs(.lfc) > lfc_cutoff) %>%
    transmute(
      !!symbol_col := as.character(.data[[symbol_col]]),
      !!lfc_col := .lfc,
      !!padj_col := .padj
    )
  
  sig2 <- df2 %>%
    mutate(
      .padj = as.numeric(.data[[padj_col]]),
      .lfc  = as.numeric(.data[[lfc_col]])
    ) %>%
    filter(.padj < padj_cutoff, abs(.lfc) > lfc_cutoff) %>%
    transmute(
      !!symbol_col := as.character(.data[[symbol_col]]),
      !!lfc_col := .lfc,
      !!padj_col := .padj
    )
  
  # ---- 3) Rank vectors ----
  rank1 <- collapse_rank_by_log2fc(df1, symbol_col = symbol_col, lfc_col = lfc_col)
  rank2 <- collapse_rank_by_log2fc(df2, symbol_col = symbol_col, lfc_col = lfc_col)
  
  stopifnot(!anyDuplicated(names(rank1)))
  stopifnot(!anyDuplicated(names(rank2)))
  
  # ---- 4) MSigDB TERM2GENE ----
  msig <- msigdbr(species = species, category = category, subcategory = subcategory) %>%
    distinct(gs_name, gene_symbol)
  
  # ---- 5) GSEA ----
  gsea1 <- clusterProfiler::GSEA(
    geneList = rank1,
    TERM2GENE = msig,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = pvalueCutoff,
    verbose = FALSE
  )
  
  gsea2 <- clusterProfiler::GSEA(
    geneList = rank2,
    TERM2GENE = msig,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = pvalueCutoff,
    verbose = FALSE
  )
  
  res1 <- as.data.frame(gsea1)
  res2 <- as.data.frame(gsea2)
  
  # If either has no results, fail gracefully with a helpful message
  if (nrow(res1) == 0) stop("GSEA returned 0 results for Dataset 1. Check gene symbols/species or MSigDB settings.")
  if (nrow(res2) == 0) stop("GSEA returned 0 results for Dataset 2. Check gene symbols/species or MSigDB settings.")
  
  # ---- 6) Overlap + common UP/DOWN ----
  common <- inner_join(
    res1,
    res2,
    by = "ID",
    suffix = c(paste0("_", label1_safe), paste0("_", label2_safe))
  )
  
  # Build column names post-join
  p1   <- paste0("p.adjust_", label1_safe)
  p2   <- paste0("p.adjust_", label2_safe)
  nes1 <- paste0("NES_", label1_safe)
  nes2 <- paste0("NES_", label2_safe)
  
  # Significant in BOTH
  common_sig <- common %>%
    filter(.data[[p1]] < gsea_padj_cutoff, .data[[p2]] < gsea_padj_cutoff)
  
  # Common direction (within significant-in-both)
  common_up <- common_sig %>%
    filter(.data[[nes1]] > 0, .data[[nes2]] > 0) %>%
    arrange(desc(pmin(.data[[nes1]], .data[[nes2]])))
  
  common_down <- common_sig %>%
    filter(.data[[nes1]] < 0, .data[[nes2]] < 0) %>%
    arrange(pmax(.data[[nes1]], .data[[nes2]])) # more negative first
  
  # ---- 7) Plot (paired points) ----
  plot_source <- if (plot_set == "drivers") common_sig else common
  
  if (nrow(plot_source) == 0) {
    warning("No pathways to plot for plot_set = '", plot_set, "'. Plot will be NULL.")
    p_common <- NULL
    plot_df <- NULL
  } else {
    # Prefer set size from dataset 1 (exists as setSize_{label1} after join)
    setsize1 <- paste0("setSize_", label1_safe)
    if (!setsize1 %in% colnames(plot_source)) setsize1 <- NULL
    
    plot_df <- plot_source %>%
      transmute(
        ID,
        pathway = ID,
        setSize_1 = if (!is.null(setsize1)) .data[[setsize1]] else NA_real_,
        NES_1 = .data[[nes1]],
        NES_2 = .data[[nes2]]
      ) %>%
      pivot_longer(cols = c(NES_1, NES_2), names_to = "Dataset", values_to = "NES") %>%
      mutate(
        Dataset = recode(
          Dataset,
          "NES_1" = label1,
          "NES_2" = label2
        ),
        NES_dir = ifelse(NES > 0, "NES > 0", "NES < 0"),
        pathway = reorder(pathway, abs(NES), FUN = max)
      )
    
    p_common <- ggplot(plot_df, aes(x = NES, y = pathway)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_point(
        aes(size = setSize_1, color = NES_dir, shape = Dataset),
        alpha = 0.9
      ) +
      # Avoid hard-coded colors for accessibility/portability: use default palette
      scale_shape_manual(values = c(16, 1), name = "Dataset") +
      scale_size_continuous(name = "setSize (from Dataset 1)") +
      labs(
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        color = "Direction",
        title = paste0("Common pathways: ", label1, " vs ", label2),
        subtitle = paste0(
          "MSigDB ", category,
          if (!is.null(subcategory)) paste0(" (", subcategory, ")") else "",
          " • plot_set=", plot_set,
          if (plot_set == "drivers") paste0(" • both padj < ", gsea_padj_cutoff) else ""
        )
      ) +
      theme_classic() +
      theme(axis.text.y = element_text(size = 10),
            legend.position = "right")
  }
  
  # ---- 8) Exports ----
  if (export_plots && !is.null(p_common)) {
    ggsave(
      filename = file.path(outdir_final, paste0("p_", plot_set, "_", prefix_safe, ".png")),
      plot = p_common, width = 10, height = 8, dpi = 600
    )
    ggsave(
      filename = file.path(outdir_final, paste0("p_", plot_set, "_", prefix_safe, ".pdf")),
      plot = p_common, width = 10, height = 8
    )
  }
  
  if (export_tables) {
    readr::write_tsv(common,      file.path(outdir_final, paste0("common_paths_all_", prefix_safe, ".tsv")))
    readr::write_tsv(common_sig,  file.path(outdir_final, paste0("common_paths_sig_both_", prefix_safe, ".tsv")))
    readr::write_tsv(common_up,   file.path(outdir_final, paste0("common_up_", prefix_safe, ".tsv")))
    readr::write_tsv(common_down, file.path(outdir_final, paste0("common_down_", prefix_safe, ".tsv")))
    readr::write_tsv(sig1,        file.path(outdir_final, paste0("sig_genes_", label1_safe, "_", prefix_safe, ".tsv")))
    readr::write_tsv(sig2,        file.path(outdir_final, paste0("sig_genes_", label2_safe, "_", prefix_safe, ".tsv")))
  }
  
  if (save_rds) {
    saveRDS(
      list(
        meta = list(
          prefix = prefix_safe,
          outdir = outdir_final,
          file1 = file1, file2 = file2,
          label1 = label1, label2 = label2,
          symbol_col = symbol_col, lfc_col = lfc_col, padj_col = padj_col,
          padj_cutoff = padj_cutoff, lfc_cutoff = lfc_cutoff,
          gsea_padj_cutoff = gsea_padj_cutoff,
          species = species, category = category, subcategory = subcategory,
          minGSSize = minGSSize, maxGSSize = maxGSSize
        ),
        df1 = df1, df2 = df2,
        sig1 = sig1, sig2 = sig2,
        rank1 = rank1, rank2 = rank2,
        term2gene = msig,
        gsea1 = gsea1, gsea2 = gsea2,
        common = common,
        common_sig = common_sig,
        common_up = common_up,
        common_down = common_down,
        plot = p_common,
        plot_df = plot_df
      ),
      file = file.path(outdir_final, paste0("run_", prefix_safe, "_objects.rds"))
    )
  }
  
  message("Done. Outputs written to: ", outdir_final)
  
  invisible(list(
    outdir = outdir_final,
    gsea1 = gsea1,
    gsea2 = gsea2,
    common = common,
    common_sig = common_sig,
    common_up = common_up,
    common_down = common_down,
    plot = p_common
  ))
}

# ---- Example usage ----
# run_common_gsea(
#   file1 = "EO.xlsx",
#   file2 = "mAst.xlsx",
#   label1 = "EO771Br_RNA",
#   label2 = "mAst_RNA",
#   prefix = "Hallmark_RNA",
#   outdir = ".",
#   species = "Mus musculus",
#   category = "H",
#   plot_set = "drivers"
# )
