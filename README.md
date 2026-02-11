# common-gsea-pathway-comparison
Compare common up/downregulated pathways across datasets using GSEA

## Requirements

R ≥ 4.1 recommended

install.packages(c(
  "dplyr", "tidyr", "tibble", "readxl", "ggplot2", "readr"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "msigdbr"))

## Input format

Each dataset must be an Excel file (.xlsx) with at least:

###Column name	Description
symbol	Gene symbol (e.g., TP53, Fgf2)
log2FoldChange	Effect size (used for ranking)
padj	Adjusted p-value (optional filters)

### Example:

symbol   log2FoldChange   padj
Fgf2     -1.25            0.001
Cxcl10    2.10            0.0003

## Usage
source("common_gsea.R")

run_common_gsea(
  file1  = "EO771Br_RNA.xlsx",
  file2  = "mAst_RNA.xlsx",
  label1 = "EO771Br_RNA",
  label2 = "mAst_RNA",
  prefix = "Hallmark_RNA",
  outdir = ".",
  species = "Mus musculus",
  category = "H",        # Hallmarks
  plot_set = "drivers"  # "drivers" or "common"
)


## This will create a timestamped output folder containing:

common_paths_all_*.tsv

common_paths_sig_both_*.tsv

common_up_*.tsv

common_down_*.tsv

p_drivers_*.png and .pdf

run_*_objects.rds

## Key parameters
Parameter	Description
label1, label2	Names for datasets in plots
species	"Mus musculus" or "Homo sapiens"
category	"H" (Hallmark), "C2", "C5" (GO)
subcategory	Optional (e.g., "CP:REACTOME")
plot_set	"drivers" (sig in both) or "common"
gsea_padj_cutoff	FDR threshold for common drivers
padj_cutoff	Gene-level filter for sig tables
lfc_cutoff	Gene-level effect size filter

## Output interpretation

NES (Normalized Enrichment Score)
Positive NES = pathway enriched in upregulated genes
Negative NES = pathway enriched in downregulated genes

common_up.tsv
Pathways upregulated in both datasets

common_down.tsv
Pathways downregulated in both datasets

Dot plot
Each pathway appears twice (one point per dataset)
X-axis = NES, Y-axis = pathway
Shape = dataset, Color = direction

## Typical use cases

RNA-seq vs Ribo-seq pathway concordance

Treatment vs control across models

Tumor vs microenvironment comparison

Cross-species pathway agreement

Validation of perturbation signatures

## Reproducibility

Each run saves:

full inputs

GSEA objects

results tables

plotting data
into a single .rds file for downstream reuse.

## Citation

If you use this script in a publication or preprint, please cite:

**This repository:**
> **Márquez-Ortiz RA**. *Common GSEA Pathway Comparison Utility*. GitHub repository, 2026.

**Gene Set Enrichment Analysis (GSEA):**
> Subramanian A, Tamayo P, Mootha VK, et al.  
> *Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.*  
> Proc Natl Acad Sci U S A (2005). doi: 10.1073/pnas.0506580102 :contentReference[oaicite:0]{index=0}

**Enrichment implementation in R (clusterProfiler):**
> Wu T, Hu E, Xu S, et al.  
> *clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.*  
> The Innovation (2021). doi: 10.1016/j.xinn.2021.100141 :contentReference[oaicite:1]{index=1}

**Gene sets (MSigDB Hallmark collection):**
> Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P.  
> *The Molecular Signatures Database (MSigDB) hallmark gene set collection.*  
> Cell Systems (2015). doi: 10.1016/j.cels.2015.12.004 :contentReference[oaicite:2]{index=2}

*MSigDB access in R used here: `{msigdbr}` (CRAN package; provides MSigDB gene sets in tidy format).* :contentReference[oaicite:3]{index=3}
