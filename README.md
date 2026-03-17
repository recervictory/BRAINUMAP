# BRAINUMAP Working Manual

This repository contains the analysis code used to build and visualize a cross-cohort RNA-seq reference map of human normal and neoplastic brain tissue (BRAIN-UMAP), and to reproduce the figure panels from the associated study.

[![DOI](https://zenodo.org/badge/584982012.svg)](https://zenodo.org/badge/latestdoi/584982012)

[Scientific Reports publication](https://www.nature.com/articles/s41598-023-31180-z)

## 1. Repository Objective

The codebase has three core goals:

1. Build a harmonized expression matrix across adult and pediatric brain cohorts.
2. Generate low-dimensional embeddings (UMAP, PCA, tSNE) for biological interpretation.
3. Reproduce all main and supplemental visualizations used in the manuscript.

In practical terms, this repository turns cohort-specific RNA-seq inputs into:

1. Integrated expression matrices.
2. Metadata-enriched embedding tables.
3. Figure-ready PDF outputs.

## 2. Scientific Method Objective

The analysis strategy is:

1. Download and process public data (TCGA, GTEx, CGGA, CBTN).
2. Convert counts to comparable log2(TPM+1)-scale matrices.
3. Align genes across cohorts.
4. Optionally remove batch effects using ComBat.
5. Compute dimensionality reductions (UMAP 2D/3D, PCA, tSNE).
6. Overlay clinical, molecular, pathway, mutation, fusion, and copy-number annotations.
7. Export publication figures.

## 3. Repository Structure

Top-level files:

1. `process_data.Rmd` - Data acquisition and preprocessing for public datasets and CBTN.
2. `combine_datasets_umap.Rmd` - Cohort integration, batch correction, embedding computation, metadata merge.
3. `Fig1.Rmd` - Figure 1 and Supplemental Figure 1 generation.
4. `Fig2.Rmd` - Figure 2 and Supplemental Figure 2 generation.
5. `Fig3.Rmd` - Figure 3 and supplemental nearest-neighbor survival overlays.
6. `Fig4.Rmd` - Figure 4 and supplemental pediatric subtype visualizations.
7. `Fig5.Rmd` - Figure 5 and supplemental GSVA pathway overlays.
8. `Fig6.Rmd` - Figure 6 and supplemental pathway gene expression overlays.
9. `Fig7.Rmd` - Figure 7 and supplemental mutation/fusion burden boxplots.
10. `Fig8.Rmd` - Figure 8 and supplemental gene-level multi-omics overlays.

## 4. File-by-File Working Guide

### 4.1 `process_data.Rmd`

Purpose:

1. Download recount/recount-brain resources.
2. Build log2(TPM+1) matrices for TCGA/GTEx/CGGA.
3. Prepare CBTN log2(TPM+1) matrix.
4. Build a combined SummarizedExperiment with glioma subtype annotations.

Key outputs referenced downstream:

1. `data/log2_tpm_big_rse_1_28_2021.Rdata`
2. `data/FULL_pca_tsne_umap_combined_analysis_1_28_2021.txt`
3. `data/log2_tpm_cbtn.rds`

Method objective:

Create harmonized expression and metadata objects with consistent gene universe and sample annotations before dimensionality reduction.

### 4.2 `combine_datasets_umap.Rmd`

Purpose:

1. Combine adult brain cohorts with CBTN.
2. Apply ComBat batch correction and keep uncorrected comparator.
3. Generate UMAP (2D/3D), PCA, and tSNE tables.
4. Create subset embeddings (GTEx+CBTN, TCGA+CBTN).
5. Add metadata to embedding tables.

Key outputs:

1. `data/log2_tpm_combatseq_cbtn_brain_umap_8_27_2021.rds`
2. `tables/umapdata_with_CBTN_8_27_2021.txt`
3. `tables/no_batch_correction_umapdata_with_CBTN_2_1_2022.txt`
4. `tables/umapdata_TCGA_with_CBTN_8_27_2021.txt`
5. `tables/umapdata_GTEX_with_CBTN_8_27_2021.txt`
6. Metadata-augmented table variants (`*_with_metadata.txt`)

Method objective:

Create the core BRAIN-UMAP coordinate system used by all figure scripts.

### 4.3 `Fig1.Rmd`

Purpose:

1. Global cohort view of BRAIN-UMAP.
2. GTEx-focused tissue overlays.
3. With-vs-without batch correction comparison.

Method objective:

Demonstrate global cohort structure and the impact of harmonization.

### 4.4 `Fig2.Rmd`

Purpose:

1. TCGA glioma-focused overlays.
2. Subtype, age, copy-number, methylation, and mutation overlays.
3. Supplemental clinical and molecular context panels.

Method objective:

Resolve adult glioma axes and molecular stratification on BRAIN-UMAP.

### 4.5 `Fig3.Rmd`

Purpose:

1. Joint TCGA-CGGA glioma molecular labels (IDH, 1p/19q, grade).
2. Survival-time overlays in UMAP space.
3. Supplemental nearest-neighbor sensitivity views.

Method objective:

Connect map position with molecular subtype and clinical outcome gradients.

### 4.6 `Fig4.Rmd`

Purpose:

1. Pediatric CBTN disease map.
2. Full BRAIN-UMAP disease legend integration.
3. Supplemental ependymoma and medulloblastoma subtype overlays.

Method objective:

Show pediatric tumor diversity and how it situates relative to adult/normal cohorts.

### 4.7 `Fig5.Rmd`

Purpose:

1. GSVA pathway scoring across cohorts.
2. Differential pathway analysis (GTEx versus glioma classes).
3. Venn-style overlap summaries and supplemental pathway panels.

Method objective:

Translate map regions into biological pathway programs.

### 4.8 `Fig6.Rmd`

Purpose:

1. Gene-level expression overlays for pathway members.
2. Main mismatch-repair set and supplemental RELA pathway set.

Method objective:

Show within-pathway gene expression heterogeneity in map space.

### 4.9 `Fig7.Rmd`

Purpose:

1. Mutation burden overlays.
2. Fusion burden overlays.
3. Copy-number gain/loss burden overlays.
4. Supplemental burden boxplots by tumor type.

Method objective:

Quantify and visualize genomic event burden across BRAIN-UMAP regions.

### 4.10 `Fig8.Rmd`

Purpose:

For selected genes, build multi-panel overlays of:

1. Expression
2. Mutation status
3. Copy-number status
4. Fusion status

Method objective:

Provide gene-centric multi-omics context tied to location on the map.

## 5. Function Dictionary (Core R Functions and Roles)

This dictionary focuses on high-impact functions repeatedly used across scripts.

### 5.1 Data I/O and containers

1. `read.delim`, `read.csv`, `read_xlsx`, `readRDS` - Load metadata, tabular outputs, and serialized matrices.
2. `write.table`, `write_xlsx`, `saveRDS`, `save` - Persist analysis artifacts for downstream scripts.
3. `SummarizedExperiment`, `assay`, `colData`, `rowRanges` - Structured storage and retrieval of assay matrix plus metadata.

### 5.2 Data wrangling and harmonization

1. `match`, `intersect`, `identical` - Enforce sample/gene alignment and integrity checks.
2. `factor`, `levels`, `table` - Controlled categorical ordering and QC summaries.
3. `split`, `lapply`, `sapply`, `do.call` - Iterative transformations and cohort-wise summaries.
4. `gsub`, `grep`, `which`, `order` - Label normalization and plotting-order control.

### 5.3 Normalization and adjustment

1. `getRPKM`, `rpkm` - Expression normalization helpers.
2. `ComBat` - Batch correction across cohorts.

### 5.4 Dimensionality reduction

1. `umap` - UMAP embeddings (2D and 3D).
2. `prcomp` - PCA projection and explained variance.
3. `Rtsne` - tSNE projection.

### 5.5 Pathway analysis

1. `gsva` - Per-sample pathway activity scoring.
2. `model.matrix`, `lmFit`, `eBayes`, `topTable` - Linear modeling for pathway-level differential analysis.
3. `makeContrasts` - Contrast definitions for modeled group comparisons.

### 5.6 Visualization

1. `ggplot`, `geom_point`, `scale_colour_manual`, `scale_colour_gradientn`, `scale_color_distiller` - Main plotting grammar.
2. `theme_void`, `theme_classic`, custom `theme(...)` - Shared visual style across figures.
3. `grid.arrange`, `arrangeGrob`, `marrangeGrob` - Multi-panel layout and pagination.
4. `pdf`, `dev.off` - Figure export lifecycle.

### 5.7 Specialized utilities

1. `new_scale_color` (ggnewscale) - Multiple color scales in one panel.
2. `venn`/Euler set plotting utilities - Pathway overlap visualization.

## 6. Typical Execution Order

Recommended run order for full reproduction:

1. `process_data.Rmd`
2. `combine_datasets_umap.Rmd`
3. `Fig1.Rmd`
4. `Fig2.Rmd`
5. `Fig3.Rmd`
6. `Fig4.Rmd`
7. `Fig5.Rmd`
8. `Fig6.Rmd`
9. `Fig7.Rmd`
10. `Fig8.Rmd`

Rationale:

1. Figure scripts depend on embedding tables and metadata produced in the first two notebooks.
2. Many figure scripts also assume specific absolute/local file paths and previously generated files.

## 7. Inputs and Outputs at a Glance

Primary input families:

1. Public cohort expression and metadata (recount, GTEx, TCGA, CGGA).
2. CBTN RNA-seq counts and pediatric metadata.
3. Optional molecular tables (mutations, fusions, CNV, clinical annotations, survival metadata).

Primary output families:

1. Harmonized assay objects in `data/`.
2. Embedding tables in `tables/`.
3. Publication and supplemental figures in `figures/`.

## 8. Environment and Dependency Notes

Core R packages used across scripts include:

1. `SummarizedExperiment`
2. `recount`
3. `rtracklayer`
4. `edgeR`
5. `umap`
6. `Rtsne`
7. `sva`
8. `ggplot2`
9. `grid`, `gridExtra`
10. `RColorBrewer`
11. `GSVA`, `GSEABase`, `GSVAdata`
12. `limma`
13. `ggnewscale`
14. `readxl`, `writexl`

Important:

1. Several scripts currently reference user-specific absolute paths (for example under `~/HollandLabShared/...`).
2. For reproducible execution in a new environment, convert these paths into a configuration block (for example `resdir`, `maindir`, `data_dir`) at the top of each Rmd.

## 9. Practical Usage Guidance

Before running:

1. Confirm all required external source files are available.
2. Standardize directory roots (`resdir`, `maindir`) for your system.
3. Create expected output subfolders (`data`, `tables`, `figures`) if missing.

While running:

1. Run each Rmd top-to-bottom in a fresh R session.
2. Preserve object naming as used in scripts to avoid downstream breaks.
3. Keep the same sample ordering checks (`identical(...)`) to prevent silent mismatch errors.

After running:

1. Verify generated `*.txt`, `*.rds`, and `*.pdf` filenames match those expected by later scripts.
2. Check the number of rows/samples in intermediate tables before moving forward.

## 10. Known Reproducibility Risks

1. Hard-coded absolute paths can fail on new machines.
2. Scripts assume availability of external files not stored in this repository.
3. Naming mismatches (for example `GTEX` vs `GTEx`) can alter legend behavior and joins.
4. Some analyses rely on exact sample ID conventions across cohorts.

## 11. Suggested Refactor Roadmap (Optional)

To improve long-term maintainability, consider:

1. Centralizing paths and constants in one config script.
2. Moving repeated theme/color definitions into helper functions.
3. Encapsulating repeated overlay logic into reusable plotting functions.
4. Adding a lockfile-driven package environment (`renv`).
5. Adding a minimal validation script that checks for all required input files before execution.

## 12. Citation

If you use this workflow, please cite:

1. The Scientific Reports publication above.
2. The Zenodo DOI for this repository.
