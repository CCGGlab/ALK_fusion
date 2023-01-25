## Environment

Analysis was performed in a Conda environment. See **ALK_fusion.yml** for details. **scripts/Rpacks** describes R packages that were installed independently.

## Experimental conditions

### RNA-Seq data

NL20 cell lines (non-tumorigenic human bronchial epithelial cell line), 2 replicates

4 fusions + ctrl with 3 colonies for each fusion
* EML-ALK-V1: C1, D1, F12
* EML-ALK-V3: A3, B8, C2
* KIF5B-ALK: E2, E5, G1
* TFG-ALK: C7, D4, F11
* Ctrl: A9, B11, C3

Always without and with  (doxycycline) induction of fusion protein

Totalling to 60 samples (2x5x3x2)

### Phosphoproteomics data

4 fusions + control, either with or without (doxycycline) induction of fusion protein as well as treatment with lorlatinib, with differing numbers of replicates:
* Ctrl, TFG-ALK, KIF5B-ALK: 0 replicates
* EML4-ALK-V1, EML4-ALK-V3: 2 replicates

## Data
Raw RNA-Seq data have been deposited to Arrayexpress
* ALK-Fusion samples: E-MTAB-11304
* NSCLC cell lines with ALK inhibitors: E-MTAB-11342(CUTO* & YU1077),E-MTAB-11342 (H2228 & H3122)
* Raw phosphoproteomics data have been deposited to PRIDE: PXD035100

* Processed data are available in data/ALK_fusion_data.RData. This file containes the following objects:
  * sample_info[["RNA"]]: sample information cell line RNA-Seq data
  * normalized_counts[["RNA"]]: DESeq2-normalized counts from cell line RNA-Seq data
  * res_diff_expr[["RNA"]]: DESeq2 output from cell line RNA-Seq data
  * sample_info[["PP"]]: sample information cell line phospho-proteomics data
  * normalized_counts[["PP"]]: Normalized counts from cell line phosphoproteomics data
  * res_diff_expr[["PP"]]: DEP output from cell line phosphoproteomics data
  * sample_info[["Proteomics"]]: sample information cell line total proteomics data (not included in the final analysis and publication)
  * normalized_counts[["Proteomics"]]: Normalized counts from cell line total proteomics data (not included in the final analysis and publication)
  * res_diff_expr[["Proteomics"]]: DEP output from cell line phosphoproteomics data (not included in the final analysis and publication)
  * geneset_ls: genesets downloaded from MSigDB v7.2
  * NSCLCL_res: DESeq2 output from cell line RNA-Seq data 
  * H2228_res: DESeq2 output from cell line RNA-Seq data 
  * H3122_res: DESeq2 output from cell line RNA-Seq data 
  * SH2andPTBfromSuppExcel2: SH2 and PTB containing proteins downloaded from InterPro

## Data processing

### RNA-Seq

#### Processing of fastq files

```{bash}
scripts/other/process_fastq.sh
```

```{bash}
scripts/other/alignment_summary.sh
```

#### DESeq2 analysis

Structuring meta data

```{r}
source("scripts/get_sample_info.R")
```

Differential expression analysis and generation of normalized counts

```{r}
source("scripts/DESeq2_analysis.R")
```

### Phosphoproteomics

#### Preprocessing of raw phosphoproteomics data and structuring meta data

```{r}
source("scripts/PP_processing.R")
```

#### DEP analysis

Differential phosphorylation analysis

```{r}
source("scripts/DEP_PP.R")
```

#### Motif generation

Using UniProt to get +/-5 amino acids around identified phosphosites

```{r}
source("scripts/get_motif_data.R")
```

## Manuscript
The main analysis, as reported in the manuscript

### RNA-Seq 

Overall RNA response to fusions

```{r}
source("scripts/manuscript_RNA.R")
```

Inflammatory response analysis.

  ```{r}
source("scripts/manuscript_inflammatory_nfkb_response_enrichment.R")
```

### Phospho-proteomics

Overall phosphorylation response to fusions and treatment with lorlatinib

  ```{r}
source("scripts/manuscript_PP.R")
```

Heatmap of phosphosites of interest.

  ```{r}
source("scripts/manuscript_PP_heatmap.R")
```


