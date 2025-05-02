DESeq2
================
Marco Tello
2025-05-01

# Differential expression analysis

This report is to identify differentially expressed genes across the
following conditions:

1.  low glucose control (LG)
2.  high glucose control (HG)
3.  high glucose high fat control (HF)
4.  high glucose high fat PTS (PTS)
5.  high glucose high fat QUE (QUE)

## Format condition names

The first step is to load the selected expression data and format the
condition names to match the experimental conditions.

``` r
# Read exp table
input_file <- file.path(data_path, "protein_coding_exp_filtered-S10-R6.tsv")
exp_df <- fread(input = input_file, drop = c("P2_S2", "Q3_S1"))


# Connect to the Ensembl mart
ensembl <- useEnsembl(
  biomart    = "genes",
  dataset    = "hsapiens_gene_ensembl")

# Query for gene symbol
bm <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = exp_df$gene_id,
    mart       = ensembl)
exp_df <- merge.data.table(x = bm, y = exp_df, 
                           by.x = "ensembl_gene_id", by.y = "gene_id", 
                           all.x = FALSE, all.y = TRUE)
exp_df <- data.table(exp_df)

# Remove genes with no gene symbol
exp_df <- exp_df[external_gene_name != ""]

# Only one gene symbol has two ENSEMBL IDs so we handle it separately
exp_df[external_gene_name == "POLR2J3", external_gene_name := paste(external_gene_name, 
                                                                    ensembl_gene_id, 
                                                                    sep = "_")]

# Convert data.table to integer matrix
exp_matrix <- exp_df %>%
  dplyr::select(!ensembl_gene_id) %>%
  dplyr::select(!external_gene_name) %>%
  as.matrix()
rownames(exp_matrix) <- exp_df$external_gene_name


# Capture metadata in data frame
condition <- colnames(exp_matrix) %>% 
  str_replace(pattern = "\\d_S\\d+", 
              replacement = "")

# Format condition to match labels
condition <- case_when(
    condition == "C" ~ "HG",
    condition == "F" ~ "HF",
    condition == "L" ~ "LG",
    condition == "P" ~ "PTS",
    condition == "Q" ~ "QUE")
# Create variables for DE analysis
annot <- data.frame(condition = factor(condition, levels = c("LG", "HG", "HF", "PTS", "QUE")),
                    Hglucose  = ifelse(condition == "LG", 0, 1),
                    Hfat      = ifelse(condition %in% c("LG", "HG"), 0, 1),
                    PTS       = ifelse(condition == "PTS", 1, 0),
                    QUE       = ifelse(condition == "QUE", 1, 0),
                    row.names = colnames(exp_matrix))

annot 
```

    ##        condition Hglucose Hfat PTS QUE
    ## C1_S44        HG        1    0   0   0
    ## C2_S45        HG        1    0   0   0
    ## C3_S46        HG        1    0   0   0
    ## F1_S47        HF        1    1   0   0
    ## F2_S48        HF        1    1   0   0
    ## F3_S49        HF        1    1   0   0
    ## F4_S50        HF        1    1   0   0
    ## L1_S57        LG        0    0   0   0
    ## L2_S58        LG        0    0   0   0
    ## L3_S59        LG        0    0   0   0
    ## P1_S54       PTS        1    1   1   0
    ## P2_S55       PTS        1    1   1   0
    ## P3_S56       PTS        1    1   1   0
    ## Q1_S51       QUE        1    1   0   1
    ## Q2_S52       QUE        1    1   0   1
    ## Q3_S53       QUE        1    1   0   1

Based on the experimental groups we created the appropriate labels and
encoded conditions as binary columns to facilitate building models.

## Set-up DEA parameters

For consistency, we’ll define the parameters for differential expression
here:

``` r
alpha <- 0.01
log2FC <- 1
```

The comparisons we’ll be looking into are the following:

1.  LG vs HG - to see how increasing glucose levels affects gene
    expression without fat
2.  HG vs HF - to isolate the impact of the high fat component under
    hyperglycemic conditions
3.  HF vs PTS - whether PTS can reverse or alter HF-induced gene changes
4.  HF vs QUE - whether QUE can reverse or alter HF-induced gene changes
5.  PTS vs QUE - to identify differences in molecular responses between
    the two compounds.
6.  LG vs PTS - to assess whether PTS brings HF-exposed cells back
    toward the low-glucose baseline
7.  LG vs QUE - to assess whether QUE brings HF-exposed cells back
    toward the low-glucose baseline

Thus, we’ll create a section for each of the target comparisons
containing:

- Model fitting and result selection
- p-value histograms
- volcano plot
- Output of DEGs

\##1. LG vs HG

how increasing glucose levels affects gene expression without fat?

``` r
# Re-define base condition
annot$condition <- relevel(annot$condition, ref = "LG")
# Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = exp_matrix,
                              colData = annot,
                              design= ~ condition)

# Run DESeq2
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
# Print comparisons made
resultsNames(dds)
```

    ## [1] "Intercept"           "condition_HG_vs_LG"  "condition_HF_vs_LG" 
    ## [4] "condition_PTS_vs_LG" "condition_QUE_vs_LG"

``` r
comparison <- "condition_HG_vs_LG"

de_res <- results(dds,
                  name=comparison,
                  pAdjustMethod = "BH", 
                  alpha = alpha) 
temp <- as.data.table(de_res)
de_res <- data.table(geneID = row.names(de_res), temp)


out_file <- file.path(out_path, paste("deseq2_", 
                                      str_replace_all(comparison, "_", ""),
                                      "_FDR", str_replace(alpha, "\\.", ""),
                                      ".tsv", sep = ""))
if(!file.exists(out_file)){
  output_table = de_res
  output_table[, log2FoldChange := round(log2FoldChange, 4)]
  output_table[, baseMean       := round(baseMean, 4)]
  output_table[, lfcSE          := round(lfcSE, 4)]
  output_table[, stat           := round(stat, 4)]
  output_table[, pvalue         := round(pvalue, 4)]
  output_table[, padj           := round(padj, 4)]
  output_table[, FoldChange := round((2**abs(log2FoldChange))*sign(log2FoldChange), 4)]
  fwrite(x = output_table, file = out_file, 
         append = FALSE, quote = FALSE, sep = '\t', 
         row.names = FALSE, col.names = TRUE)
}

temp <- visualize_degs(de_res, alpha, log2FC)
temp
```

    ## [[1]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

    ## 
    ## [[2]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

\##2. HG vs HF

``` r
# Re-define base condition
annot$condition <- relevel(annot$condition, ref = "HG")
# Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = exp_matrix,
                              colData = annot,
                              design= ~ condition)

# Run DESeq2
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
# Print comparisons made
resultsNames(dds)
```

    ## [1] "Intercept"           "condition_LG_vs_HG"  "condition_HF_vs_HG" 
    ## [4] "condition_PTS_vs_HG" "condition_QUE_vs_HG"

``` r
comparison <- "condition_HF_vs_HG"

de_res <- results(dds,
                  name=comparison,
                  pAdjustMethod = "BH", 
                  alpha = alpha) 
temp <- as.data.table(de_res)
de_res <- data.table(geneID = row.names(de_res), temp)


out_file <- file.path(out_path, paste("deseq2_", 
                                      str_replace_all(comparison, "_", ""),
                                      "_FDR", str_replace(alpha, "\\.", ""),
                                      ".tsv", sep = ""))
if(!file.exists(out_file)){
  output_table = de_res
  output_table[, log2FoldChange := round(log2FoldChange, 4)]
  output_table[, baseMean       := round(baseMean, 4)]
  output_table[, lfcSE          := round(lfcSE, 4)]
  output_table[, stat           := round(stat, 4)]
  output_table[, pvalue         := round(pvalue, 4)]
  output_table[, padj           := round(padj, 4)]
  output_table[, FoldChange := round((2**abs(log2FoldChange))*sign(log2FoldChange), 4)]
  fwrite(x = output_table, file = out_file, 
         append = FALSE, quote = FALSE, sep = '\t', 
         row.names = FALSE, col.names = TRUE)
}

temp <- visualize_degs(de_res, alpha, log2FC)
temp
```

    ## [[1]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

    ## 
    ## [[2]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

\##3. HF vs PTS

``` r
# Re-define base condition
annot$condition <- relevel(annot$condition, ref = "HF")
# Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = exp_matrix,
                              colData = annot,
                              design= ~ condition)

# Run DESeq2
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
# Print comparisons made
resultsNames(dds)
```

    ## [1] "Intercept"           "condition_HG_vs_HF"  "condition_LG_vs_HF" 
    ## [4] "condition_PTS_vs_HF" "condition_QUE_vs_HF"

``` r
comparison <- "condition_PTS_vs_HF"

de_res <- results(dds,
                  name=comparison,
                  pAdjustMethod = "BH", 
                  alpha = alpha) 
temp <- as.data.table(de_res)
de_res <- data.table(geneID = row.names(de_res), temp)


out_file <- file.path(out_path, paste("deseq2_", 
                                      str_replace_all(comparison, "_", ""),
                                      "_FDR", str_replace(alpha, "\\.", ""),
                                      ".tsv", sep = ""))
if(!file.exists(out_file)){
  output_table = de_res
  output_table[, log2FoldChange := round(log2FoldChange, 4)]
  output_table[, baseMean       := round(baseMean, 4)]
  output_table[, lfcSE          := round(lfcSE, 4)]
  output_table[, stat           := round(stat, 4)]
  output_table[, pvalue         := round(pvalue, 4)]
  output_table[, padj           := round(padj, 4)]
  output_table[, FoldChange := round((2**abs(log2FoldChange))*sign(log2FoldChange), 4)]
  fwrite(x = output_table, file = out_file, 
         append = FALSE, quote = FALSE, sep = '\t', 
         row.names = FALSE, col.names = TRUE)
}

temp <- visualize_degs(de_res, alpha, log2FC)
temp
```

    ## [[1]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    ## 
    ## [[2]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

\##4. HF vs QUE

``` r
# Re-define base condition
annot$condition <- relevel(annot$condition, ref = "HF")
# Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = exp_matrix,
                              colData = annot,
                              design= ~ condition)

# Run DESeq2
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
# Print comparisons made
resultsNames(dds)
```

    ## [1] "Intercept"           "condition_HG_vs_HF"  "condition_LG_vs_HF" 
    ## [4] "condition_PTS_vs_HF" "condition_QUE_vs_HF"

``` r
comparison <- "condition_QUE_vs_HF"

de_res <- results(dds,
                  name=comparison,
                  pAdjustMethod = "BH", 
                  alpha = alpha) 
temp <- as.data.table(de_res)
de_res <- data.table(geneID = row.names(de_res), temp)


out_file <- file.path(out_path, paste("deseq2_", 
                                      str_replace_all(comparison, "_", ""),
                                      "_FDR", str_replace(alpha, "\\.", ""),
                                      ".tsv", sep = ""))
if(!file.exists(out_file)){
  output_table = de_res
  output_table[, log2FoldChange := round(log2FoldChange, 4)]
  output_table[, baseMean       := round(baseMean, 4)]
  output_table[, lfcSE          := round(lfcSE, 4)]
  output_table[, stat           := round(stat, 4)]
  output_table[, pvalue         := round(pvalue, 4)]
  output_table[, padj           := round(padj, 4)]
  output_table[, FoldChange := round((2**abs(log2FoldChange))*sign(log2FoldChange), 4)]
  fwrite(x = output_table, file = out_file, 
         append = FALSE, quote = FALSE, sep = '\t', 
         row.names = FALSE, col.names = TRUE)
}

temp <- visualize_degs(de_res, alpha, log2FC)
temp
```

    ## [[1]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

    ## 
    ## [[2]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

\##5. PTS vs QUE

``` r
# Re-define base condition
annot$condition <- relevel(annot$condition, ref = "PTS")
# Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = exp_matrix,
                              colData = annot,
                              design= ~ condition)

# Run DESeq2
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
# Print comparisons made
resultsNames(dds)
```

    ## [1] "Intercept"            "condition_HF_vs_PTS"  "condition_HG_vs_PTS" 
    ## [4] "condition_LG_vs_PTS"  "condition_QUE_vs_PTS"

``` r
de_res <- results(dds,
                  name="condition_QUE_vs_PTS",
                  pAdjustMethod = "BH", 
                  alpha = alpha) 
summary(de_res)
```

    ## 
    ## out of 13281 with nonzero total read count
    ## adjusted p-value < 0.01
    ## LFC > 0 (up)       : 1730, 13%
    ## LFC < 0 (down)     : 2684, 20%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 515, 3.9%
    ## (mean count < 12)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
temp <- as.data.table(de_res)
temp <- data.table(geneID = row.names(de_res), temp)
head(temp[padj <= alpha, ])
```

    ##     geneID  baseMean log2FoldChange      lfcSE      stat       pvalue
    ##     <char>     <num>          <num>      <num>     <num>        <num>
    ## 1:    DPM1 1360.9142     -0.2529432 0.05768660 -4.384783 1.161016e-05
    ## 2:   FUCA2 4589.7983     -0.2742815 0.04125285 -6.648789 2.955148e-11
    ## 3:    GCLC 3256.8692      0.1591416 0.03807268  4.179942 2.915839e-05
    ## 4:    NFYA 1508.5609     -0.4520951 0.05749965 -7.862571 3.763281e-15
    ## 5: CYP51A1  306.7181     -0.5148257 0.13320897 -3.864797 1.111814e-04
    ## 6:    LAP3 3429.7680     -0.3003817 0.04629009 -6.489113 8.634300e-11
    ##            padj
    ##           <num>
    ## 1: 6.054548e-05
    ## 2: 3.724128e-10
    ## 3: 1.409451e-04
    ## 4: 7.577610e-14
    ## 5: 4.766091e-04
    ## 6: 1.025353e-09

``` r
temp <- visualize_degs(as.data.table(de_res), alpha, log2FC)
temp
```

    ## [[1]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

    ## 
    ## [[2]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

\##6. LG vs PTS

``` r
# Re-define base condition
annot$condition <- relevel(annot$condition, ref = "LG")
# Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = exp_matrix,
                              colData = annot,
                              design= ~ condition)

# Run DESeq2
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
# Print comparisons made
resultsNames(dds)
```

    ## [1] "Intercept"           "condition_PTS_vs_LG" "condition_HF_vs_LG" 
    ## [4] "condition_HG_vs_LG"  "condition_QUE_vs_LG"

``` r
de_res <- results(dds,
                  name="condition_PTS_vs_LG",
                  pAdjustMethod = "BH", 
                  alpha = alpha) 
summary(de_res)
```

    ## 
    ## out of 13281 with nonzero total read count
    ## adjusted p-value < 0.01
    ## LFC > 0 (up)       : 3684, 28%
    ## LFC < 0 (down)     : 3450, 26%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
temp <- as.data.table(de_res)
temp <- data.table(geneID = row.names(de_res), temp)
head(temp[padj <= alpha, ])
```

    ##    geneID  baseMean log2FoldChange      lfcSE      stat       pvalue
    ##    <char>     <num>          <num>      <num>     <num>        <num>
    ## 1: TSPAN6 3210.3993      0.3215855 0.04722124  6.810187 9.747211e-12
    ## 2:   DPM1 1360.9142      0.1612024 0.05672144  2.842002 4.483126e-03
    ## 3:  SCYL3  274.8483      0.4645241 0.10961638  4.237726 2.257953e-05
    ## 4:  FIRRM  329.4765      1.5800806 0.11339021 13.934894 3.887755e-44
    ## 5: NIPAL3  578.2322     -0.6593507 0.07323424 -9.003312 2.190086e-19
    ## 6:  LAS1L  894.3796      0.3249288 0.06676452  4.866788 1.134264e-06
    ##            padj
    ##           <num>
    ## 1: 4.075967e-11
    ## 2: 8.450241e-03
    ## 3: 5.713064e-05
    ## 4: 5.417973e-43
    ## 5: 1.370067e-18
    ## 6: 3.247285e-06

``` r
temp <- visualize_degs(as.data.table(de_res), alpha, log2FC)
temp
```

    ## [[1]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

    ## 
    ## [[2]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

\##7. LG vs QUE

``` r
# Re-define base condition
annot$condition <- relevel(annot$condition, ref = "LG")
# Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = exp_matrix,
                              colData = annot,
                              design= ~ condition)

# Run DESeq2
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
# Print comparisons made
resultsNames(dds)
```

    ## [1] "Intercept"           "condition_PTS_vs_LG" "condition_HF_vs_LG" 
    ## [4] "condition_HG_vs_LG"  "condition_QUE_vs_LG"

``` r
de_res <- results(dds,
                  name="condition_QUE_vs_LG",
                  pAdjustMethod = "BH", 
                  alpha = alpha) 
summary(de_res)
```

    ## 
    ## out of 13281 with nonzero total read count
    ## adjusted p-value < 0.01
    ## LFC > 0 (up)       : 3645, 27%
    ## LFC < 0 (down)     : 3822, 29%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
temp <- as.data.table(de_res)
temp <- data.table(geneID = row.names(de_res), temp)
head(temp[padj <= alpha, ])
```

    ##    geneID  baseMean log2FoldChange      lfcSE       stat       pvalue
    ##    <char>     <num>          <num>      <num>      <num>        <num>
    ## 1: TSPAN6 3210.3993      0.4003866 0.04676212   8.562200 1.107335e-17
    ## 2:  SCYL3  274.8483      0.4599131 0.10824956   4.248637 2.150753e-05
    ## 3:  FIRRM  329.4765      1.4166342 0.11333911  12.499076 7.552400e-36
    ## 4:  FUCA2 4589.7983     -0.1636830 0.04092797  -3.999295 6.353158e-05
    ## 5:   NFYA 1508.5609     -0.5773876 0.05607340 -10.296997 7.269539e-25
    ## 6: NIPAL3  578.2322     -0.5740807 0.07146889  -8.032595 9.543203e-16
    ##            padj
    ##           <num>
    ## 1: 5.967648e-17
    ## 2: 5.191398e-05
    ## 3: 7.570070e-35
    ## 4: 1.464612e-04
    ## 5: 5.204677e-24
    ## 6: 4.692457e-15

``` r
temp <- visualize_degs(as.data.table(de_res), alpha, log2FC)
temp
```

    ## [[1]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

    ## 
    ## [[2]]

![](ribodepleted_DEGs_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->
