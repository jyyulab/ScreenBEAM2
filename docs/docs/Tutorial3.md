---
layout: default
title: 3.DR
parent: Guided tutorials
nav_order: 5

---

# Tutorial-3 Differential representative analysis
{:.no_toc}
ScreenBEAM2 is a R based tool which consists of three major parts for processing steps: 1. mapping long read sequence to short read libraries, 2. Quality control, data cleanning and data preprocessing for mapped raw counts data; 3. Differential representative analysis on gene level or shRNA level.

**This is the third part of the whole tutorial, which is focused on differential representative analysis on gene and shRNA.**


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## Step3.Perform pairwise comparisons
First, define your comparison pairs and make sure all group names are correct. Here as sample analysis, we used the combination of all unique groups as our comparison. 
```R
case_group <- 'Hi';control_group <- 'Lo' 
```

Here we use function `ScreenBEAM.Pairwise` in ScreenBEAM2 to perform pairwise comparison

```
## Gene level
analysis.par <- ScreenBEAM.Pairwise(analysis.par,choose_level = 'gene', use_index = 1,data.type = 'NGS',
                case_group = case_group, control_group = control_group,
                do.normalization = FALSE)
## RNA level
analysis.par <- ScreenBEAM.Pairwise(analysis.par,choose_level = 'RNA', use_index = 1,
                case_group = case_group, control_group = control_group,
                do.normalization = FALSE, count.cutoff = 16, pooling = 'full')

```

If you have NetBID2 installed, you could use built-in volcano plot function to generate volcano plot for you DR analysis.

```R
# Take 1 comparison as example, to return the list of significant genes and draw a plot
sig_gene <- draw.volcanoPlot(dat=DR.GENE.DF.sel, label_col = "geneID", logFC_col = names(DR.GENE.DF.sel)[2], Pv_col = names(DR.GENE.DF.sel)[5], logFC_thre = 2, Pv_thre = 1e-4, main = names(de.rna.list), show_label = T, label_cex = 1, pdf_file =  paste0(analysis.par$out.dir.output.DR, names(de.gene.list),"_GENE.pdf"))
```

