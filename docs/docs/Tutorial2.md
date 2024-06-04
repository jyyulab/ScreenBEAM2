---
layout: default
title: 2.Cleanning
parent: Guided tutorials
nav_order: 5

---

# Tutorial-2 Data annotation and data cleanning
{:.no_toc}
ScreenBEAM2 is a R based tool which consists of three major parts for processing steps: 1. mapping long read sequence to short read libraries, 2. Quality control, data cleanning and data preprocessing for mapped raw counts data; 3. Differential representative analysis on gene level or shRNA level.

**This is the second part of the whole tutorial, which is focused on expression set cleaning and meta data integration.**

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## Step2 Data annotation and data cleanning
### Step2.1 Data normalization and create Expressionset

We could perform data normalization and store different normalization method into eset. Total number of normalization is data dependent, usually we use 1e6. Save your expression matrix based on number of mismatch, this optimal mismatch number could be found in previous mapping QC report, as an example here we use 9. This function will also save an expression set with maximum mismatch and 0 mismatch, for your reference. Finnaly, the function will save normalized data and ids (gene/shRNA) in a tsv file for downstream analysis.

```R
analysis.par <- ScreenBEAM.createEset(analysis.par, normalize.total = 1e6, n.mismatch = 9)
```

### Step2.2 Qualtiy control for wrapped expressionsets
Here we have a build-in QC function imported from NetBID2 -- `draw.eset.QC` to facilitate quality control purpose. Each `draw.eset.QC` function will output a quality control report for following expression set.

```R
load(analysis.par$RData.eset.path$filepath[1])
draw.eset.QC(count.Nmm.normalized.eset, outdir = analysis.par$out.dir.output.QC, intgroup = 'group',
            do.logtransform = T, prefix = 'raw.9mm_', choose_plot = c("heatmap", "pca","density","correlation","meansd"))
```

Here we could also choose to perform more quality control for different numbers of mismatch, as a reference.

