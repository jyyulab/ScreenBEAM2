---
layout: default
title: Installation
nav_order: 2

---

# Installation

## Install blat

You can install blat from conda:

```shell
conda install -c bioconda blat
which blat
```

## Install R package from github

```R
devtools::install(pkg='.', dependencies = T)
devtools::install_deps(pkg = ".", dependencies = TRUE) ## Install package dependencies if needed.
```