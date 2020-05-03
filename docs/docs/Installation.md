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
If you are running ScreenBEAM on local node, please make sure your blat has been linked to your local bin.

## Install ScreenBEAM2 package 

### Install R package from source

In your R session, please run

```R
install.packages(path_to_file, type="source",repo=NULL)
```
or 

```R
devtools::install_local(path_to_file)
```
### R package from github

```R
devtools::install(pkg='.', dependencies = T)
devtools::install_deps(pkg = ".", dependencies = TRUE) ## Install package dependencies if needed.
```