---
layout: default
title: Running on cluster
nav_order: 6

---

# Tutorial-4 Running DR analysis for large data
{:.no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}


## Perform pairwise comparisons on Gene level by jobs in parallel
For whole genome libraries, it may cost sometime to run these pairwise comparisons sequentially. So we recommend to generate one R script for each comparisons and running them in parallel.

Define dependencies, function path (if package is not installed) , and input files.
```R
library(ScreenBEAM2)

input.file<-"your_normalized_eset.tsv"
input.file<-normalizePath(input.file)
dir.tmp<-"DR_hpc/" # directory to store individual results
```

```R
load([your_stored_eset])

# define your comparisons
d<-openxlsx::read.xlsx("your_CRISPR_comparisons.xlsx")
compare.pairs<-t(d) 

# turn off filtering 
for(i in 1:dim(compare.pairs)[2]){
  
  case.groupname<-compare.pairs[1,i]
  control.groupname<-compare.pairs[2,i]
  
  compare.name<-paste0(case.groupname,".vs.",control.groupname)

  case.sample.id<-which(eset$group==compare.pairs[1,i])
  control.sample.id<-which(eset$group==compare.pairs[2,i])
  
  control.samples<-eset$group[control.sample.id]
  case.samples<-eset$group[case.sample.id]
  
  inputs<-paste0("control.samples ='",control.samples ,"',case.samples = '", case.samples, 
              "',control.groupname='", control.groupname,"',case.groupname='",case.groupname,"'")

  other.args<-"data.type = 'NGS',gene.columnId=2, do.normalization=FALSE, filterLowCount=FALSE,
              count.cutoff=0, rna.size = 6, sample.rna.time = 100, pooling = 'partial', method='Bayesian'"
  
  fileR<-paste0('DR_hpc/',compare.name,".R")
  sink(fileR)
  cat('rm(list=ls())\n',sep='')
  cat('setwd(\'',normalizePath(dir.tmp),'\')\n',sep='')
  cat("source('dependencies.R')\n",sep='')
  cat('source(\'',func.path,'\')\n',sep='')
  cat('source(\'',func.util.path,'\')\n',sep='')
  cat('de.gene.list<-list()\n',sep='')
  cat('print(\'',compare.name,'\')\n',sep='')
  cat('de<-ScreenBEAM.gene.level(\'',input.file,'\',',inputs, ',', other.args,')\n',sep='')
  cat('names(de)[1]<-"geneID"\n',sep='')
  cat('de.gene.list[[\'',compare.name,'\']]<-de', '\n',sep='')
  cat('save(de.gene.list, file=\'DR.gene.c',i,'\')\n',sep = '')
  sink()
  
  compare.name<-gsub(".","_",compare.name,fixed = TRUE)
  
  file.sh<-file.path(dir.tmp,paste(compare.name,'sh',sep='.'))
  project.name<-compare.name
  memory<-paste('#BSUB -R \"rusage[mem=8000]\"')
  project<-paste('#BSUB -P ',project.name,'')
  
  sh.header<-paste(
    '#!/bin/bash',
    project,
    memory,
    '#BSUB -q standard',
    sep='\n')
  
  sink(file.sh)
  cat(sh.header,'\n', sep = '')
  cat('Rscript', normalizePath(fileR), sep=' ')
  sink()
  
  cmd<-paste('bsub < ',normalizePath(file.sh),'\n',sep='')
  system(cmd)
  system(paste('sleep ',1,'\n'))
  cat(cmd)
}

```

