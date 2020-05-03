---
layout: default
title: Tutorial
nav_order: 4

---

# Tutorial-3 Differential representative analysis
{:.no_toc}
ScreenBEAM2 is a R based tool which consists of three major parts for processing steps: 1. mapping long read sequence to short read libraries, 2. Quality control, data cleanning and data preprocessing for mapped raw counts data; 3. Differential representative analysis on gene level or shRNA level.




## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}


## Step3.Perform pairwise comparisons on Gene level
For whole genome libraries, it may cost sometime to run these pairwise comparisons sequentially. So we recommend to generate one R script for each comparisons and running them in parallel.

Define dependencies, function path (if package is not installed) , and input files.
```R
#source("dependencies.R")
func.path<-[your_path_to_functions] #if package is not installed
func.util.path<-"/research/projects/yu3grp/FuncGenomics/yu3grp/training/ScreenBEAM2_Rpackage/ScreenBEAM2/R/ScreenBEAM2_functions.R"
input.file<-"4cl_combined_withlib_normalized_50M.tsv"
input.file<-normalizePath(input.file)
dir.tmp<-"DR_hpc/"
```


```R
load([your_stored_eset])
d<-openxlsx::read.xlsx("20191020_CRISPR_comparisons.xlsx")
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



## Step4.Perform pairwise comparisons on shRNA level
shRNA level comparison also give another layer of information when gene level comparisons are not that informative. Statistical significance could only be computed when there are 2 or more replicates in experinment design, otherwise only folder changes (FC) will be computed successfully.

First of all, double check the group label in your comparisons are in `eset$group` info.
```R
unique.group<-unique(eset$group) # group label in your eset 
d<-openxlsx::read.xlsx("20191020_CRISPR_comparisons.xlsx") # comparisons that you want to perform
all(unlist(d)%in%unique.group) # make sure all group label exists in unique.group

compare.pairs<-t(d) # your comparisons are required to be as follows
```

Specify your input file:
```R
input.file<-"4cl_combined_withlib_normalized_50M.tsv"
```

Do comparisons by function `ScreenBEAM.rna.level`.

```R
de.rna.list<-list() # create an empty list to store following analysis results
```

```R
for(i in 1:dim(compare.pairs)[2]){
  case.groupname<-compare.pairs[1,i]
  control.groupname<-compare.pairs[2,i]
  
  case.sample.id<-which(eset$group==compare.pairs[1,i])
  control.sample.id<-which(eset$group==compare.pairs[2,i])
  
  case.postfix<-LETTERS[seq(from=1, to=length(case.sample.id))]
  control.postfix<-LETTERS[seq(from=1, to=length(control.sample.id))]
  
  control.samples<-eset$group[control.sample.id]
  case.samples<-eset$group[case.sample.id]
  
  de<-ScreenBEAM.rna.level(input.file, control.samples = control.samples, case.samples = case.samples,
                           control.groupname=control.groupname, case.groupname=case.groupname,
                           gene.columnId=2, do.log2=TRUE, do.normalization=FALSE, filterLowCount=TRUE,
                           filterBy='control',count.cutoff=8,family=gaussian,estimation.method='Bayesian')
  
  compare.name<-paste0(case.groupname,".vs.",control.groupname)
  names(de)<-paste(names(de), compare.name, sep = ".")
  names(de)[1]<-"rnaID"
  de.rna.list[[compare.name]]<-de
}

save.image("DE.sgRNA.RData") #save your data
```

Merge your results by "rnaID" and combine feature information from fData(eset).
```R
if(length(de.rna.list)>1){
  DR.RNA.DF<-de.rna.list[[1]]
  for(i in 2:length(de.rna.list)){
    DR.RNA.DF<-merge(DR.RNA.DF, de.rna.list[[i]], by="rnaID",all.x = TRUE)
  }
}else{
  DR.RNA.DF<-de.rna.list[[1]]
}

DR.RNA.DF<-merge(DR.RNA.DF,fData(eset),by.x="rnaID",by.y="id")
```


Clean up your data frame using function `select`. Save your analysis results.
```R
DR.RNA.DF.sel<-dplyr::select(DR.RNA.DF,
                      rnaID,
                      gRNA.sequence:nchar,
                      dplyr::starts_with('log2FC'),
                      dplyr::starts_with('Ave.')
)
save.image("DE.sgRNA.RData")
```


Output your result as an excel file.
```R
df<-exprs(eset)
DR.RNA.DF.sel<-merge(DR.RNA.DF.sel,df,by.x="rnaID",by.y="row.names")

openxlsx::write.xlsx(DR.RNA.DF.sel, "4CL.2Lib.DR.sgRNA.xlsx")

```



