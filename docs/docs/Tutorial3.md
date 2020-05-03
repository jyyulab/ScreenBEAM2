---
layout: default
title: Tutorial
nav_order: 4

---

# Tutorial-3 Differential representative analysis
{:.no_toc}
ScreenBEAM2 is a R based tool which consists of three major parts for processing steps: 1. mapping long read sequence to short read libraries, 2. Quality control, data cleanning and data preprocessing for mapped raw counts data; 3. Differential representative analysis on gene level or shRNA level.

**This is the third part of the whole tutorial, which is focused on differential representative analysis on gene and shRNA.**


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## Step3.Perform pairwise comparisons on Gene level
First, define your comparison pairs and make sure all group names are correct. Here as sample analysis, we used a 
```R
unique.group<-unique(m$samples$group)
compare.pairs<-combn(unique.group,2)
compare.pairs
```

Here we use function `ScreenBEAM.gene.level` in ScreenBEAM2 to perform gene level 
```
# Define your input file
input.file<-paste0(analysis.par$out.dir.output.DR,"YOUR_PROJECT_9mm_normalized.tsv")
de.gene.list<-list()

for(i in 1:dim(compare.pairs)[2]){
  
  case.groupname<-compare.pairs[1,i]
  control.groupname<-compare.pairs[2,i]
  
  case.sample.id<-which(m$samples$group==compare.pairs[1,i])
  control.sample.id<-which(m$samples$group==compare.pairs[2,i])
  
  case.postfix<-LETTERS[seq(from=1, to=length(case.sample.id))]
  control.postfix<-LETTERS[seq(from=1, to=length(control.sample.id))]
  
  control.samples<-paste(m$samples$group[control.sample.id], control.postfix, sep='_')
  case.samples<-paste(m$samples$group[case.sample.id], case.postfix, sep='_')
  compare.name<-paste0(case.groupname,".vs.",control.groupname)
  print(compare.name)
  
  de <- ScreenBEAM.gene.level(input.file, control.samples = control.samples, case.samples = case.samples, data.type = "NGS",control.groupname=control.groupname, case.groupname=case.groupname,gene.columnId=2, do.normalization=FALSE, filterLowCount=TRUE,filterBy='control',count.cutoff=4, rna.size = 6, sample.rna.time = 100, pooling = "partial", method='Bayesian')
  
  names(de)[1]<-"geneID"
  de.gene.list[[compare.name]]<-de
}


if(length(de.gene.list)>1){
  DR.GENE.DF<-de.gene.list[[1]]
  for(i in 2:length(de.gene.list)){
    DR.GENE.DF<-merge(DR.GENE.DF, de.gene.list[[i]], by="ID")
  }
}else{
  DR.GENE.DF<-de.gene.list[[1]]
}


DR.GENE.DF.sel<-dplyr::select(DR.GENE.DF,
                      geneID,
                      starts_with('log2FC.'),
                      starts_with('z.'),
                      starts_with('pval.'),
                      starts_with('FDR.'),
                      starts_with('B.'),
                      starts_with('B.sd'),
                      starts_with('n.sh_sgRNAs.passFilter.')
)

write.xlsx(DR.GENE.DF.sel, paste0(analysis.par$out.dir, lib.name,".DR.GENE.xlsx"))
DR.GENE.DF.sel$geneID<-as.character(DR.GENE.DF.sel$geneID)

# Take 1 comparison as example, to return the list of significant genes and draw a plot
sig_gene <- draw.volcanoPlot(dat=DR.GENE.DF.sel, label_col = "geneID", logFC_col = names(DR.GENE.DF.sel)[2], Pv_col = names(DR.GENE.DF.sel)[5], logFC_thre = 2, Pv_thre = 1e-4, main = names(de.rna.list), show_label = T, label_cex = 1, pdf_file =  paste0(analysis.par$out.dir.output.DR, names(de.gene.list),"_GENE.pdf"))
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
                           gene.columnId=2, do.log2=TRUE, do.normalization=FALSE, total=1e6, 
                           filterLowCount=TRUE,
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



