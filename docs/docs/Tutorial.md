---
layout: default
title: Tutorial
nav_order: 3

---

# Tutorial
{:.no_toc}
ScreenBEAM2 is a R based tool which consists of three major parts for processing steps: 1. mapping long read sequence to short read libraries, 2. Quality control, data cleanning and data preprocessing for mapped raw counts data; 3. Differential representative analysis on gene level or shRNA level.


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## Step1.Mapping to library 
### Step1.1 Prepare raw fastq and library files
First of all, define your project name and create a project R object. This could be achieved by function `ScreenBEAM.dir.create`. 

**note: Project object should be library based! One library per run of ScreenBEAM2!**

```R
lib.name<-'[your_library_name_+_project_name]'
analysis.par<-ScreenBEAM.dir.create(project_main_dir = './', lib_name = lib.name, DATE = T)

analysis.par$par.path<-analysis.par$out.dir
analysis.par$par.name<-paste0("analysis_par_",lib.name, ".RData")

load(paste0(par.path, par.name))
```
Create a table with all raw fastq files that are mapped to this library
```R
f<-list.files(path="../201901/rawdata/Das_fastq/Lib_B/",patter=".fastq$",recursive=TRUE,full.names=TRUE)
f<-normalizePath(f)
write.table(f,file=file.path(par.path,"fastq_paths.txt"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

```

Next, soft link your fastq file to `analysis.par$out.dir.fastq` folder. This step is **ESSENTIAL**, please **DON'T** move files directly, in case of losing your data. You can go to your `analysis.par$out.dir.fastq` folder, then run command:

```Shell
ln -s filepaths
```
Put library csv file with file name "JJ_Das_LibB.csv" in the `analysis.par$out.dir.library` folder. The 1st column must be gRNA name and 2nd column must be sequence. Then go to `analysis.par$out.dir.library` folder and run following command to create fasta file for library.

```Shell
awk FNR-1 your_library.csv | awk -F "," '{print ">"$1"\n"$2}' > your_library.fa" 
```
Then you could move forward to mapping and collecting raw counts step, which could be executed on cluster.


### Step1.2 Mapping and collecting raw counts

first of all, please ensure `blat` was installed in your working environment, since we will be utilizing `blat` for mapping. (You can use blat version here on conda)

Next, we wil need to execute `ScreenBEAM.raw.count` function in ScreenBEAM2 R package. This function will help mapping your fastq files to library and collect raw counts by calling BLAT from R.

```R

analysis.par<-ScreenBEAM.raw.count(analysis.par)

# suggest to save analysis.par now, because a lot of data has been created.
save(analysis.par, file = paste0(par.path, par.name))

```
**NOTE: Some fastq files are ultra big, you may need to chop them into small fastq files in order to execute them successfully.**

### Step1.3 Mapping qualtiy control

After mapping is done, please proceed your analysis back in interactive environment (such as R studio), and run command as follows:
```R
ScreenBEAM.mapping.QC(analysis.par)
```
This will output a QC report for your library mapping, including mapping rate, mismatch rate, count boxplots, etc.



## Step2 Data annotation and data cleanning
### Step2.1 Create Metadata 
First of all, please have your Metadata information collected in an excel file, and saved it to `analysis.par$out.dir.metadata`. You could take following excel file as an example. Please note that, ID that could be matched to sampleID on fastq files is required,as well as proper group informations.

```R
meta <-read.xlsx(file.path(analysis.par$out.dir.metadata,"Metadata.xlsx"), sheet = 1)
dim(meta)
raw<-analysis.par$raw.count.table
s.cur<-data.frame(sampleID.full=names(raw))
#chop sample name into sectors
s.cur<-mutate(s.cur,
              sampleID=gsub('(.*)_(.*)_(.*)_(.*)_R1_001','\\2',sampleID.full),
              HartwellID=as.numeric(gsub('(.*)_(.*)_(.*)_(.*)_R1_001','\\1',sampleID.full)),
              HiseqSampleID=gsub('(.*)_(.*)_(.*)_(.*)_R1_001','\\3',sampleID.full),
              HiseqLaneID=gsub('(.*)_(.*)_(.*)_(.*)_R1_001','\\4',sampleID.full)
)#this could be adjusted accordingly

meta<-left_join(s.cur,meta,by='sampleID')
dim(meta)
```

Check if there are duplicates in your metadata -- this step is crucial because duplicates will affect following data manipulation.
```
meta<-mutate(meta,sampleName=paste(sampleID,HiseqLaneID,sep='_'),index=paste(HiseqSampleID,HiseqLaneID,sep='_'))
filter(meta,duplicated(sampleName))
filter(meta,duplicated(index))
dim(meta)
```

### Step2.2 Combine mapping info to your metadata
Read mapping info and combine them in meta data:
```R
mapping<-analysis.par$raw.summary
mapping<-mapping[,c('sample','total','n.matched')]
names(mapping)[1]<-'sampleID.full'
mapping<-mutate(mapping,mappingRate=n.matched/total,coverage.seq=round(n.matched/nrow(analysis.par$raw.count.table)))
head(mapping) # calculate mapping statistics
meta<-left_join(meta,mapping,by="sampleID.full")
unlink(file.path(analysis.par$out.dir.output.mapping,'mapping.summary.xlsx'))
write.xlsx(meta,file=file.path(analysis.par$out.dir.output.mapping,'mapping.summary.xlsx'),sheet='mapping')
```

Add n.mismatch filter and Save eset raw data 
```R
count.dist<-as.data.frame(analysis.par$raw.count.dist)
count.dist$sample<-meta$sampleName[match(count.dist$sample,meta$sampleID.full)]
names(count.dist)[2]<-'sampleName'
```

### Step2.3 Data normalization and create Expressionset
For this step you will need to perform data normalization and store different normalization method into eset.

```R
normalize.total<-1e6

m<-list(samples=meta,features=analysis.par$lib, count.dist=count.dist)

save(m, file = paste0(par.path, "raw.list.Rdata"))
saveCountEset(m,
	save.path=analysis.par$out.dir.output.mapping, 
	n.mismatch = 9) # n.mismatch should be adjusted according to mapping qc

```

### Step2.4 Qualtiy control for wrapped expressionsets
Here we have a build-in QC function `draw.eset.QC` to facilitate quality control purpose. Each `draw.eset.QC` function will output a quality control report for following expression set.

```R
load(paste0(analysis.par$out.dir.output.mapping,"count.9mm.raw.eset"));
raw.9mm<-count.Nmm.raw.eset # assign your eset with a more specific name
draw.eset.QC(raw.9mm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'raw.9mm_',
             choose_plot = c("heatmap", "pca","density","correlation","meansd"))

load(paste0(analysis.par$out.dir.output.mapping,"count.9mm.normalized.eset.1M.eset"));
norm.9mm<-count.Nmm.normalized.eset
draw.eset.QC(norm.9mm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'norm.9mm_',
             choose_plot = c("heatmap", "pca","density","correlation","meansd"))

load(paste0(analysis.par$out.dir.output.mapping,"count.maxmm.raw.eset"));
raw.maxmm<-count.maxmm.raw.eset
draw.eset.QC(raw.maxmm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'raw.maxmm_',
             choose_plot = c("heatmap", "pca","density","correlation","meansd"))

load(paste0(analysis.par$out.dir.output.mapping,"count.maxmm.normalized.1M.eset"));
norm.maxmm<-count.maxmm.normalized.eset
draw.eset.QC(norm.maxmm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'norm.maxmm_',
             choose_plot = c("heatmap", "pca","density","correlation","meansd"))

```
### Step2.5 Prepare tsv file for differential representation analysis
Save normalized data and ids (gene/shRNA) in a tsv file for downstream analysis.

```R
first.2.column<-data.frame(RNAid=analysis.par$lib$id, geneid=analysis.par$lib$gene)
count.column<-as.data.frame(exprs(norm.9mm))
count.column$RNAid<-rownames(count.column)

final.table<-merge(first.2.column, count.column, by="RNAid")

colnames(final.table)<-c("rnaID","geneID", paste(meta$group,meta$replicate, sep = "_"))
write.table(final.table, file = paste0(analysis.par$out.dir.output.DR,"YOUR_PROJECT_9mm_normalized.tsv"), quote = F, row.names = F, sep='\t')

```

## Step3.Perform pairwise comparisons on Gene level
For whole genome libraries, it may cost sometime to run these pairwise comparisons sequentially. So we recommend to generate one R script for each comparisons and running them in parallel.











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



