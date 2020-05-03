---
layout: default
title: Tutorial - Cleanning
nav_order: 4

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
### Step2.1 Create Metadata 
First of all, please have your Metadata information collected in an excel file, and saved it to `analysis.par$out.dir.metadata`. Please note that, ID that could be matched to sampleID on fastq files is required, as well as proper group informations. A sample metadata file is shown as follows:

**sampleID**|**group**|description|replicate|expName|
| ------------- |:-------------:| -----:|-----:|-----:|
|F171| 0.25Hi| 0.25ug nacl-high | A | 2nd round screening in cells using library A
|F172| 0.25Lo| 0.25ug nacl-low  | A | 2nd round screening in cells using library A
|F173| 0.25Hi| 0.25ug nacl-high | B | 2nd round screening in cells using library A


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
)#this could be adjusted according to names of your fastq files

meta<-left_join(s.cur,meta,by='sampleID') # merge sample info from fastq with your meta data
dim(meta)
```

Create unique sample name and index for each sample by pasting sampleID and lane ID together. This could be adjusted accordingly.
```R
meta<-mutate(meta,sampleName=paste(sampleID,HiseqLaneID,sep='_'),  
	index=paste(HiseqSampleID,HiseqLaneID,sep='_')) 
```

Check if there are duplicates in your metadata -- this step is crucial because duplicates will affect following data manipulation. 

```R	
filter(meta,duplicated(sampleName))
filter(meta,duplicated(index))
dim(meta)
```

### Step2.2 Combine mapping info to your metadata
After all mapping jobs are completed and mapping statistics are generated and stored in designated location `analysis.par$raw.summary`. We could calculate mapping rate, and sequencing coverages from mapping statistics and store them in meta data.

```R
mapping<-analysis.par$raw.summary
mapping<-mapping[,c('sample','total','n.matched')]
names(mapping)[1]<-'sampleID.full'

#calculate mapping rate and sequencing coverages
mapping<-mutate(mapping,mappingRate=n.matched/total,coverage.seq=round(n.matched/nrow(analysis.par$raw.count.table)))
head(mapping) 

#combine mapping statistics with meta data, save files
meta<-left_join(meta,mapping,by="sampleID.full")
unlink(file.path(analysis.par$out.dir.output.mapping,'mapping.summary.xlsx'))
write.xlsx(meta,file=file.path(analysis.par$out.dir.output.mapping,'mapping.summary.xlsx'),sheet='mapping')
```

Add sample information to count distribution matrix.
```R
count.dist<-as.data.frame(analysis.par$raw.count.dist)
count.dist$sample<-meta$sampleName[match(count.dist$sample,meta$sampleID.full)]
names(count.dist)[2]<-'sampleName'
```
Now we are ready to deal with the count data and create expression sets.


### Step2.3 Data normalization and create Expressionset
Then, we could perform data normalization and store different normalization method into eset. Total number of normalization is data dependent, usually we use 1e6.

```R
normalize.total<-1e6
m<-list(samples=meta,features=analysis.par$lib, count.dist=count.dist)
save(m, file = paste0(par.path, "raw.list.Rdata"))
```

Save your expression matrix based on number of mismatch, this optimal mismatch number could be found in previous mapping QC report, as an example here we use 9. This function will also save an expression set with maximum mismatch and 0 mismatch, for your reference.

```R
saveCountEset(m,
	save.path=analysis.par$out.dir.output.mapping, 
	n.mismatch = 9) # n.mismatch should be adjusted according to mapping qc
```

### Step2.4 Qualtiy control for wrapped expressionsets
Here we have a build-in QC function imported from NetBID2 -- `draw.eset.QC` to facilitate quality control purpose. Each `draw.eset.QC` function will output a quality control report for following expression set.

```R
load(paste0(analysis.par$out.dir.output.mapping,"count.9mm.raw.eset"));
raw.9mm<-count.Nmm.raw.eset # assign your eset with a more specific name
draw.eset.QC(raw.9mm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'raw.9mm_', choose_plot = c("heatmap", "pca","density","correlation","meansd"))
```

Here we could also choose to perform more quality control for different numbers of mismatch, as a reference.

```R
load(paste0(analysis.par$out.dir.output.mapping,"count.9mm.normalized.eset.1M.eset"));
norm.9mm<-count.Nmm.normalized.eset
draw.eset.QC(norm.9mm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'norm.9mm_', choose_plot = c("heatmap", "pca","density","correlation","meansd"))

load(paste0(analysis.par$out.dir.output.mapping,"count.maxmm.raw.eset"));
raw.maxmm<-count.maxmm.raw.eset
draw.eset.QC(raw.maxmm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'raw.maxmm_', choose_plot = c("heatmap", "pca","density","correlation","meansd"))

load(paste0(analysis.par$out.dir.output.mapping,"count.maxmm.normalized.1M.eset"));
norm.maxmm<-count.maxmm.normalized.eset
draw.eset.QC(norm.maxmm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'norm.maxmm_', choose_plot = c("heatmap", "pca","density","correlation","meansd"))
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
After writing table, you will be able to quit current working environment and run following analysis on cluster. Or if your data size is moderate, you could stay in your r environment and move forward.


