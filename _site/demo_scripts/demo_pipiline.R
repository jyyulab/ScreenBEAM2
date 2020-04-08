### Author: Jiyang Yu, Xinge Wang
### 2019/08/02

### This demo only use 6 samples from 201805 folder of the lingyun data
### Sample ID: F171, F173, F175 (3 replicates in 0.25.pS6.Hi group); F172, F174, F176 (3 replicates in 0.25.pS6.Lo group)
### For internal use only now

library(readxl)
library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(data.table)
library(NetBID2)
library(RColorBrewer)
library(NetBID2)

# If you don't install ScreenBEAM2 package, please source the following code
source("/rgs01/project_space/yu3grp/FuncGenomics/yu3grp/training/ScreenBEAM2_Rpackage/demo_data/ScreenBEAM2_functions.R")
source("/rgs01/project_space/yu3grp/FuncGenomics/yu3grp/training/ScreenBEAM2_Rpackage/demo_data/ScreenBEAM2_utilities.R")


################ Make directory for library#############
lib.name<-'lingyun6sample'
analysis.par<-ScreenBEAM.dir.create(project_main_dir = './', lib_name = lib.name, DATE = T)
par.path<-analysis.par$out.dir
par.name<-paste0("analysis_par_",lib.name, ".RData")
save(analysis.par, file = paste0(par.path, par.name))
load(paste0(par.path, par.name))

################### Prepare fastq and library files #####################
### This step is ESSENTIAL, we suggest to COPY or SOFT LINK the original data files. DON'T move files directly, in case of lose your data.
# 1. put unzipped fastq files in analysis.par$out.dir.fastq folder
# (First, go to analysis.par$out.dir.fastq folder, then run the following commands in linux:
#ln -s /rgs01/project_space/yu3grp/FuncGenomics/yu3grp/ChiLab/CRISPRscreens/201805/rawdata/chigrp_131847_crispr_tagged_amplicon/1382255/1382255_F171_S3_L003_R1_001.fastq
#ln -s /rgs01/project_space/yu3grp/FuncGenomics/yu3grp/ChiLab/CRISPRscreens/201805/rawdata/chigrp_131847_crispr_tagged_amplicon/1382256/1382256_F172_S4_L003_R1_001.fastq
#ln -s /rgs01/project_space/yu3grp/FuncGenomics/yu3grp/ChiLab/CRISPRscreens/201805/rawdata/chigrp_131847_crispr_tagged_amplicon/1382257/1382257_F173_S5_L003_R1_001.fastq
#ln -s /rgs01/project_space/yu3grp/FuncGenomics/yu3grp/ChiLab/CRISPRscreens/201805/rawdata/chigrp_131847_crispr_tagged_amplicon/1382258/1382258_F174_S6_L003_R1_001.fastq
#ln -s /rgs01/project_space/yu3grp/FuncGenomics/yu3grp/ChiLab/CRISPRscreens/201805/rawdata/chigrp_131847_crispr_tagged_amplicon/1382259/1382259_F175_S7_L003_R1_001.fastq
#ln -s /rgs01/project_space/yu3grp/FuncGenomics/yu3grp/ChiLab/CRISPRscreens/201805/rawdata/chigrp_131847_crispr_tagged_amplicon/1382260/1382260_F176_S8_L003_R1_001.fastq
# 2. put library csv file with file name "lingyun6sample" in the analysis.par$out.dir.library folder. The 1st column must be gRNA name and 2nd column must be sequence. Suggested column ids are = c("id", "seqRNA", "Gene"). In Yulab cluster, go to analysis.par$out.dir.library folder first, then run "cp /rgs01/project_space/yu3grp/FuncGenomics/yu3grp/training/ScreenBEAM2_Rpackage/demo_data/lingyun.csv ./"
# 3. Convert the csv library file to .fa file for blat. cd to analysis.par$out.dir.library folder, and run command "awk FNR-1 lingyun.csv | awk -F "," '{print ">"$1"\n"$2}' > lingyun.fa' to create fasta file for library. And check if the lingyun6sample.fa file has shRNA as fasta read name, and shRNA sequence as sequence.
# 4. meta data, please put meta data excel file in the analysis.par$out.dir.metadata folder. The meta data will be used after the raw count table is prepared. In Yulab cluster,
#go to run: "cp /rgs01/project_space/yu3grp/FuncGenomics/yu3grp/training/ScreenBEAM2_Rpackage/demo_data/Metadata.xlsx ./"

################### Mapping and collecting raw counts ####################
load(paste0(par.path, par.name))
# Can run this step on cluster, and save it as analysis.par object, and run the rest analysis on local Rstudio
analysis.par<-ScreenBEAM.raw.count(analysis.par)
# suggest to save analysis.par now, because a lot of data has been created.
save(analysis.par, file = paste0(par.path, par.name))
load(paste0(par.path, par.name))

################### QC report for n.mismatch threshold ####################
## if you don't have ScreenBEAM2 package installed, use the following code
#output_rmd_file <- sprintf("%s/mapping_%s_QC.Rmd", analysis.par$out.dir.output.QC, analysis.par$lib.name)
#file.copy(from = "path/to/mapping_QC_report.Rmd", to = output_rmd_file, overwrite = T)
#rmarkdown::render(output_rmd_file, html_document(toc = TRUE))

## If you have ScreenBEAM2 package stored already
ScreenBEAM.mapping.QC(analysis.par, QC.Rmd.path = "/rgs01/project_space/yu3grp/FuncGenomics/yu3grp/training/ScreenBEAM2_Rpackage/demo_run/mapping_QC_report.Rmd") # Will create QC report of mapping in analysis.dir$out.dir.QC folder

################### Add meta data, this pipeline need to be modified case to case, because sample names are different ####################
# pheno type data
meta<-read.xlsx(file.path(analysis.par$out.dir.metadata,"Metadata.xlsx"), sheet = 1)
dim(meta)
raw<-analysis.par$raw.count.table
s.cur<-data.frame(sampleID.full=names(raw))
#chop sample name into sectors
s.cur<-mutate(s.cur,
              sampleID=gsub('(.*)_(.*)_(.*)_(.*)_R1_001','\\2',sampleID.full),
              HartwellID=as.numeric(gsub('(.*)_(.*)_(.*)_(.*)_R1_001','\\1',sampleID.full)),
              HiseqSampleID=gsub('(.*)_(.*)_(.*)_(.*)_R1_001','\\3',sampleID.full),
              HiseqLaneID=gsub('(.*)_(.*)_(.*)_(.*)_R1_001','\\4',sampleID.full)
)
meta<-left_join(s.cur,meta,by='sampleID')
dim(meta)
#check if there are duplicates
meta<-mutate(meta,sampleName=paste(sampleID,HiseqLaneID,sep='_'),index=paste(HiseqSampleID,HiseqLaneID,sep='_'))
filter(meta,duplicated(sampleName))
filter(meta,duplicated(index))
dim(meta)
###read mapping info
mapping<-analysis.par$raw.summary
mapping<-mapping[,c('sample','total','n.matched')]
names(mapping)[1]<-'sampleID.full'
mapping<-mutate(mapping,mappingRate=n.matched/total,coverage.seq=round(n.matched/nrow(analysis.par$raw.count.table)))
head(mapping)
meta<-left_join(meta,mapping,by="sampleID.full")
unlink(file.path(analysis.par$out.dir.output.mapping,'mapping.summary.xlsx'))
write.xlsx(meta,file=file.path(analysis.par$out.dir.output.mapping,'mapping.summary.xlsx'),sheet='mapping')

############# Add n.mismatch filter and Save eset raw data and normalized data #################
count.dist<-as.data.frame(analysis.par$raw.count.dist)
count.dist$sample<-meta$sampleName[match(count.dist$sample,meta$sampleID.full)]
names(count.dist)[2]<-'sampleName'
normalize.total<-1e6
m<-list(samples=meta,features=analysis.par$lib, count.dist=count.dist)
save(m, file = paste0(analysis.par$out.dir, "raw.list.Rdata"))
saveCountEset(m,save.path=analysis.par$out.dir.output.mapping, n.mismatch = 9)

################### QC report for ESET ####################
load(paste0(analysis.par$out.dir.output.mapping,"count.9mm.raw.eset"));
raw.9mm<-count.Nmm.raw.eset
draw.eset.QC(raw.9mm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'raw.9mm_',
             choose_plot = c("heatmap", "pca","density","correlation"))

load(paste0(analysis.par$out.dir.output.mapping,"count.9mm.normalized.eset.1M.eset"));
norm.9mm<-count.Nmm.normalized.eset
draw.eset.QC(norm.9mm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'norm.9mm_',
             choose_plot = c("heatmap", "pca","density","correlation"))

load(paste0(analysis.par$out.dir.output.mapping,"count.maxmm.raw.eset"));
raw.maxmm<-count.maxmm.raw.eset
draw.eset.QC(raw.maxmm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'raw.maxmm_',
             choose_plot = c("heatmap", "pca","density","correlation"))

load(paste0(analysis.par$out.dir.output.mapping,"count.maxmm.normalized.1M.eset"));
norm.maxmm<-count.maxmm.normalized.eset
draw.eset.QC(norm.maxmm, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', do.logtransform = T, prefix = 'norm.maxmm_',
             choose_plot = c("heatmap", "pca","density","correlation"))

################### Prepare tsv file for runnint ScreenBEAM gene/rna level data ####################
#### This step match how the ScreenBEAM1 works
first.2.column<-data.frame(RNAid=analysis.par$lib[,1], geneid=analysis.par$lib[,3])
count.column<-as.data.frame(exprs(norm.9mm))
count.column$RNAid<-rownames(count.column)
final.table<-merge(first.2.column, count.column, by="RNAid")
colnames(final.table)<-c("rnaID","geneID", paste(meta$group,meta$replicate, sep = "_"))
write.table(final.table, file = paste0(analysis.par$out.dir.output.DR,"lingyun_9mm_normalized.tsv"), quote = F, row.names = F, sep='\t')

################## Pairwise comparison, RNA level ############################3

unique.group<-unique(m$samples$group)
compare.pairs<-combn(unique.group,2)
compare.pairs
input.file<-paste0(analysis.par$out.dir.output.DR,"lingyun_9mm_normalized.tsv")

de.rna.list<-list()
for(i in 1:dim(compare.pairs)[2]){
  case.groupname<-compare.pairs[1,i]
  control.groupname<-compare.pairs[2,i]
  case.sample.id<-which(m$samples$group==compare.pairs[1,i])
  control.sample.id<-which(m$samples$group==compare.pairs[2,i])
  case.postfix<-LETTERS[seq(from=1, to=length(case.sample.id))]
  control.postfix<-LETTERS[seq(from=1, to=length(control.sample.id))]
  control.samples<-paste(m$samples$group[control.sample.id], control.postfix, sep='_')
  case.samples<-paste(m$samples$group[case.sample.id], case.postfix, sep='_')
  de<-ScreenBEAM.rna.level(input.file, control.samples = control.samples, case.samples = case.samples,
                           control.groupname=control.groupname, case.groupname=case.groupname,
                           gene.columnId=2, logTransformed=TRUE, do.log2=TRUE, do.normalization=TRUE, total=1e6, filterLowCount=TRUE,
                           filterBy='control',count.cutoff=32,family=gaussian,estimation.method='Bayesian')
  compare.name<-paste0(case.groupname,".vs.",control.groupname)
  names(de)<-paste(names(de), compare.name, sep = ".")
  names(de)[1]<-"rnaID"
  de.rna.list[[compare.name]]<-de
}

if(length(de.rna.list)>1){
  DR.RNA.DF<-de.rna.list[[1]]
  for(i in 2:length(de.rna.list)){
    DR.RNA.DF<-merge(DR.RNA.DF, de.rna.list[[i]], by="ID")
  }
}else{
  DR.RNA.DF<-de.rna.list[[1]]
}

DR.RNA.DF.sel<-dplyr::select(DR.RNA.DF,
                      rnaID,
                      starts_with('logFC'),
                      starts_with('Z-statistics'),
                      starts_with('AveExpr'),
                      starts_with('P.Value.'),
                      starts_with('adj.P.Val'),
                      starts_with('Ave.')
)
write.xlsx(DR.RNA.DF.sel, paste0(analysis.par$out.dir, lib.name,".DR.sgRNA.xlsx"))
DR.RNA.DF.sel$rnaID<-as.character(DR.RNA.DF.sel$rnaID)
# Take 1 comparison as example, to return the list of significant RNAs and draw a plot
sig_rna<-draw.volcanoPlot(dat=DR.RNA.DF.sel, label_col = "rnaID", logFC_col = paste0("logFC.",names(de.rna.list)),
                          Pv_col = paste0("adj.P.Val.",names(de.rna.list)), logFC_thre = 3, Pv_thre = 1e-2,
                          main = names(de.rna.list), show_label = T, label_cex = 1, pdf_file = paste0(analysis.par$out.dir.output.DR, names(de.rna.list),"_RNA.pdf"))

################## Pairwise comparison, GENE level ############################
unique.group<-unique(m$samples$group)
compare.pairs<-combn(unique.group,2)
compare.pairs
input.file<-paste0(analysis.par$out.dir.output.DR,"lingyun_9mm_normalized.tsv")
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
  de<-ScreenBEAM.gene.level(input.file, control.samples = control.samples, case.samples = case.samples, data.type = "NGS",
                            control.groupname=control.groupname, case.groupname=case.groupname,
                            gene.columnId=2, do.normalization=TRUE, filterLowCount=TRUE, filterBy='control',count.cutoff=4, rna.size = 6, sample.rna.time = 100, pooling = "partial",
                            method='Bayesian')
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
sig_gene<-draw.volcanoPlot(dat=DR.GENE.DF.sel, label_col = "geneID", logFC_col = names(DR.GENE.DF.sel)[2],
                          Pv_col = names(DR.GENE.DF.sel)[5], logFC_thre = 2, Pv_thre = 1e-4,
                          main = names(de.rna.list), show_label = T, label_cex = 1, pdf_file = paste0(analysis.par$out.dir.output.DR, names(de.gene.list),"_GENE.pdf"))


########## Trim helper plot function test ############
ScreenBEAM.trim.helper(fastq.path = paste0(analysis.par$out.dir.fastq, "1382255_F171_S3_L003_R1_001.fastq"), sample_num = 1000, pdf.file = paste0(analysis.par$out.dir,"trim_plot.pdf"))

