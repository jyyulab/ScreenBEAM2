library(MAGeCKFlute)
####################
source('src/functions.R')
demoInput_dir  <- '/home/dongxinran/project/screenBEAM2/dataset/2013NBT/'
demoOutput_dir <- '/home/dongxinran/project/screenBEAM2/dataset_result/D4/MAGeCKFlute/'
if (!file.exists(demoOutput_dir)) dir.create(demoOutput_dir, recursive = TRUE)
project_name <- 'mouseESC' # patient-derived GBM stem-like cells (GSCs) 
demoInput_fastq_dir <- demoInput_dir
library_csv <- sprintf('%s/library.csv',demoInput_dir)
Metadata_csv <- sprintf('%s/Metadata.csv',demoInput_dir)
output_prefix <- sprintf('%s/%s',demoOutput_dir,project_name)
org <- 'mmu';case_group <- 'Hi';control_group <- 'Lo'

########## A. Process CRISPR screen data step by step with MAGeCK
##### (i) Download and unzip the test data for both datasets
##### (ii) Generate a count table for Dataset 1
##### (iii) run the mageck count on Dataset
meta_data <- read.csv(Metadata_csv)
sample_label <- paste(meta_data$sampleLabel,collapse = ',')
fastq_file <- paste(meta_data$fastqFile,collapse = ' ')
fastq_file <- paste(sprintf('%s/%s',demoInput_fastq_dir,meta_data$fastqFile),collapse = ' ')
cmd <- sprintf('mageck count -l %s -n %s --sample-label %s --fastq %s',
				library_csv, output_prefix, sample_label, fastq_file)
print(cmd)
#system(cmd)

##### (iv) Batch effect removal
##### (v) Identify screen hits using MAGeCK RRA
count_file <- sprintf('%s.count.txt',output_prefix)
sample_label_case <- paste(meta_data$sampleLabel[which(meta_data$group == case_group)],collapse = ',')
sample_label_control <- paste(meta_data$sampleLabel[which(meta_data$group == control_group)],collapse = ',')
cmd <- sprintf('mageck test -k %s -t %s -c %s -n %s_rra --remove-zero both --remove-zero-threshold 0',
			   count_file, sample_label_case, sample_label_control,output_prefix)
print(cmd)
#system(cmd)

##### (vi/vii) Identify screen hits using MAGeCK MLE; Identify screen hits using MAGeCK MLE. If an experiment contains more than two conditions, for example, a three-condition design: day 0, drug treatment and DMSO treatment, we recommend using MAGeCK MLE or MAGeCK-NEST (Box 2) instead of MAGeCK RRA to perform the previous step

##### (viii/ix) Correct copy-number bias. MAGeCK RRA and MAGeCK MLE contain an optional method to correct copy-number biases in the calculated RRA scores and beta scores, respectively; if the CNV information is available for the cell line. 
# SYMBOL  HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
# A1BG    0.1105
# NAT2    0.0104

########## B. Process CRISPR screen data with MAGeCK-VISPR
###### A. Functional analysis for MAGeCK RRA results
gene_summary_file <- sprintf('%s_rra.gene_summary.txt',output_prefix)
FluteRRA(gene_summary = gene_summary_file,organism=org, outdir = demoOutput_dir, proj = project_name)

###### B. Functional analysis of MAGeCK MLE results
###### generate report
currentdir <- getwd()
setwd(demoOutput_dir)
Rspt <- list.files('.',pattern = '*R$')
for(eachRspt in Rspt) system(sprintf('Rscript %s',eachRspt))
setwd(currentdir)









