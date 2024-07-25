library(ScreenBEAM2)
####################
source('src/functions.R')
demoInput_dir  <- '/home/dongxinran/project/screenBEAM2/dataset/demo/mageck_cnv/'
demoOutput_dir <- '/home/dongxinran/project/screenBEAM2/dataset_result/D3.4/screenBEAM2/'
if (!file.exists(demoOutput_dir)) dir.create(demoOutput_dir, recursive = TRUE)
project_name <- 'HL60' #  The HL60 dataset is a genome-wide CRISPR screen using the acute myelocytic leukemia HL60 cell line42 with copy-number variation (CNV) information. Cells were collected at two time points: HL60_initial (initial time point of the screen) and HL60_final (after 12 doubling times). We use this screen to demonstrate copy-number bias correction
count_file <- sprintf('%s/rawcount.txt',demoInput_dir) # sgRNA   Gene    HL60.initial    KBM7.initial    HL60.final  KBM7.final
cnv_file <- sprintf('%s/cnv_data.txt',demoInput_dir) # SYMBOL  HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
Metadata_csv <- sprintf('%s/Metadata.csv',demoInput_dir)
output_prefix <- sprintf('%s/%s',demoOutput_dir,project_name)
org <- 'hsa';case_group <- 'HL60.final';control_group <- 'HL60.initial'
cellline_name <- 'HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'

####################
#### Step1.Mapping to library
### Step1.1 Prepare raw fastq and library files
# rm(analysis.par) # if re-run, need to remove analysis.par
analysis.par <- ScreenBEAM.dir.create(project_main_dir = demoOutput_dir, lib_name = project_name, DATE = T,
				library_file = NULL, metadata_file = Metadata_csv, count_file = count_file,
				startFromFastq = FALSE)

### Step1.2 Mapping and collecting raw counts
### Step1.3 Mapping qualtiy control
#### Step2 Data annotation and data cleanning
### Step2.1 Data normalization, create Expressionset and Prepare tsv file for differential representation analysis
### Step2.2 Qualtiy control for wrapped expressionsets

#### Step3.Perform pairwise comparisons + CNV information correction 
use_index <- 1 # check analysis.par$norm.path to select
## Gene level
analysis.par <- ScreenBEAM.Pairwise(analysis.par,choose_level = 'gene', 
				use_index = use_index,
				gene.columnId = 2, data.type = 'NGS',
				case_group = case_group, control_group = control_group,
				do.normalization = TRUE, total = 1e6, filterLowCount = TRUE,
				filterBy = "control", count.cutoff = 4, nitt = 15000, burnin = 5000, thin=10,
				rna.size=6, sample.rna.time=100, method = "Bayesian", pooling="partial")
# get sig gene
DR.GENE.DF.sel <- read.xlsx(analysis.par$norm.path$DR_gene_filepath[use_index])
compare.name <- analysis.par$norm.path$DR_compare[use_index]
sig_gene <- draw.volcanoPlot(dat=DR.GENE.DF.sel,
			label_col = "geneID", 
			logFC_col = names(DR.GENE.DF.sel)[2], 
			Pv_col = names(DR.GENE.DF.sel)[4], 
			logFC_thre = 0.25, Pv_thre = 0.05, 
			main = compare.name, show_label = F, 
			label_cex = 1, 
			pdf_file =  paste0(analysis.par$out.dir.output.DR, compare.name,"_GENE.pdf"))


## RNA level
analysis.par <- ScreenBEAM.Pairwise(analysis.par,choose_level = 'RNA', use_index = use_index, 
				case_group = case_group, control_group = control_group,
				do.normalization = FALSE, count.cutoff = 16, pooling = 'full', method = 'MLE')
# get sig RNA
DR.RNA.DF.sel <- read.xlsx(analysis.par$norm.path$DR_rna_filepath[use_index])
compare.name <- analysis.par$norm.path$DR_compare[use_index]
sig_rna <- draw.volcanoPlot(dat=DR.RNA.DF.sel,
            label_col = "rnaID",
            logFC_col = names(DR.RNA.DF.sel)[2],
            Pv_col = names(DR.RNA.DF.sel)[4],
            logFC_thre = 0.25, Pv_thre = 0.05,
            main = compare.name, show_label = F,
            label_cex = 1,
            pdf_file =  paste0(analysis.par$out.dir.output.DR, compare.name,"_RNA.pdf"))






