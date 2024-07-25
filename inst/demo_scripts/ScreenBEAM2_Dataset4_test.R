library(ScreenBEAM2)
####################
source('src/functions.R')
demoInput_dir  <- '/home/dongxinran/project/screenBEAM2/dataset/2013NBT/'
demoOutput_dir <- '/home/dongxinran/project/screenBEAM2/dataset_result/D4/screenBEAM2/'
if (!file.exists(demoOutput_dir)) dir.create(demoOutput_dir, recursive = TRUE)
project_name <- 'mouseESC' # patient-derived GBM stem-like cells (GSCs) 
library_csv <- sprintf('%s/library.csv',demoInput_dir)
Metadata_csv <- sprintf('%s/Metadata.csv',demoInput_dir)
case_group <- 'Hi';control_group <- 'Lo' # group column

####################
#### Step1.Mapping to library
### Step1.1 Prepare raw fastq and library files
# rm(analysis.par) # if re-run, need to remove analysis.par
analysis.par <- ScreenBEAM.dir.create(project_main_dir = demoOutput_dir, lib_name = project_name, DATE = T,
				library_file = library_csv, metadata_file = Metadata_csv)

### Step1.2 Mapping and collecting raw counts
analysis.par <- ScreenBEAM.raw.count(analysis.par)

### Step1.3 Mapping qualtiy control
ScreenBEAM.mapping.QC(analysis.par)

#### Step2 Data annotation and data cleanning
### Step2.1 Data normalization, create Expressionset and Prepare tsv file for differential representation analysis
analysis.par <- ScreenBEAM.createEset(analysis.par, normalize.total = 1e6, n.mismatch = 9)
analysis.par <- ScreenBEAM.createEset(analysis.par, normalize.total = 1e6, n.mismatch = 12)
analysis.par <- ScreenBEAM.createEset(analysis.par, normalize.total = 1e6, n.mismatch = 0)

### Step2.2 Qualtiy control for wrapped expressionsets
count.Nmm.raw.eset <- readRDS(analysis.par$norm.path$rawRData_filepath[1])
draw.eset.QC(count.Nmm.raw.eset, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', 
			do.logtransform = T, prefix = 'raw.9mm_', 
			choose_plot = c("heatmap", "pca","density","correlation","meansd"))
count.Nmm.normalized.eset <- readRDS(analysis.par$norm.path$RData_filepath[1])
draw.eset.QC(count.Nmm.normalized.eset, outdir = analysis.par$out.dir.output.QC, intgroup = 'group', 
			do.logtransform = T, prefix = 'norm.9mm_', 
			choose_plot = c("heatmap", "pca","density","correlation","meansd"))

#### Step3.Perform pairwise comparisons 
use_index = 1 # check analysis.par$norm.path to select
use_index <- 4
## Gene level
analysis.par <- ScreenBEAM.Pairwise(analysis.par,choose_level = 'gene', 
				use_index = use_index,
				gene.columnId = 2, data.type = 'NGS',
				case_group = case_group, control_group = control_group,
				do.normalization = FALSE, total = 1e6, filterLowCount = TRUE,
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






