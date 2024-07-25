D1. https://github.com/sbodapati/CRISPR_Benchmarking_Algorithms
   2020. A benchmark of algorithms for the analysis of pooled CRISPR screens
   -- Simulation
     * Number of guides per gene
     * Number of control guides
     * Sequencing depth
     * Gene effect size
     * Guide binding efficiency
   -- Real dataset
     * TKO results

D2. https://github.com/DepMap-Analytics/HT29benchmark
R package for benchmarking genome-wide CRISPR-Cas9 knock-out viability screening pipelines making use of a reference dataset from six high-quality screens of the HT29 cell line
c903(rep 6) c904(rep 3) c905(rep 9) c906(rep 6) c907(rep 3) c908(rep 3)

D3. http://cistrome.org/MAGeCKFlute/demo.tar.gz (fastq)

D3.1 A CRISPR screen dataset generated from patient-derived GBM sGSCs (Gene Expression Omnibus (GEO) accession number: GSE70038)20 (Dataset 1). In this screen, cells were collected at two time points: day 0 (initial time point of the screen) and day 23 (after 23 d of culture). Replicate 1 and replicate 2 are biological replicates. This dataset is in FASTQ format. These data are used to demonstrate how to analyze screen data using MAGeCK RRA.

D3.2 The second dataset is a CRISPR screen in a melanoma cancer cell line, A375, treated with PLX7(Dataset 2). In this case, the cells were collected at two time points, 7 and 14 d after treatment, and were compared to a control (DMSO-treated) condition. The A375 dataset provides a raw readcount table for each sgRNA. These data contain three conditions (day 0, DMSO treatment and drug treatment) and are used to demonstrate how to analyze a screen with more than two conditions using MAGeCK MLE.

D3.3 The HCT116 dataset is a genome-wide CRISPR screen using HCT116 colorectal carcinoma cells8. The sgRNAs from several time points were collected and sequenced. This dataset was generated in different batches, and a read-count table is included in the demo data. We use this dataset to demonstrate how to perform batch effect removal.

D3.4 The HL60 dataset is a genome-wide CRISPR screen using the acute myelocytic leukemia HL60 cell line42 with copy-number variation (CNV) information. Cells were collected at two time points: HL60_initial (initial time point of the screen) and HL60_final (after 12 doubling times). We use this screen to demonstrate copy-number bias correction. (Genetic Screens in Human Cells Using the CRISPR-Cas9 System)
D3.5 The LNCap dataset (Supplementary Data 4) is a genome-wide CRISPR screen dataset from two cell lines, LNCap95 and LNCap abl43. Both contain three conditions (day 0, DMSO treatment and drug treatment) and include AAVS1-targeting sgRNAs as negative controls, so they are used to demonstrate the normalization with AAVS1-targeting sgRNAs.

D4. lingyun6sample
mouseESC
They designed 87,897 guide RNAs (gRNAs) targeting 19,150 mouse protein-coding genes and screening the resulting ESC mutant libraries for resistance to either Clostridium septicum alpha-toxin or 6-thioguanine, identified 27 known and 4 previously unknown genes implicated in these phenotypes. 

