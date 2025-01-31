# load project-specific parameters
config <- yaml::read_yaml(file = 'config.yaml')

# set knitr options
knitr::opts_chunk$set(comment = NA, fig.width = 7, fig.height = 4, out.width = '70%',
                      warning = TRUE, error = TRUE, echo = TRUE, message = TRUE,
                      dpi = 100)

library(data.table)
library(ggplot2)
library(simpleCache)
library(cowplot)
library(patchwork)

# set some other package-specific options
options(ggrepel.max.overlaps = Inf)

options(future.plan = 'sequential', 
        future.globals.maxSize = 8 * 1024 ^ 3)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 8, font_family = "Helvetica", rel_small = 1, rel_tiny = 6/8, rel_large = 1))

# set a random seed
set.seed(993751)

# store the current CPU and real time
SETUP_TIME <- proc.time()

# source any scripts with commonly used functions
source('R/utilities.R')

refGenome <- config$implied_attributes$organism$mouse$genome
refGenomeChrom <- config$implied_attributes$organism$mouse$genome_short
refGtf <- file.path(config$resources_root, "genomes", refGenome, paste0(refGenome, ".gtf")) 

dA <- loadAnnot()[protocol=="ATAC-seq", ] 
dA[, reads_file:=resultsDir("pipeline/", sample_name, "/aligned_", refGenomeChrom, "/", sample_name, "_sort_dedup.bam")] 
dA[, peaks_file:=resultsDir("pipeline/", sample_name, "/peak_calling_", refGenomeChrom, "/", sample_name, "_peaks.narrowPeak")] 
dA[, reads_file_basename:=basename(reads_file)]
dA[, batch := c("A","B","C")[as.numeric(as.factor(flowcell))]]
dA[, sample_name_short := paste_(gsub("^S_","",gsub("_S\\d+$","",sample_name)), batch)]
dA[, condition:=ifelse(grepl("Wt",sample_group), "Wt", "EF")]
dA[, treatment:=gsub("_xeno","",gsub("(Wt|EF)_?","",sample_group))]
dA[treatment=="", treatment:="none"]

colorPalettes <- list() 
colorPalettes <- c(colorPalettes, sapply(setdiff(colnames(dA), names(colorPalettes)), function(annotType) {
  lvls <- dA[,sort(unique(get(annotType)))]
  lvls[is.na(lvls) | lvls==""] <- "NA"
  getColors(lvls)
}, simplify=F))

colorPalettes$FRiP <- RColorBrewer::brewer.pal(9, "RdYlGn")
colorPalettes$pass_qc <- c("FALSE"=colorPalettes$FRiP[1],"TRUE"=colorPalettes$FRiP[length(colorPalettes$FRiP)])
colorPalettes$sample_group <- c(EF_Prx1CreErt2="#dd1c77",Prx1Cre="#222222", structure(RColorBrewer::brewer.pal(5, "YlOrBr")[3:5],names=c("EF","EF_IGF1","EF_IGF1_xeno")), structure(RColorBrewer::brewer.pal(4, "Greys")[3:4],names=c("Wt","Wt_IGF1")) )
colorPalettes$genotype <- c(MSCLC=colorPalettes$sample_group[["Wt"]], "Prx1-MSCLC"=colorPalettes$sample_group[["Prx1Cre"]], "EFPrx1-MSCLC"=colorPalettes$sample_group[["EF"]], "EF Prx1CreErt2"=colorPalettes$sample_group[["EF_Prx1CreErt2"]])
colorPalettes$batch <- getColors(dA[,unique(batch)])
colorPalettes$treatment <- c("none"="white", "IGF1"="red", "Prx1CreErt2"="pink")
colorPalettes$condition <- c("Wt"=colorPalettes$sample_group[["Wt"]], "EF"=colorPalettes$sample_group[["EF"]])
colorPalettes$sig_dir <- c("-1"="red","0"="white","1"="blue")
colorPalettes$module_id <- getColors(paste0("M", 1:5), "Pastel2")

sexGenes <- tolower(sort(unique(c("UTY", "XIST","DDX3Y","KDM5D","ZFY","EIF2S3Y","EIF1AY","LOC110255320","LOC110257894","LOC396706","LOC100625207","LOC110255257"))))

MAX_DIST_GENE <- 100000 
MICROSAT_EF_MOTIF <- "GGAA"	# rev comp = TTCC
MICROSAT_REPEATS <- 8
MICROSAT_MISMATCHES <- 0.1
MICROSAT_RANK_BY <- "stat"
MICROSAT_REGEX <- "agrep_([ATGCN]+)_(\\d+)_(\\d+(.\\d+)?)$"
GENE_ENRICH_PADJ <- 0.005
GENE_MIN_COUNT <- 10
GENE_ENRICH_ODDS <- log2(2)

options("RCACHE.DIR"=file.path(config$out_root,"rcache"))
