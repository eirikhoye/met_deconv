# Load libraries
library(immunedeconv)
library(tidyverse)

# Get snakemake input and outputs
base_dir <- '/Users/eirikhoy/Dropbox/projects/met_deconv'
infile   <- 'data/count_matrix_normalScale.xlsx'

method   <- 'quantiseq'

outfile  <- 'results/count_matrix_normalscale_qts.tsv'

set_cibersort_binary(paste(base_dir, "src/CIBERSORT.R", sep=""))
set_cibersort_mat(paste(base_dir, "src/LM22.txt", sep=""))

# Load data
setwd(base_dir)
gene_expression_matrix <- readxl::read_excel(paste0(base_dir, '/', infile))
gene_expression_matrix <- gene_expression_matrix %>% select(!ENSG)
gene_expression_matrix <- as.data.frame(gene_expression_matrix)
rownames(gene_expression_matrix) <- gene_expression_matrix$Gene
gene_expression_matrix$Gene <- NULL

print(gene_expression_matrix)

# define function
run_deconvolute <- function(gene_expression_matrix, method, outfile){
    res <- immunedeconv::deconvolute(gene_expression_matrix, method)
    readr::write_tsv(res, file = paste0(base_dir, '/', outfile))
}

run_deconvolute(gene_expression_matrix, method, outfile)