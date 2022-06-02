# Load libraries
library(immunedeconv)
library(tidyverse)

# Get snakemake input and outputs
base_dir <- snakemake@input[['base_dir']]
infile   <- snakemake@input[['infile']]
method   <- snakemake@params[['method']]
outfile  <- snakemake@output[['outfile']]

set_cibersort_binary(paste0(base_dir, "src/CIBERSORT.R"))
set_cibersort_mat(paste0(base_dir, "src/LM22.txt"))

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