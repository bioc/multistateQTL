# Subset GTEx data
# Amelia Dunstone
# 2024-10-07

# For GTEx, I downloaded the "GTEx_Analysis_v8_eQTL.tar" data from 
# https://www.gtexportal.org/home/downloads/adult-gtex/qtl.
# I unzipped the tar file using 7zip. 

root_dir <- "GTEx/data/GTEx_Analysis_v8_eQTL/"

# All file paths in the directory
filenames <- list.files(root_dir)

# State names
states <- sapply(strsplit(filenames, "\\."), "[[", 1)

# States of interest
selected_states <- c("Muscle_Skeletal", "Whole_Blood", "Skin_Sun_Exposed_Lower_leg", 
                     "Adipose_Subcutaneous", "Thyroid", "Lung", "Cells_Cultured_fibroblasts", 
                     "Stomach", "Pancreas", "Pituitary", "Spleen", "Liver", "Brain_Cerebellum", "Brain_Cortex")

# Short names for these states?
short_names <- c("muscle", "blood", "skin", "adipose", "thyroid", "lung", "fibroblast",
                 "stomach", "pancreas", "pituitary", "spleen", "liver", "brain (cerebellum)",
                 "brain (cortex)")

# Selected file paths
files <- filenames[states %in% selected_states & grepl("signif", filenames)]

input <- data.frame(
    filenames = paste0(root_dir, selected_states, ".v8.egenes.txt.gz"),
    destnames = paste0("GTEx/data/", selected_states, ".v8.egenes.txt"), 
    states = selected_states)

library(R.utils)

# Loop through and extract all the files
for (i in 1:nrow(input)) {
    gunzip(input$filenames[i], destname = input$destnames[i], remove = FALSE) 
}


brain <- read.delim("GTEx/data/Brain_Cerebellum.v8.egenes.txt", sep = "\t")

head(brain)

# We have the columns: gene_id, gene_chr, variant_id, slope, slope_se, qval

# Q: Which p value?

#             gene_id         gene_name gene_chr gene_start gene_end strand num_var beta_shape1 beta_shape2
# 1 ENSG00000227232.5            WASH7P     chr1      14410    29553      -    1364     1.02506     283.802
# 2 ENSG00000268903.1     RP11-34P13.15     chr1     135141   135895      -    1863     1.05390     304.825
# 3 ENSG00000269981.1     RP11-34P13.16     chr1     137682   137965      -    1868     1.03909     326.232
# 4 ENSG00000241860.6     RP11-34P13.13     chr1     141474   173862      -    2066     1.04385     341.834
# 5 ENSG00000279928.2 ABC7-43046700E7.1     chr1     182696   184174      +    2135     1.04488     335.959
# 6 ENSG00000279457.4     RP11-34P13.18     chr1     185217   195411      -    2234     1.04828     363.268
#   true_df pval_true_df            variant_id tss_distance  chr variant_pos ref alt num_alt_per_site
# 1 153.767  8.91297e-04   chr1_985691_G_A_b38       956138 chr1      985691   G   A                1
# 2 149.080  9.83353e-05   chr1_108826_G_C_b38       -27069 chr1      108826   G   C                1
# 3 151.968  1.50509e-05   chr1_108826_G_C_b38       -29139 chr1      108826   G   C                1
# 4 150.697  2.18732e-04    chr1_14677_G_A_b38      -159185 chr1       14677   G   A                1
# 5 147.892  3.29576e-03 chr1_1081644_C_CG_b38       898948 chr1     1081644   C  CG                1
# 6 149.789  4.35221e-04   chr1_602064_G_A_b38       406653 chr1      602064   G   A                1
#   rs_id_dbSNP151_GRCh38p7 minor_allele_samples minor_allele_count       maf ref_factor pval_nominal    slope
# 1              rs78164078                   25                 26 0.0622010          1  4.94027e-04 0.530302
# 2              rs62642117                   15                 15 0.0358852          1  3.36089e-05 1.119820
# 3              rs62642117                   15                 15 0.0358852          1  5.00013e-06 1.237960
# 4             rs201327123                   23                 23 0.0550239          1  9.03769e-05 0.919054
# 5             rs113385670                    9                 10 0.0273224          1  1.67739e-03 0.867407
# 6             rs375649215                    7                  7 0.0167464          1  1.86060e-04 1.477630
#   slope_se pval_perm  pval_beta       qval pval_nominal_threshold log2_aFC log2_aFC_lower log2_aFC_upper
# 1 0.149260 0.2140260 0.21336300 0.12744800            0.000243956 0.683372       0.470040       0.877250
# 2 0.262750 0.0210979 0.02386460 0.02262380            0.000248231 2.087923       1.619353       2.486133
# 3 0.262474 0.0039996 0.00391188 0.00459953            0.000221721 2.419963       2.062900       3.005948
# 4 0.229093 0.0611939 0.06303990 0.05069870            0.000214717 2.275564       1.750418       2.815456
# 5 0.271626 0.6525420 0.65184900 0.26981500            0.000219160 1.982902       0.633628       2.901933
# 6 0.386648 0.1257220 0.13076600 0.08925470            0.000204795 0.779340       0.386629       1.035745


