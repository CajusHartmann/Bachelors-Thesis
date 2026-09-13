#TASK: filter the csv files (with the information about the genotypes) created by viva and only keep 
#the variants with MAF (minor allel frequency) >= 0.005 and missing rate <= 0.5
library(tidyverse)

names <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/regions_miRNAs.csv", header = TRUE, sep = ",")[,1]
changes <- data.frame()

for(name in names){
  #read genotype matrix created by VIVA
  matrix <- read.table(file = paste("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/results_all/all/", name, "/Genotype_miRNA_positions.vcf_matrix.csv", sep = ""), header = TRUE, sep = ",") 
  matrix <- matrix %>% separate(chr.position, into = c("chr", "pos"), sep = ",", convert = TRUE)
  #check for every variant (row) if MAF >= 0.005 and missing rate <= 0.5
  #meaning of the numbers in matrix: 0 = no data, 1 = homozygous reference, 2 = heterozygous variant, 3 = homozygous variant
  filtered <- matrix[rowMeans(matrix[ ,c(-1,-2)] == 3) >= 0.005 & rowMeans(matrix[ ,c(-1,-2)] == 0) <= 0.5, ]
  #save new matrix
  write.table(filtered, file = paste("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/results_all/filtered/", name, "_Genotype_MAF0.005_MR0.5_vcf.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE)
  #save information about the changes
  new <- data.frame(name = name, variants_before_filtering = length(matrix[[1]]), variants_after_filtering = length(filtered[[1]]), num_of_removed_variants = length(matrix[[1]]) - length(filtered[[1]]))
  changes <- rbind(changes, new)
}

write.table(changes, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/results_all/filtered/changes.csv", quote = FALSE, sep = ",", row.names = FALSE)

#see how the different parameters (set per hand) influence the filtering
MAF0.005_MR0.5 <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/results_all/changes_MAF0.005_MR0.5.csv", header = TRUE, sep = ",")
MAF0.005_MR0.8 <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/results_all/changes_MAF0.005_MR0.8.csv", header = TRUE, sep = ",")
MAF0.01_MR0.5 <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/results_all/changes_MAF0.01_MR0.5.csv", header = TRUE, sep = ",")
MAF0.01_MR0.8 <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/results_all/changes_MAF0.01_MR0.8.csv", header = TRUE, sep = ",")

means_MAF0.005_MR0.8 <- colMeans(MAF0.005_MR0.8[ ,c(-1)])
means_MAF0.005_MR0.5 <- colMeans(MAF0.005_MR0.5[ ,c(-1)])
means_MAF0.01_MR0.8 <- colMeans(MAF0.01_MR0.8[ ,c(-1)])
means_MAF0.01_MR0.5 <- colMeans(MAF0.01_MR0.5[ ,c(-1)])


