#TASK: use the genotype matrices created by viva to create heatmaps of the whole miRNAs (for every miRNA with 3p and 5p)
library(tidyverse)
source("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/functions.R")

positions <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/regions_miRNAs.csv", header = TRUE, sep = ",")
#iterate over all miRNAs
for(i in 1:length(positions$name)){
  #start positions of the different regions of the miRNA (when the miRNA starts at 1 and ends at length(miRNA))
  #we set the length of the base to 25bp
  base_one <- 1
  miRNA_one <- abs(positions[i, "miRNA_one"] - positions[i, "base_one"]) + 1
  loop <- abs(positions[i, "loop"] - positions[i, "base_one"]) + 1
  miRNA_two <- abs(positions[i, "miRNA_two"] - positions[i, "base_one"]) + 1
  base_two <- abs(positions[i, "base_two"] - positions[i, "base_one"]) + 1
  
  #read genotype matrix created by VIVA
  matrix <- read.table(file = paste("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/results_all/all/", positions[i,1], "/Genotype_miRNA_positions.vcf_matrix.csv", sep = ""), header = TRUE, sep = ",")  
  matrix <- matrix %>% separate(chr.position, into = c("chr", "pos"), sep = ",", convert = TRUE)
  #filter out the positions outside of the miRNA (originally the base was set to 50bp, now it is 25)
  valid_pos <- read.table(file = paste0("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/csv_all/", positions[i,1], ".csv"), header = FALSE, sep = ",")
  #+25 and -25 because originally viva was started with the length of the base set to 50bp
  matrix <- matrix[matrix$pos >= valid_pos[1,2]+25 & matrix$pos <= valid_pos[length(valid_pos[ ,1]),2]-25, ]
  #convert positions
  matrix$pos <- abs(matrix$pos - positions[i, "base_one"]) + 1
  #create new matrix with all positions of the miRNA (rows) and all accessions (columns)
  matrix_all_positions <- matrix(data = NA, nrow = base_two + 24, ncol = ncol(matrix) - 2)
  #fill new matrix with the values of the genotype matrix
  for(x in 1:nrow(matrix)){
    matrix_all_positions[matrix$pos[x], ] <- as.numeric(matrix[x,3:ncol(matrix)])
  }
  
  name <- sub("mi", "MI", positions[i,1])
  #call function to create a heatmap
  create_heatmap(t(matrix_all_positions), base_one, miRNA_one, loop, miRNA_two, base_two, 0, 3, name, FALSE, colorscale = list(list(0.0, "white"), list(1/3, "lightblue"), list(2/3, "green"), list(1.0, "red")), path = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/heatmaps")
}


