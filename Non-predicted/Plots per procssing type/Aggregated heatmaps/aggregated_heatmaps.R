#TASK: recreate the recreated aggregated heatmaps with the length of the bases set to 25 and
#without scaling the important regions of the loop (20bp at the start and 20bp at the end)
#and with the 5p being always miRNA_one
library(tidyverse)
source("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/functions.R")

####################################################################
####################################################################
####################################################################
#This part was originally in another script:

#First Step: Group the miRNAs (based on processing type) with the file from Yan etal, 2024 nature plants (doi: 10.1038/s41477-024-01725-9)
#Therefore I need to check which table to use from the file (I check if the file with the 147 DCL1 dependent miRNAs contains all miRNAs that are in the other files)
file_147 <- read.table("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data/processing modes/147.csv", header = TRUE, sep = ";")
file_52 <- read.table("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data/processing modes/52.csv", header = TRUE, sep = ";")
file_43 <- read.table("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data/processing modes/43.csv", header = TRUE, sep = ";")
file_SE_HYL1 <- read.table("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data/processing modes/SE_HYL1.csv", header = TRUE, sep = ";")

length(intersect(file_147$name, file_52$name)) #result: 52 --> no new miRNAs in this file
length(intersect(file_147$name, file_43$name)) #result: 43 --> no new miRNAs in this file
length(intersect(file_147$name, file_SE_HYL1$name)) #result: 103 --> 13 new miRNAs in this file because it contains 116 miRNAs
#find the 13 new miRNAs from this file
new <- file_SE_HYL1[!file_SE_HYL1$name %in% intersect(file_147$name, file_SE_HYL1$name), ]
new$patterns <- new$processing.pattern
#types like "bidirection-BTL", "bidirection-LTB" and "bidirection-SBTL" are replaced with Bidirectional
new$patterns <- sub("bidirection-.*", "Bidirectional", new$patterns)
#create new data frame with all miRNAs from Yan et. al (147 + 13)
complete <- rbind(file_147[ ,c(1,3)], new[ ,c(1,5)])
length(complete[[1]])

#names of the miRNAs for which I know 3p and 5p
names <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/regions_miRNAs.csv", header = TRUE, sep = ",")[,1]
length(intersect(names, complete[[1]])) #result: 61 --> for 61 of the 95 miRNAs with 3p and 5p we know the processing type

#read file with positions of the regions of the miRNAs
positions <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/regions_miRNAs.csv", header = TRUE, sep = ",")
#assign the miRNAs to their group
merged <- merge(positions, complete, by.x = "name", by.y = "name")
all <- split(merged, merged$pattern)
length(all$SLTB[[1]])
write.table(merged, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/overlap_Yan_and_93.csv", sep = ",", row.names = FALSE, quote = FALSE)
####################################################################
####################################################################
####################################################################

#read file with overlaps of the 160 miRNAs from Yan and the 95 for which we know 3p and 5p
merged <- read.csv("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/overlap_Yan_and_95_5pFirst.csv")
all <- split(merged, merged$pattern)

#medians of the loops
medians <- list()
medians$Bidirectional <- 61; medians$BTL <- 51; medians$LTB <- 42; medians$SBTL <- 62; medians$SLTB <- 114

all_mats <- list()
#create a matrix for every miRNA with the scaled regions
for(i in 1:length(merged$name)){
  #start positions of the different regions of the miRNA (when the miRNA starts at 1 and ends at length(miRNA))
  #we set the length of the base to 25bp
  base_one <- 1
  miRNA_one <- 26
  loop <- abs(merged[i, "loop"] - merged[i, "base_one"]) + 1
  miRNA_two <- abs(merged[i, "miRNA_two"] - merged[i, "base_one"]) + 1
  base_two <- abs(merged[i, "base_two"] - merged[i, "base_one"]) + 1
  
  #read genotype matrix created by VIVA
  matrix <- read.table(file = paste("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/results_all/all/",merged[i,1], "/", "Genotype_miRNA_positions.vcf_matrix.csv", sep = ""), header = TRUE, sep = ",")  
  matrix <- matrix %>% separate(chr.position, into = c("chr", "pos"), sep = ",", convert = TRUE)
  #filter out the positions outside of the miRNA (originally the base was set to 50bp, now it is 25)
  valid_pos <- read.table(file = paste0("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/csv_all/", merged[i,1], ".csv"), header = FALSE, sep = ",")
  #+25 and -25 because originally viva was started with the length of the base set to 50bp
  matrix <- matrix[matrix$pos >= valid_pos[1,2]+25 & matrix$pos <= valid_pos[length(valid_pos[ ,1]),2]-25, ]
  #convert positions
  matrix$pos <- abs(matrix$pos - merged[i, "base_one"]) + 1
  #create new matrix with all positions of the miRNA (rows) and all accessions (columns)
  matrix_all_positions <- matrix(data = NA, nrow = base_two + 24, ncol = ncol(matrix) - 2)
  #fill new matrix with the values of the genotype matrix
  for(x in 1:nrow(matrix)){
      matrix_all_positions[matrix$pos[x], ] <- as.numeric(matrix[x,3:ncol(matrix)])
  }
  matrix_all_positions <- t(matrix_all_positions)
  #extract the regions from the matrix
  mat_base_one <- matrix_all_positions[ ,1:25]
  mat_miRNA_one <- matrix_all_positions[ ,26:(loop-1)]
  mat_loop <- matrix_all_positions[ ,loop:(miRNA_two-1)]
  mat_miRNA_two <- matrix_all_positions[ ,miRNA_two:(base_two-1)]
  mat_base_two <- matrix_all_positions[ ,base_two:(base_two + 24)]
  
  #scale the regions
  mat_miRNA_one <- scaler(mat_miRNA_one, 21)
  mat_miRNA_two <- scaler(mat_miRNA_two, 21)
  #only scale the middle of the loop (don't scale the first and last 20bp) to median(loop) - 40
  #when there is nothing in the middle to be scaled --> scale the whole loop to the median
  if(ncol(mat_loop) < 42){
    mat_loop_scaled <- scaler(mat_loop, medians[[merged[i, "patterns"]]])
    #combine the scaled regions (and the un-scaled regions) to one matrix
    matrix_scaled <- cbind(mat_base_one, mat_miRNA_one, mat_loop_scaled, mat_miRNA_two, mat_base_two)
  }
  #when there is something in the middle to be scaled
  else{
    mat_loop_scaled <- scaler(mat_loop[ , 21:(ncol(mat_loop)-20)], medians[[merged[i, "patterns"]]]-40)
    #combine the scaled regions (and the un-scaled regions) to one matrix
    matrix_scaled <- cbind(mat_base_one, mat_miRNA_one, mat_loop[ ,1:20], mat_loop_scaled, mat_loop[ ,(ncol(mat_loop)-20+1):ncol(mat_loop)], mat_miRNA_two, mat_base_two)
  }
  #save matrix
  all_mats[[merged[i,"patterns"]]][[merged[i,"name"]]] <- matrix_scaled
}

#add up the matrices of the groups (so that I receive one matrix for every group that I can use for the heatmap)
#0.0 -> no variants (either no data or reference), 0.25 -> 25% of the miRNAs have a variant, 0.5 -> 50% of the...
combined_mats <- list()
for(name in names(all_mats)){
  #convert all 1 to 0 and all 3 to 1
  mats <- lapply(all_mats[[name]],function(m){m[m==0] <- NA; m[m==1] <- 0; m[m==2] <- 1; m[m==3] <- 1; m})
  #sum over all matrices (NA + NA = NA, NA + number = number)
  sum_mat <- Reduce(function(x,y){all_na <- is.na(x) & is.na(y); x[is.na(x)] <- 0; y[is.na(y)] <- 0; out <- x+y; out[all_na] <- NA; out}, mats)
  #save sum/group size
  combined_mats[[name]] <- sum_mat / length(mats)
  print(length(mats))
  
  #create heatmap
  base_one <- 1
  miRNA_one <- 26
  loop <- 26 + 21
  miRNA_two <- loop + medians[[name]]
  base_two <- miRNA_two + 21
  
  #call function to create a heatmap
  create_heatmap(combined_mats[[name]], base_one, miRNA_one, loop, miRNA_two, base_two, 0, 1, name, TRUE, colorscale = "RdBu", path = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/recreate_recreated_heatmaps_grouped_group_size")
} 