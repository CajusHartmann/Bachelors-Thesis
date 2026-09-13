#TASK: create one heatmap for every miRNA-group (based on processing type) that visualizes the results from all the
#heatmaps from the predicted miRNAs from the specific group
#red regions should indicate that there are lots of SNPs/indels on that specific position
#blue regions should indicate that there are few/no SNPs/indels on that specific position
library(tidyverse)
source("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/functions.R")

#First Step: Group the miRNAs (based on processing type) with the files from Yan et. al, 2024 nature plants (doi: 10.1038/s41477-024-01725-9)
#I already saw for the non-predicted miRNAs that I need those two files (because they contain all miRNAs available)
file_147 <- read.table("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data/processing modes/147.csv", header = TRUE, sep = ";")
file_SE_HYL1 <- read.table("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data/processing modes/SE_HYL1.csv", header = TRUE, sep = ";")

#find the 13 new miRNAs from file_SE_HYL1
new <- file_SE_HYL1[!file_SE_HYL1$name %in% intersect(file_147$name, file_SE_HYL1$name), ]
new$patterns <- new$processing.pattern
#types like "bidirection-BTL", "bidirection-LTB" and "bidirection-SBTL" are replaced with Bidirectional
new$patterns <- sub("bidirection-.*", "Bidirectional", new$patterns)
#create new data frame with all miRNAs from Yan et. al (147 + 13)
complete <- rbind(file_147[ ,c(1,3)], new[ ,c(1,5)])
length(complete[[1]]) #160

#names of the predicted miRNAs 
names_pred <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Create heatmaps/regions_predicted_miRNAs.csv", header = TRUE, sep = ",")[,1]
length(intersect(names_pred, complete[[1]])) #result: 71 --> for 71 of the 273 predicted miRNAs we know the processing type
#names of the miRNAs for which I know 3p and 5p
names <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/regions_miRNAs.csv", header = TRUE, sep = ",")[,1]
#check the miRNAs from Yan et. al, that have not been covered yet (normally every miRNA from the 160 should be found in the non-predicted or in the predicted)
complete[[1]][!(complete[[1]] %in% intersect(names_pred, complete[[1]])) & !(complete[[1]] %in% names)]
#result:
#[1] "miR161"  "miR780a" "miR851a" "miR857a" "miR169L" "miR779"  "miR848a" "miR849a" "miR835a" "miR839a"
#[11] "miR840a" "miR869a" "miR472a" "miR771a" "miR844a" "miR775a" "miR825a" "miR831a" "miR864a" "miR824a"
#[21] "miR827a" "miR868a" "miR776a" "miR779a" "miR822a" "miR823a" "miR828a" "miR867a"

#add 161.1 and 161.2 to complete
new <- data.frame(name = c("miR161.1", "miR161.2"), patterns = "BTL")
complete <- rbind(complete, new)
#add 780.1 and 780.2
new <- data.frame(name = c("miR780.1", "miR780.2"), patterns = "BTL")
complete <- rbind(complete, new)
#add miR857
new <- data.frame(name = "miR857", patterns = "BTL")
complete <- rbind(complete, new)
#add 169l
new <- data.frame(name = "miR169l", patterns = "SBTL")
complete <- rbind(complete, new)
#add 779.1 and 779.2
new <- data.frame(name = c("miR779.1","miR779.2"), patterns = "LTB")
complete <- rbind(complete, new)
#add 848
new <- data.frame(name = "miR848", patterns = "LTB")
complete <- rbind(complete, new)
#add 849
new <- data.frame(name = "miR849", patterns = "LTB")
complete <- rbind(complete, new)
#add 839
new <- data.frame(name = "miR839", patterns = "SLTB")
complete <- rbind(complete, new)
#add 869.1 and 869.2
new <- data.frame(name = c("miR869.1","miR869.2"), patterns = "SLTB")
complete <- rbind(complete, new)
#add 771
new <- data.frame(name = "miR771", patterns = "Bidirectional")
complete <- rbind(complete, new)

#no need to add 775 because it is already there (775a is additional)

#no need to add 825 because it is already there (825a is additional)

#no need to add 827 because it is already there (827a is additional)

#no need to add 776 because it is already there (776a is additional)

#no need to add 779 because it is already there (779a is additional)

#no need to add 822 because it is already there (822a is additional)

#no need to add 823 because it is already there (823a is additional)

#no need to add 828 because it is already there (828a is additional)

#no need to add 867 because it is already there (867a is additional)

length(intersect(names_pred, complete[[1]])) #85

#I need to scale the different regions of the miRNAs to the same length for the heatmaps
#Therefore I need to calculate the median value per processing type and region

#read file with positions of the regions of the miRNAs
positions <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Create heatmaps/regions_predicted_miRNAs.csv", header = TRUE, sep = ",")
#assign the miRNAs to their group
merged <- merge(positions, complete, by.x = "name", by.y = "name")
all <- split(merged, merged$pattern)
write.table(merged, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Create heatmaps/overlap_Yan_and_predicted_miRNAs.csv", sep = ",", row.names = FALSE, quote = FALSE)

#added after creating the file
merged <- read.csv("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Create heatmaps/overlap_Yan_and_predicted_miRNAs.csv")
all <- split(merged, merged$pattern) 

#median values per group and region
medians <- list()
#iterate over all processing types
for(type in names(all)){
  miRNA_one <- c()
  loop <- c()
  miRNA_two <- c()
  #iterate over all miRNAs of the specific type
  for(i in 1:length(all[[type]][[1]])){
    miRNA_one <- c(miRNA_one, abs(all[[type]][i,"loop"] - all[[type]][i,"miRNA_one"]))
    loop <- c(loop, abs(all[[type]][i, "miRNA_two"] - all[[type]][i, "loop"]))
    miRNA_two <- c(miRNA_two, abs(all[[type]][i,"base_two"] - all[[type]][i,"miRNA_two"]))
  }
  medians[[type]]["miRNA_one"] <- median(miRNA_one)
  medians[[type]]["loop"] <- median(loop)
  medians[[type]]["miRNA_two"] <-  median(miRNA_two)
}
medians$Bidirectional
medians$BTL
medians$LTB
medians$SBTL
medians$SLTB
#median for all miRNA_one and miRNA_two is 21 (except for SBTL miRNA_two but I will set it to 21, so that is is consistent)
medians$SBTL["miRNA_two"] <- 21
medians$SLTB["loop"] <- 132


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
  matrix <- read.table(file = paste("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/VIVA/results_all/",merged[i,1], "/", "Genotype_miRNA_positions_pred.vcf_matrix.csv", sep = ""), header = TRUE, sep = ",")  
  matrix <- matrix %>% separate(chr.position, into = c("chr", "pos"), sep = ",", convert = TRUE)
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
    mat_loop_scaled <- scaler(mat_loop, medians[[merged[i, "patterns"]]][["loop"]])
    #combine the scaled regions (and the un-scaled regions) to one matrix
    matrix_scaled <- cbind(mat_base_one, mat_miRNA_one, mat_loop_scaled, mat_miRNA_two, mat_base_two)
  }
  #when there is something in the middle to be scaled
  else{
    mat_loop_scaled <- scaler(mat_loop[ , 21:(ncol(mat_loop)-20)], medians[[merged[i, "patterns"]]][["loop"]]-40)
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
  miRNA_two <- loop + medians[[name]][["loop"]]
  base_two <- miRNA_two + 21
  
  #call function to create a heatmap
  create_heatmap(combined_mats[[name]], base_one, miRNA_one, loop, miRNA_two, base_two, 0, 1, name, TRUE, colorscale = "RdBu", path = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Create heatmaps/heatmaps_grouped_group_size")
} 
