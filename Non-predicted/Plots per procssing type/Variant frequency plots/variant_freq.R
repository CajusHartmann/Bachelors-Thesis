#TASK: create one plot per group, that shows the SNP/indel frequency on the y-axis 
#and the position of the miRNA on the x-axis (scaling of the regions the same as in recreation of aggregated heatmaps)
#5p miRNA is always miRNA_one
library(tidyverse)
source("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/functions.R")

#read file with overlaps of the 160 miRNAs from Yan and the 95 for which we know 3p and 5p
merged <- read.csv("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/overlap_Yan_and_95_5pFirst.csv")
all <- split(merged, merged$pattern)

#medians of the loops
medians <- list()
medians$Bidirectional <- 61; medians$BTL <- 51; medians$LTB <- 42; medians$SBTL <- 62; medians$SLTB <- 114

all_mats <- list()
#lists to create the tables 
numbers <- list()
#list of data frames that will be used for statistical tests
stat <- list(base_one = data.frame(), miRNA_one = data.frame(), loop = data.frame(), miRNA_two = data.frame(), base_two = data.frame())
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
  #convert all 1 to 0 and all 3 to 1
  matrix_all_positions[matrix_all_positions == 1] <- 0
  matrix_all_positions[matrix_all_positions == 2] <- 1
  matrix_all_positions[matrix_all_positions == 3] <- 1
  #extract the regions from the matrix
  mat_base_one <- matrix_all_positions[ ,1:25]
  mat_miRNA_one <- matrix_all_positions[ ,26:(loop-1)]
  mat_loop <- matrix_all_positions[ ,loop:(miRNA_two-1)]
  mat_miRNA_two <- matrix_all_positions[ ,miRNA_two:(base_two-1)]
  mat_base_two <- matrix_all_positions[ ,base_two:(base_two + 24)]
  matrix_all_positions <- t(matrix_all_positions)

  #create table that shows the processing type, the region and the SNP-freq in percent (and another table with numbers, not percentages)
  #therefore I look for every miRNA of the type separately where the regions start
  #here I calculate the total number of SNPs of a group and region
  numbers[[merged[i, "patterns"]]][["base_one"]] <- sum(mat_base_one, numbers[[merged[i, "patterns"]]][["base_one"]], na.rm = TRUE)
  numbers[[merged[i, "patterns"]]][["miRNA_one"]] <- sum(mat_miRNA_one, numbers[[merged[i, "patterns"]]][["miRNA_one"]], na.rm = TRUE)
  numbers[[merged[i, "patterns"]]][["loop"]] <- sum(mat_loop, numbers[[merged[i, "patterns"]]][["loop"]], na.rm = TRUE)
  numbers[[merged[i, "patterns"]]][["miRNA_two"]] <- sum(mat_miRNA_two, numbers[[merged[i, "patterns"]]][["miRNA_two"]], na.rm = TRUE)
  numbers[[merged[i, "patterns"]]][["base_two"]] <- sum(mat_base_two, numbers[[merged[i, "patterns"]]][["base_two"]], na.rm = TRUE)
  #fill list of data frames for statistical tests
  for(reg in c("base_one", "miRNA_one", "loop", "miRNA_two", "base_two")){
    len <- switch(reg,
                  "base_one" = 25,
                  "miRNA_one" = abs(merged[i, "miRNA_one"] - merged[i, "loop"]),
                  "loop" = abs(merged[i, "miRNA_two"] - merged[i, "loop"]),
                  "miRNA_two" = abs(merged[i, "base_two"] - merged[i, "miRNA_two"]),
                  "base_two" = 25)
    
    new <- data.frame(name = merged[i, "name"], region = reg, num_SNPs = sum(get(paste0("mat_", reg)), na.rm = TRUE), type = merged[i, "patterns"], length = len)
    stat[[reg]] <- rbind(stat[[reg]], new)
  }
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
  
  #for the line plot I need the SNP/indel frequency at each position
  SNP_freq <- colSums(matrix_scaled, na.rm = TRUE) / 1135
  #save SNP_freq
  all_mats[[merged[i,"patterns"]]][[merged[i,"name"]]] <- SNP_freq
}

#create file for statistical tests
write.table(stat$base_one, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Other plots per group/statistical_analysis_non_predicted.csv", sep = ",", quote = FALSE, row.names = FALSE)
write.table(stat$miRNA_one, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Other plots per group/statistical_analysis_non_predicted.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
write.table(stat$loop, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Other plots per group/statistical_analysis_non_predicted.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
write.table(stat$miRNA_two, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Other plots per group/statistical_analysis_non_predicted.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
write.table(stat$base_two, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Other plots per group/statistical_analysis_non_predicted.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

#create tables with the help of the list named numbers
df_numbers <- t(data.frame(LTB = unlist(numbers$LTB), BTL = unlist(numbers$BTL), Bidirectional = unlist(numbers$Bidirectional), SBTL = unlist(numbers$SBTL), SLTB = unlist(numbers$SLTB)))
df_percentages <- t(data.frame(LTB = unlist(numbers$LTB) / sum(unlist(numbers$LTB)) * 100, BTL = unlist(numbers$BTL) / sum(unlist(numbers$BTL)) * 100, Bidirectional = unlist(numbers$Bidirectional) / sum(unlist(numbers$Bidirectional)) * 100, SBTL = unlist(numbers$SBTL) / sum(unlist(numbers$SBTL)) * 100, SLTB = unlist(numbers$SLTB) / sum(unlist(numbers$SLTB)) * 100))
write.table(round(df_percentages,2), file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Other plots per group/recreate_recreated_percentages.csv", sep = ",", quote = FALSE)
write.table(df_numbers, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Other plots per group/recreate_recreated_numbers.csv", sep = ",", quote = FALSE)

#create plots
#iterate over all processing types
for(name in names(all_mats)){
  png(paste("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Other plots per group/recreate_recreated_SNP-freq_pos/", name, ".png", sep=""), width = 14, height = 12, units = "cm", res=300)
  
  par(tcl = -0.2, mgp = c(2, 0.3, 0), cex.axis = 0.5, cex.lab = 0.8)
  my_col <- adjustcolor("steelblue", alpha.f = 0.3)
  #define regions for x-axis
  base_one <- 1; miRNA_one <- 26; loop <- 26 + 21; miRNA_two <- loop + medians[[name]]; base_two <- miRNA_two + 21
  regions <- data.frame(start = c(base_one, miRNA_one, loop, miRNA_two, base_two), end = c(miRNA_one-1, loop-1, miRNA_two-1, base_two-1, base_two+24), name = c("base", "miRNA", "loop", "miRNA", "base"))
  
  first_el <- all_mats[[name]][[names(all_mats[[name]])[1]]]
  x <- 1:length(first_el)
  #plot first element of the processing type
  plot(x, first_el, xaxt = "n", ylab = "SNP/indel frequency", xlab = "", type = "l", ylim = c(0,1), col = my_col)
  #add all other elements to the plot
  for(name_miRNA in names(all_mats[[name]])[2:length(all_mats[[name]])]){
    points(x, all_mats[[name]][[name_miRNA]], col = my_col, type = "l")
  }
  #add upper axis with labels
  axis(side = 3, at = seq(0, length(first_el), by = 40))
  mtext("Position", side = 3, line = 1.5, cex = 0.8)
  #add lower axis with regions
  axis(side = 1, at = (regions$start + regions$end) / 2, labels = regions$name, tick = FALSE, line = -0.5)
  axis(side = 1, at = c(regions$start, max(regions$end)), labels = FALSE)
  mtext("Region", side = 1, line = 1, cex = 0.8)
  
  #create legend
  #function to calculate transparency, when 'x' lines with transparency 'a' are overlapping
  transparency <- function(x, a){1 - (1 - a)^x}
  #calculate colors for the legend
  x_vals <- c(1, 2, 4, 8, 16)
  alphas <- sapply(x_vals, function(x){transparency(x, 0.3)})
  colours <- sapply(alphas, function(a){adjustcolor("steelblue", alpha.f = a)})
  
  legend(x = "topright", legend = c("1", "2", "4", "8", "16"), col = colours, lwd = 1.5, bty = "n", cex = 0.5)
  dev.off()
}