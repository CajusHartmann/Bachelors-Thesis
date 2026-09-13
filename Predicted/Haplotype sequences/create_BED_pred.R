#create bed file with CHR START END NAME of each miRNA (predicted)

reg <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Create heatmaps/regions_predicted_miRNAs.csv", header = TRUE, sep = ",")
bed <- data.frame()

for(name in reg[ ,1]){
  #extract chromosome
  pos <- read.table(file = paste0("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/VIVA/csv_all_predicted/", name, ".csv"), header = FALSE, sep = ",")
  #check if the miRNA is on the - or + strand
  if(reg[reg$name == name,2] > reg[reg$name == name,3]){
    #Attention: bed file is 0-based --> start = start - 1
    new <- data.frame(CHR = pos[1,1], START = pos[1,2]-1, END = pos[length(pos[[1]]),2], NAME = name, STRAND = "-")
    bed <- rbind(bed, new)
  }
  else{
    new <- data.frame(CHR = pos[1,1], START = pos[1,2]-1, END = pos[length(pos[[1]]),2], NAME = name, STRAND = "+")
    bed <- rbind(bed, new)
  }
  
}

write.table(bed, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Haplotype sequences/coordinates_pred.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
