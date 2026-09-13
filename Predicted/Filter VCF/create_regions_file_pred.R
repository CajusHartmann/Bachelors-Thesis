setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Create heatmaps")

reg <- read.csv("regions_predicted_miRNAs.csv")

df <- data.frame()

setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/VIVA/csv_all_predicted")

for(name in reg[ ,1]){
  chr_csv <- read.csv(paste0(name, ".csv"))
  start_pos <- min(reg[reg$name==name, 2], reg[reg$name==name, 6])
  end_pos <- max(reg[reg$name==name, 2], reg[reg$name==name, 6])
  
  if(start_pos == reg[reg$name==name, 2]){
    new <- data.frame(chr=chr_csv[1,1],start=start_pos, end=end_pos+25)
  }
  else{
    new <- data.frame(chr=chr_csv[1,1],start=start_pos-25, end=end_pos)
  }
 
  df <- rbind(df,new)
}

setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Filter VCF")
write.table(df, file = "regions_pred.tsv", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)


