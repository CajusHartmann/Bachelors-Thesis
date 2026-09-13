#TASK: calculate the median length of the different regions with the new overlap_Yan_and_95_5pFirst.csv file from the 1.7.2026

#read file with overlaps of the 160 miRNAs from Yan and the 95 for which we know 3p and 5p
merged <- read.csv("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/overlap_Yan_and_95_5pFirst.csv")
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
