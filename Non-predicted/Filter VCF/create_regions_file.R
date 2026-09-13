#TASK: create tsv file with the positions of the miRNA (including loop and 50bp upstream and downstream) for bfctools (-R argument)
library(dplyr)
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data")
#Read miRBase file
ath <- read.table("ath.gff3", sep="\t", quote="")
#only keep miRNA
ath <- ath[ath$V3 == "miRNA",]
paste("miRBase: ", length(ath[[1]]), sep="")
#Add column with the name of the miRNA (from the 9th column) for the grouping (see below)
ath$name <- sub(".*Name=([^;]+);.*", "\\1", ath$V9)
ath$name <- sub("-[35]p$", "", ath$name)
#Group the miRNAs by name and add start and endposition (min and max)
result <- ath %>% group_by(name) %>% summarise(chr = first(V1), start = min(V4)-50, end = max(V5)+50, .groups = "drop")
result <- result[,2:4]
result$chr <- substr(result$chr,4,4)

#create tsv from result
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/bcftools")
write.table(result, file = "regions.tsv", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)




