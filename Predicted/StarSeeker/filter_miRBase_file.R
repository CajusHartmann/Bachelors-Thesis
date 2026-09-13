#for the use of StarSeeker two input files are required (precursor sequences and mature sequences)
#to get the relevant miRNAs (the ones without 3p AND 5p) from A. thaliana I have to filter the hairpin.fa and mature.fa files
library(seqinr)
#read files from miRBase
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data/predict_missing")
hairpin <- read.fasta(file = "hairpin.fa", as.string = TRUE, whole.header = TRUE, forceDNAtolower = FALSE)
mature <- read.fasta(file = "mature.fa", as.string = TRUE, whole.header = TRUE, forceDNAtolower = FALSE)

#filter for Arabidopsis thaliana
idx_h <- grepl(pattern = ".*Arabidopsis thaliana.*", x = names(hairpin))
hairpin <- hairpin[idx_h]
names(hairpin)

idx_m <- grepl(pattern = ".*Arabidopsis thaliana.*", x = names(mature))
mature <- mature[idx_m]
names(mature)

#read file with the names of the miRNAs for which we know 3p AND 5p (these ones will be filtered out)
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA")
names <- read.table(file = "names.csv", sep = ",")
#create regular expression
names <- paste0(".*", names[,1], ".*", collapse = "|")

#filter out the miRNAs where we already know 3p AND 5p
idx_name_m <- grepl(pattern= names, x = names(mature))
mature <- mature[!idx_name_m]
length(names(mature)) #242

names <- gsub("mi", "MI", names) #change mi to MI because the miRNAs are written as MIR156e in the hairpin file
idx_name_h <- grepl(pattern = names, x = names(hairpin))
hairpin <- hairpin[!idx_name_h]
length(names(hairpin))  #233

#save the files
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data/predict_missing")
write.fasta(sequences = mature, gsub("ath-","", gsub(" .*", "", names(mature))), file.out = "mature_filtered.fa")
write.fasta(sequences = hairpin, gsub("ath-","", gsub(" .*", "", names(hairpin))), file.out = "hairpin_filtered.fa")



