#create csv file (chr,position) for every miRNA with 3p and 5p for VIVA that contain all the positions of the specific miRNA (one file per miRNA)

setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data")
#Read miRBase file
ath <- read.table("ath.gff3", sep="\t", quote="")
#filter out all miRNAs without 3p OR 5p
ath <- ath[grepl("-[35]p", ath$V9), ]
#save names of miRNAs
names <- sub(".*Name=([^;]+);.*", "\\1", ath$V9)
names <- sub("-[35]p$", "", names)
names <- sub("ath-", "", names)
#filter out all miRNAs that do not have 3p AND 5p
tab <- table(names)
names <- names(tab[tab == 2])
length(names)
length(unique(names))
#save names as csv file for VIVA (array job)
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA")
write.table(names, file = "names.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

#wd for csv files
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/csv_all")
#iterate over all names and create csv files
for(name in names){
  miR <- ath[grepl(name, ath$V9), ]
  df <- data.frame()
  lower <- miR[which.min(miR$V4), ]
  upper <- miR[which.max(miR$V5), ]
  #-50 and +50 for the base
  for(i in (lower$V4-50):(upper$V5+50)){
    new <- data.frame(chr = substr(lower$V1,4,4), position = i)
    df <- rbind(df, new)
  }
  
  write.table(df, file = paste(name, ".csv", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
}

#I also want to create a csv file with the information about the start position of base_one,miRNA_one,loop,miRNA_two and base_two
#miRNA_one is always the 5p
df_starting_pos <- data.frame()
for(name in names){
  #rows in ath.gff3 for name
  miR <- ath[grepl(pattern = name, x = ath$V9), ]
  #making sure every name was found exactly two times in the ath.gff3
  if(length(miR$V1) != 2){
    print(paste("Error1!!!:",name))
  }
  
  #on the plus strand the 5p is the one with the smaller starting base
  if(miR$V7[1] == "+"){
    miRNA_5p <- miR[which.min(miR$V4), ]
    miRNA_3p <- miR[which.max(miR$V4), ]

    new_starting_pos <- data.frame(name = name, base_one = miRNA_5p$V4-25, miRNA_one = miRNA_5p$V4, loop = miRNA_5p$V5+1, miRNA_two = miRNA_3p$V4, base_two = miRNA_3p$V5+1)
    df_starting_pos <- rbind(df_starting_pos, new_starting_pos)
  }
  #on the minus strand the 5p is the one with the bigger starting base
  else{
    miRNA_5p <- miR[which.max(miR$V4), ]
    miRNA_3p <- miR[which.min(miR$V4), ]
    
    new_starting_pos <- data.frame(name = name, base_one = miRNA_5p$V5+25, miRNA_one = miRNA_5p$V5, loop = miRNA_5p$V4-1, miRNA_two = miRNA_3p$V5, base_two = miRNA_3p$V4-1)
    df_starting_pos <- rbind(df_starting_pos, new_starting_pos)
  }
}

#add the ath-MIR829 (3p.2 and 5p match) and the ath-MIR1886 (.2 and .3 match)
#they were not detected before because of the different names
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data")
#Read miRBase file
ath <- read.table("ath.gff3", sep="\t", quote="")
#rows in ath.gff3 for ath-MIR829 (3p.2 and 5p)
miR <- ath[grepl(pattern = "miR829-3p.2|miR829-5p", x = ath$V9), ]
#ath-MIR829 is on the - strand
miRNA_5p <- miR[which.max(miR$V4), ]
miRNA_3p <- miR[which.min(miR$V4), ]

new_starting_pos <- data.frame(name = "miR829-5p_miR829-3p.2", base_one = miRNA_5p$V5+25, miRNA_one = miRNA_5p$V5, loop = miRNA_5p$V4-1, miRNA_two = miRNA_3p$V5, base_two = miRNA_3p$V4-1)
df_starting_pos <- rbind(df_starting_pos, new_starting_pos)

#csv input for VIVA
df <- data.frame()
#-50 and +50 for the base
for(i in (miRNA_3p$V4-50):(miRNA_5p$V5+50)){
  new <- data.frame(chr = substr(miRNA_3p$V1,4,4), position = i)
  df <- rbind(df, new)
}
write.table(df, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/csv_all/miR829-5p_miR829-3p.2.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") 

#rows in ath.gff3 for ath-MIR1886 (.2 and .3)
miR <- ath[grepl(pattern = "miR1886.2|miR1886.3", x = ath$V9), ]
#ath-MIR1886 is on the + strand
miRNA_5p <- miR[which.min(miR$V4), ]
miRNA_3p <- miR[which.max(miR$V4), ]

new_starting_pos <- data.frame(name = "miR1886.2_miR1886.3", base_one = miRNA_5p$V4-25, miRNA_one = miRNA_5p$V4, loop = miRNA_5p$V5+1, miRNA_two = miRNA_3p$V4, base_two = miRNA_3p$V5+1)
df_starting_pos <- rbind(df_starting_pos, new_starting_pos)

#csv input for VIVA
df <- data.frame()
#-50 and +50 for the base
for(i in (miRNA_5p$V4-50):(miRNA_3p$V5+50)){
  new <- data.frame(chr = substr(miRNA_3p$V1,4,4), position = i)
  df <- rbind(df, new)
}
write.table(df, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/csv_all/miR1886.2_miR1886.3.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

#save file
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps")
write.table(df_starting_pos, file = "regions_miRNAs.csv", quote = FALSE, row.names = FALSE, sep = ",")

