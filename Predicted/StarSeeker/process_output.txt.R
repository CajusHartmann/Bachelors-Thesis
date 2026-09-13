#The output from Starseeker is a file with the sequences of the miRNA*s in the precursors
#For my analysis I need the positions of the miRNA*s (so that I can give VIVA the positions of the precursors)

library(seqinr)
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/StarSeeker")
#read output from Starseeker
output <- read.fasta("output.txt", as.string = TRUE, whole.header = TRUE, forceDNAtolower = FALSE)
#read fasta with hairpin structures
hairpin <- read.fasta(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data/predict_missing/hairpin_filtered.fa", as.string = TRUE, whole.header = TRUE, forceDNAtolower = FALSE)

length(output) #237 (8 less than the original output because I removed a sequence with an error and 7 sequences that are already in miRBase or not that good like the other prediction that was made)
length(hairpin) #233

#calculate the positions of the mature sequences in the precursors (hairpin)
positions_in_hairpin <- list()
num_dup <- list()
sequence <- list()
for(i in 1:length(output)){
  if(is.null(positions_in_hairpin[[names(output)[i]]])){
    positions_in_hairpin[[names(output)[i]]] <- regexpr(pattern = output[[i]][1], text = hairpin[[sub(" \\*star\\*", "", names(output)[i])]][1])[1]
    sequence[[names(output)[i]]] <- output[[i]][1]
    num_dup[[names(output)[i]]] <- 1
  }
  else{
    positions_in_hairpin[[paste0(names(output)[i], ".", num_dup[names(output)[i]])]] <- regexpr(pattern = output[[i]][1], text = hairpin[[sub(" \\*star\\*", "", names(output)[i])]][1])[1]
    sequence[[paste0(names(output)[i], ".", num_dup[names(output)[i]])]] <- output[[i]][1]
    num_dup[[names(output)[i]]] <- num_dup[[names(output)[i]]] + 1
  }
}

#check which MIR ocures twice in the output.txt (for some it is okay to occure two times because there are two mature sequences (not 3p and 5p but on the same strand))
for(i in 1:length(output)){
  if(num_dup[[names(output)[i]]] != 1){
    print(paste(num_dup[[names(output)[i]]], names(output)[[i]]))
  }
}

#now the names(positions_in_hairpin) are e.g. "MIR161 *star*" and "MIR161 *star*.1", but I need miR161.1 and miR161.2 like in the ath.gff3
#therefore I cut the *star*
names(positions_in_hairpin) <- sub(" \\*star\\*", "", names(positions_in_hairpin))
names(sequence) <- sub(" \\*star\\*", "", names(sequence))
#for the ath-MIR161, ath-MIR447a, ath-MIR779, ath-MIR780 and ath-MIR869 are two mature sequences (not 3p and 5p but on the same strand)
#I have to rename the names(positions_in_hairpin) of this MIRs to have .1 and .2 (or no .x) the same way as in the ath.gff3
#therefore I check manually which sequence corresponds to which miRNA (.1 or .2)
sequence$MIR161
sequence$MIR161.1
names(positions_in_hairpin)[names(positions_in_hairpin) == "MIR161"] <- "MIR161.2"
names(sequence)[names(sequence) == "MIR161"] <- "MIR161.2"

sequence$MIR447a
sequence$MIR447a.1
names(positions_in_hairpin)[names(positions_in_hairpin) == "MIR447a.1"] <- "MIR447a.2"
names(positions_in_hairpin)[names(positions_in_hairpin) == "MIR447a"] <- "MIR447a-3p"
names(sequence)[names(sequence) == "MIR447a.1"] <- "MIR447a.2"
names(sequence)[names(sequence) == "MIR447a"] <- "MIR447a-3p"

sequence$MIR779
sequence$MIR779.1
names(positions_in_hairpin)[names(positions_in_hairpin) == "MIR779"] <- "MIR779.2"
names(sequence)[names(sequence) == "MIR779"] <- "MIR779.2"

sequence$MIR780
sequence$MIR780.1
names(positions_in_hairpin)[names(positions_in_hairpin) == "MIR780.1"] <- "MIR780.2"
names(positions_in_hairpin)[names(positions_in_hairpin) == "MIR780"] <- "MIR780.1"
names(sequence)[names(sequence) == "MIR780.1"] <- "MIR780.2"
names(sequence)[names(sequence) == "MIR780"] <- "MIR780.1"

sequence$MIR869
sequence$MIR869.1
names(positions_in_hairpin)[names(positions_in_hairpin) == "MIR869.1"] <- "MIR869.2"
names(positions_in_hairpin)[names(positions_in_hairpin) == "MIR869"] <- "MIR869.1"
names(sequence)[names(sequence) == "MIR869.1"] <- "MIR869.2"
names(sequence)[names(sequence) == "MIR869"] <- "MIR869.1"

#create files for VIVA

#therefore I need to read the ath.gff3 file from miRBase
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data")
#Read miRBase gff3 file
ath <- read.table("ath.gff3", sep="\t", quote="")
#filter out the miRNAs with 3p and 5p for which we already did the analysis
#names of the miRNAs for which we already did the analysis
names <- read.table("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/VIVA/names.csv", sep = ",")
names <- paste0(names[ ,1], collapse = "|")
ath <- ath[!grepl(pattern = names, x = ath$V9, ignore.case = TRUE), ]
#create list of names for the new miRNAs for which we did not know the 3p and 5p miRNA
names_new <- sub(".*Name=([^;]+);.*", "\\1", ath[ath$V3 == "miRNA", ]$V9) #only use the mature sequences to extract the names
names_new <- sub("-[35]p$", "", names_new)
names_new <- sub("ath-", "", names_new)
#delete the miR829-3p.1, miR829-3p.2, miR1886.2 and miR1886.3 from the names because they are 3p and 5p and belong to the group of the 93 for which we know 3p and 5p
names_new <- names_new[names_new != "miR829-3p.1" & names_new != "miR829-3p.2" & names_new != "miR1886.2" & names_new != "miR1886.3"]
#delete "miR156i" because no prediction was made from StarSeeker
names_new <- names_new[names_new != "miR156i"]
#I have 237 names for the miRNAs for which we don't know 3p and 5p (names_new)
#I have 93 miRNAs for which we know 3p and 5p
#237 + 93 = 330
#There are only 326 MIRs in miRBase --> 330 - 5 = 325 (I can subtract 5 because for the hairpins ath-MIR161, ath-MIR447a, ath-MIR779, ath-MIR780 and ath-MIR869 there are two entries in names_new because there are two mature sequences (not 3p and 5p but on the same strand))
#Attention: In the end I will have 325 + 5 + 2 = 332 heatmaps in total
length(names_new) #237

#change name from "miR447a" in names_new to "miR447a-3p" to identify it correctly in the grepl() statement in the for loop
names_new[names_new == "miR447a"] <- "miR447a-3p"
#change name from "miR829" in names_new to "miR829-3p.1"
names_new[names_new == "miR829"] <- "miR829-3p.1"
#also change it in the positions_in_hairpin
names(positions_in_hairpin)[names(positions_in_hairpin) == "MIR829"] <- "MIR829-3p.1"
names(sequence)[names(sequence) == "MIR829"] <- "MIR829-3p.1"
#change "miR1886" to "miR1886.1" because it appears like this in names_new
names(positions_in_hairpin)[names(positions_in_hairpin) == "MIR1886"] <- "MIR1886.1"
names(sequence)[names(sequence) == "MIR1886"] <- "MIR1886.1"
#change the MI to mi
names(positions_in_hairpin) <- sub("MI", "mi", names(positions_in_hairpin))
names(sequence) <- sub("MI", "mi", names(sequence))

#save names of the hairpins without 3p and 5p miRNA as csv file for VIVA (array job)
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/VIVA")
write.table(names_new, file = "predictions_names.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

#create csv file for VIVA (with chr,position)
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/miRNA_data")
#Read miRBase gff3 file
ath <- read.table("ath.gff3", sep="\t", quote="")

#wd for csv files
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/VIVA/csv_all_predicted")
#I also want to create a csv file with the information about the start position of base_one,miRNA_one,loop,miRNA_two and base_two
#miRNA_one is always the 5p
df_starting_pos <- data.frame()
#iterate over all names and create csv files
for(name in names_new){
  #row in ath.gff3 for name
  miR <- ath[grepl(pattern = name, x = ath$V9), ]
  #making sure every name was found exactly one time in the ath.gff3
  if(length(miR$V1) != 1){
    print(paste("Error1!!!:",name))
  }
  #position of the predicted miRNA in the hairpin
  pos_hair <- positions_in_hairpin[[name]]
  #making sure every name was found in positions_in_hairpin
  if(is.null(pos_hair)){
    print(paste("Error2!!!:",name))
  }
  #find the line for the hairpin in the ath.gff3
  pattern <- sub("mi","MI",sub("\\.1|\\.2|-3p|-3p\\.1","",name))
  hairpin_start <- ath[grepl(pattern = pattern, x = ath$V9), ]$V4
  hairpin_end <- ath[grepl(pattern = pattern, x = ath$V9), ]$V5
  #making sure one exactly one line was found for the hairpin
  if(length(hairpin_start) != 1){
    print(paste("Error3!!!:",pattern, name))
  }
  
  #on the plus strand the 5p is the one with the smaller starting base
  if(miR$V7 == "+"){
    #start position of the predicted miRNA in the genome
    pos_genome <- hairpin_start + pos_hair - 1
    
    miRNA_5p <- min(miR$V4, pos_genome)
    miRNA_3p <- max(miR$V4, pos_genome)
    
    #for the csv file with all starting positions
    if(miRNA_5p == miR$V4){
      new_starting_pos <- data.frame(name = name, base_one = miRNA_5p-25, miRNA_one = miRNA_5p, loop = miR$V5+1, miRNA_two = miRNA_3p, base_two = miRNA_3p + nchar(sequence[[name]]))
      
      #input for VIVA
      df <- data.frame()
      #-25 and +25 for the base
      for(i in (miRNA_5p-25):(miRNA_3p+nchar(sequence[[name]])+24)){
        new <- data.frame(chr = substr(miR$V1,4,4), position = i)
        df <- rbind(df, new)
      }
    }
    if(miRNA_5p == pos_genome){
      new_starting_pos <- data.frame(name = name, base_one = miRNA_5p-25, miRNA_one = miRNA_5p, loop = miRNA_5p + nchar(sequence[[name]]), miRNA_two = miRNA_3p, base_two = miR$V5+1)
      
      #input for VIVA
      df <- data.frame()
      #-25 and +25 for the base
      for(i in (miRNA_5p-25):(miR$V5+25)){
        new <- data.frame(chr = substr(miR$V1,4,4), position = i)
        df <- rbind(df, new)
      }
    }
    df_starting_pos <- rbind(df_starting_pos, new_starting_pos)
  }
  #on the minus strand the 5p is the one with the bigger starting base
  else{
    #start position of the predicted miRNA in the genome
    pos_genome <- hairpin_end - pos_hair + 1 - nchar(sequence[[name]])
    
    miRNA_5p <- max(miR$V4, pos_genome)
    miRNA_3p <- min(miR$V4, pos_genome)
    
    #for the csv file with all starting positions
    if(miRNA_5p == miR$V4){
      new_starting_pos <- data.frame(name = name, base_one = miR$V5+25, miRNA_one = miR$V5, loop = miR$V4-1, miRNA_two = miRNA_3p+nchar(sequence[[name]])-1, base_two = miRNA_3p - 1)
      
      #input for VIVA
      df <- data.frame()
      #-25 and +25 for the base
      for(i in (miRNA_3p-25):(miR$V5+25)){
        new <- data.frame(chr = substr(miR$V1,4,4), position = i)
        df <- rbind(df, new)
      }
    }
    if(miRNA_5p == pos_genome){
      new_starting_pos <- data.frame(name = name, base_one = miRNA_5p + nchar(sequence[[name]])-1+25, miRNA_one = miRNA_5p + nchar(sequence[[name]])-1, loop = miRNA_5p-1, miRNA_two = miR$V5, base_two = miRNA_3p - 1)
      
      #input for VIVA
      df <- data.frame()
      #-25 and +25 for the base
      for(i in (miRNA_3p-25):(miRNA_5p+nchar(sequence[[name]])+24)){
        new <- data.frame(chr = substr(miR$V1,4,4), position = i)
        df <- rbind(df, new)
      }
    }
    df_starting_pos <- rbind(df_starting_pos, new_starting_pos)
  }
  write.table(df, file = paste(name, ".csv", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
}
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Create heatmaps")
write.table(df_starting_pos, file = "regions_predicted_miRNAs.csv", quote = FALSE, row.names = FALSE, sep = ",")

