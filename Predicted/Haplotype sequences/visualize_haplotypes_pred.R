#TASK: For every MIR we have the haplotypes. Now I need to extract their frequency and visualize the results
library(Biostrings)
library(msa)

source("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Create heatmaps/functions.R")
setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Haplotype sequences")

#read information about name and strand of all MIRs
MIRs <- read.table("coordinates_pred.bed")
#read information about positions of the regions
positions <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Create heatmaps/regions_predicted_miRNAs.csv", header = TRUE, sep = ",")
#iterate over all MIRs
for(name in MIRs$V4){
  setwd("C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Haplotype sequences")
  #read reverse complement when on the minus strand
  if(MIRs[MIRs$V4 == name, 5] == "-"){
    fa <- readDNAStringSet(paste0("seqkit_rc/", name, "_rc.fa"))
  }
  else{
    fa <- readDNAStringSet(paste0("bcftools_consensus/", name, ".fa"))
  }
  
  #extract the frequencies of the sequences
  seq_freq <- sort(table(fa), decreasing = TRUE)
  char_vector <- setNames(names(seq_freq), paste0(seq_freq, "x"))
  #save unique sequences with frequencies in fa file
  fa_ <- DNAStringSet(char_vector)
  dir.create(paste0("Results/", name))
  writeXStringSet(fa_, file = paste0("Results/", name, "/frequencies.fa"))
  
  #create multiple sequence alignment when there is more than one sequence
  if(length(seq_freq) > 1){
    msa <- as.matrix(msa(char_vector, type = "dna", order = "input"))
  }
  else {
    msa <- matrix(strsplit(char_vector[[1]], "")[[1]], nrow = 1)
    rownames(msa) <- names(char_vector)
  }
  #what will be called "reference sequence" is the sequence with the highest frequency
  #extract positions of the reference sequence
  pos <- c()
  count <- 1
  for(i in 1:ncol(msa)){
    #indel
    if(msa[1, i] == "-"){
      pos <- c(pos, paste0(count, "_indel"))
    }
    else{
      pos <- c(pos, count)
      count <- count + 1
    }
  }
  colnames(msa) <- pos
  
  #extract positions with variants
  variant_pos <- apply(msa, MARGIN = 2, FUN = function(col){length(unique(col)) > 1})
  
  #create matrix that has a "." when the base is identical with the reference base
  #change values in msa to numerical for the heatmaps
  if(length(msa[,1]) > 1){
    mat <- matrix(data = NA, nrow = nrow(msa), ncol = ncol(msa))
    mat[1, ] <- msa[1, ]
    colnames(mat) <- pos
    rownames(mat) <- names(fa_)
    for(i in 1:ncol(msa)){
      for(j in 2:nrow(msa)){
        #if the base is the same as in the reference sequence
        if(msa[j,i] == msa[1,i]){
          msa[j,i] <- "1"
          mat[j,i] <- "."
        }
        else{
          mat[j,i] <- msa[j,i]
          msa[j,i] <- "3"
        }
        
      }
    }
    #save whole msa 
    write.table(mat, file = paste0("Results/", name, "/msa_whole.tsv"), sep = "\t", col.names = NA, quote = FALSE)
    #save only positions with variants
    write.table(mat[ ,variant_pos, drop = FALSE], file = paste0("Results/", name, "/msa_variants.tsv"), sep = "\t", col.names = NA, quote = FALSE)
    
    rownames(msa) <- make.unique(rownames(mat))
  }
  else{
    #just save the one sequence
    write.table(msa, file = paste0("Results/", name, "/msa_whole.tsv"), sep = "\t", col.names = NA, quote = FALSE)
  }
  
  msa[1, ] <- 1
  
  #start positions of the different regions of the miRNA (when the miRNA starts at 1 and ends at length(miRNA))
  i <- positions$name == name
  base_one <- 1
  miRNA_one <- abs(positions[i, "miRNA_one"] - positions[i, "base_one"]) + 1
  loop <- abs(positions[i, "loop"] - positions[i, "base_one"]) + 1
  miRNA_two <- abs(positions[i, "miRNA_two"] - positions[i, "base_one"]) + 1
  base_two <- abs(positions[i, "base_two"] - positions[i, "base_one"]) + 1
  
  create_heatmap_haplo(msa, base_one, miRNA_one, loop, miRNA_two, base_two, 0, 3, name, FALSE, colorscale = list(list(0.0, "white"), list(1/3, "darkgrey"), list(2/3, "green"), list(1.0, "red")), path = paste0("Results/", name))                  
}
