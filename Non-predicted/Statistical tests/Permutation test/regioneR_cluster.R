library(GenomicRanges)
library(regioneR)

#read the file with the positions of the SNPs (I created this file with bcftools query)
positions <- read.table(file = "/work/hartmaca/bcftools_query/chrom_pos_all.vcf", sep = " ", col.names = c("CHROM", "POS"))
#create GRanges
gr <- GRanges(seqnames = positions$CHROM, ranges = IRanges(positions$POS, width = 1))

#get job id
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#names of the 95 non predicted miRNAs
names <- read.table(file = "/work/hartmaca/regioneR/regions_miRNAs.csv", header = TRUE, sep = ",")[,1]
#name for the current job
name <- names[task_id]
pos <- read.table(file = paste0("/data/laub-rna/hiwi/cajus/thesis/data/csv_all/", name, ".csv"), header = FALSE, sep = ",")
el <- GRanges(seqnames = pos[1,1], ranges = IRanges(start = pos[1,2]+25, end = pos[length(pos[,1]),2]-25), names = name)

#define TAIR10 genome (https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana/dna/)
genome <- GRanges(seqnames = 1:5, ranges = IRanges(start = 1, end = c(30427671, 19698289, 23459830, 18585056, 26975502)))
#This test asks if the miRNA region (el) overlaps SNPs (gr) more often than expected by chance
#SNPs are fixed and miRNA regions are randomly redistributed on the genome
#real overlaps of SNPs and miRNA regions are then compared with the overlaps from the permutations to calculate pvalues and zscores
res <- overlapPermTest(A = el, B = gr, genome = genome, alternative = "auto", ntimes = 1000, per.chromosome = TRUE, verbose = TRUE)
summ <- summary(res)
#write in file
new <- data.frame(name = name, pvalue = summ$pvalue, zscore = summ$zscore, test = summ$test)
write.table(new, file = paste0("/work/hartmaca/regioneR/res_", name, ".csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
