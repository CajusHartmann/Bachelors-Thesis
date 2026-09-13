library(magick)
library(pdftools)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

coords <- read.table("/data/laub-rna/hiwi/cajus/thesis/Haplotpye_Sequences/predicted_repeated/bcftools_consensus_pred_repeated/coordinates_pred.bed", sep = "\t")
name <- coords[task_id, 4]

setwd(paste0("/work/hartmaca/RNAfold_predicted_repeated/", name))

#list all pdf files in current dir
files <- list.files(pattern = ".pdf",full.names = TRUE)
#convert pdf to png
pngs <- unlist(lapply(files, function(f) {pdf_convert(f, format = "png", dpi = 300)}))
#delete pngs at the end of the script
on.exit(file.remove(pngs))
#read pngs
imgs <- image_read(pngs)
#cut off white background around the secondary structure
imgs_trimmed <- image_trim(imgs, fuzz = 20)
#add a small white background, to add distance between the structures
imgs_padded <- image_border(imgs_trimmed, color = "white", geometry = "80x80")
#extract frequencies of haplotypes
freqs <- gsub("_rss\\.pdf", "", basename(files))
#add frequencies to the structures
imgs_annotated <- image_annotate(imgs_padded,text = freqs, size = 30, color = "black", gravity = "north")
#merge images next to each other
result <- image_append(imgs_annotated)
#save
image_write(result, paste0("/work/hartmaca/RNAfold_predicted_repeated/merged_pdfs/", name, ".png"))
