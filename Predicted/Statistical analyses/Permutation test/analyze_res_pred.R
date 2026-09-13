#analyze the results
results_all <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Statistical tests/results_all_pred.csv", sep = ",", header = TRUE)
results_filtered <- results_all[results_all$pvalue < 0.05, ]
write.table(results_filtered, file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Predict miRNAs/Statistical tests/results_pvalue_0.05_pred.csv", quote = FALSE, sep = ",", row.names=FALSE)