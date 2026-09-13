#examine if there is a statistical significant difference between the processing types (looking at the number of SNPs in the different regions)
tab <- read.table(file = "C:/Users/cajus/OneDrive/Dokumente/Bildung/Studium/6. Semester/Bachelorarbeit/Code/Other plots per group/statistical_analysis_non_predicted.csv", header = TRUE, sep = ",")
#test if the data is normally distributed
by(tab$num_SNPs, tab$region, shapiro.test) #all pvalues are smaller than 0.01 --> data is not normally distributed
#because the data is not normally distributed and not paired I apply the Kruskal Wallis test (for every region separately)
res_b1 <- kruskal.test(num_SNPs~type, data = subset(tab, region == "base_one"))
res_m1 <- kruskal.test(num_SNPs~type, data = subset(tab, region == "miRNA_one"))
res_loop <- kruskal.test((num_SNPs/length)~type, data = subset(tab, region == "loop")) #I have to divide by length because the loops can have different lengths
res_m2 <- kruskal.test(num_SNPs~type, data = subset(tab, region == "miRNA_two"))
res_b2 <- kruskal.test(num_SNPs~type, data = subset(tab, region == "base_two"))
#I have to adjust the pvalues because I have 5 tests
p.adjust(c(0.5891, 0.3984, 0.3153, 0.4437, 0.6395), method = "bonferroni")
#result: no statistical significant difference between the processing types

#Now I want so see if it would make a difference to get rid of the outliers (SNPs that occur in more than 700 accession)
tab <- tab[tab$num_SNPs < 700, ]
#Now I just repeat everything
#test if the data is normally distributed
by(tab$num_SNPs, tab$region, shapiro.test) #all pvalues are smaller than 0.01 --> data is not normally distributed
#because the data is not normally distributed and not paired I apply the Kruskal Wallis test (for every region separately)
res_b1 <- kruskal.test(num_SNPs~type, data = subset(tab, region == "base_one"))
res_m1 <- kruskal.test(num_SNPs~type, data = subset(tab, region == "miRNA_one"))
res_loop <- kruskal.test((num_SNPs/length)~type, data = subset(tab, region == "loop")) #I have to divide by length because the loops can have different lengths
res_m2 <- kruskal.test(num_SNPs~type, data = subset(tab, region == "miRNA_two"))
res_b2 <- kruskal.test(num_SNPs~type, data = subset(tab, region == "base_two"))
#I have to adjust the pvalues because I have 5 tests
p.adjust(c(0.2559, 0.4971, 0.5216, 0.7165, 0.8113), method = "bonferroni")
#result: still no statistical significant difference between the processing types

