setwd("./Google Drive/")
setwd("./vcfs/") #Set directory to vcf folder
library(VariantAnnotation)
VCF = readVcf("192.samples.snps.only.raw.vcf", "hg19") #read VCF file into R

VCF_GT <- geno(VCF)$GT
write.table(VCF_GT, file="VCF_GT.txt", sep =" ", row.names = TRUE, col.names = TRUE)

VCF_GT <- gsub("0/0", "0", VCF_GT)
VCF_GT <- gsub("1/1", "1", VCF_GT) #Recode
VCF_GT <- t(VCF_GT) #Transpose matrix


VCF_GT <- as.data.frame(VCF_GT)

write.csv(VCF_GT, file = "Mydata.csv")

setwd("/Users/Jialin/Google Drive/phenotypes/")

Phenotype <- read.delim("phenotypes.txt", header = T, sep ='\t')

total <- cbind(VCF_GT, Phenotype)
str(total)

total$SHLP2 <- as.numeric(total$SHLP2)

X2 <- function(x) {aggregate(SHLP2 ~ x, total, mean)} 

Average_SHLP2 <- apply(total, 2, X2)

X3 <- function(x){
if(length(total$SHLP2[total[[x]] == "0"]) > 1 && length(total$SHLP2[total[[x]] == "1"]) >1){
t.test(total$SHLP2[total[[x]] == "0"], total$SHLP2[total[[x]] == "1"])} else "Only one observation"
}


X4 <- function(x){
if(length(total$SHLP2[total[[x]] == "0"]) > 1 && length(total$SHLP2[total[[x]] == "1"]) >1){
tests[[x]] <- t.test(total$SHLP2[total[[x]] == "0"], total$SHLP2[total[[x]] == "1"])$p.value} else "Only one observation"
}
lapply(1:1054, X4)

SHLP2_pvalue <- lapply(1:1054, X4)
write.csv(SHLP2_pvalue, file = "SHLP2_pvalue.csv")

means <- list()
X5 <- function(x){
means[[x]] <- mean(total$SHLP2[total[[x]] == "0"])}

X6 <- function(x){
means[[x]] <- mean(total$SHLP2[total[[x]] == "1"])}

means_zero <- lapply(1:1054, X5)
means_one <- lapply(1:1054, X6)
write.csv(means_one, file="1means.csv")
write.csv(means_zero, file="0means.csv")
Merged <- rbind(means_zero, means_one, SHLP2_pvalue)
Merged <- t(Merged)
write.csv(Merged, file = "Final_data.csv")


