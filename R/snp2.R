#analysis of qPCR results for Dilution and Preservation study

library(ggplot2)
library(reshape2)
library(scales)
library(stringr)

#set working directory for PT95 folder
setwd("~/Documents/LT2_genome_analysis/vcfUG/")

#Import of datafiles
vcf <- list.files()[grep("+.vcf", list.files(), value = F)]
snp <- list.files()[grep("+.snp", list.files(), value = F)]

autoTable = data.frame()
for(i in 1:length(vcf)){
  file <- readLines(vcf[i])
  lines <- file[!(1:length(file) %in% grep("#", file))]
  fileTable <- data.frame(colsplit(lines,"\t", 
                                   names = c("CHROM",  "POS",  "ID",	"REF",	"ALT",	"QUAL",	
                                             "FILTER",	"INFO",	"FORMAT",	"LT2")))
  fileTable$data <- vcf[i]
  autoTable <- rbind(autoTable, fileTable)
}

autoTable <- cbind(autoTable, colsplit(autoTable$data, pattern = "_", names = c("mapper", "reads", "ref")))

autoTable$reads <- factor(autoTable$reads, levels = c("Ref", "mut", "mut5k"))
autoTable$ref <- factor(autoTable$ref, levels = c("raw.vcf", "prinseq.vcf", "coral.vcf"))
ggplot(autoTable[autoTable$CHROM != "gi|17233403|ref|NC_003277.1|" & autoTable$mapper != "its",]) + 
  geom_bar(aes(x = QUAL)) + facet_grid(reads~ref, scale = "free")
+ scale_y_log10()
ggplot(autoTable[autoTable$CHROM != "gi|17233403|ref|NC_003277.1|",]) + 
  geom_bar(aes(x = POS, fill = mapper), position = "dodge") + facet_grid(reads~ref, scale = "free") 
+ scale_y_log10()

mutTable = data.frame()
for(i in 1:length(snp)){
  file <- readLines(snp[i])
  #lines <- file[!(1:length(file) %in% grep("#", file))]
  fileTable <- read.table(snp[i], col.names = c("CHROM","POS","REF","MUT","NUM"))
  fileTable$data <- snp[i]
  mutTable <- rbind(mutTable, fileTable)
}

mutTable$CP <- paste(mutTable$CHROM,mutTable$POS, sep = "_")
nonMut <- autoTable[!(autoTable$CP %in% mutTable$CP),]
sink("snp.txt")
for( i in unique(nonMut$data)){
  print(i)
  print(sort(nonMut$POS[nonMut$data == i]))
}
sink()


#venn diagram
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)
autoTable$CP <- paste(autoTable$CHROM,autoTable$POS, sep = "_")
venn <- unique(rbind(autoTable[, 1:2],mutTable[, 1:2]))
venn$CP <- paste(venn$CHROM,venn$POS, sep = "_")

for(i in vcf){
  name <- as.character(strsplit(i, '\\.')[[1]][1])
  venn[[name]] <- (venn$CP %in% autoTable$CP[autoTable$data == i])
}

#new approach
length(unique(venn$CP))
test <- vennCounts(venn[,4:15])
summary(test)
testdf <- as.data.frame(test[,1:13])
testdfSort <- arrange(testdf[testdf$Counts > 5,], Counts)
testdfSort$ID <- rownames(testdfSort)
testdfSortM <- melt(testdfSort, id = c("Counts", "ID"))
testdfSortM$CountF <- as.factor(testdfSortM$Counts)
testdfSortM$value <- as.factor(testdfSortM$value)
ggplot(testdfSortM) + geom_raster(aes(x = variable, y = ID, fill = value))
ggplot(testdfSortM) + geom_bar(stat = "identity", aes(x = ID, y = Counts))
ggplot(testdfSortM[grep("Ref", testdfSortM$variable),]) + geom_raster(aes(x = variable, y = ID, fill = value))
ggplot(testdfSortM) + geom_bar(aes(x = CountF))



pdf("venneulerPlots.pdf")
RefLogical <- as.matrix(venn[,grep("Ref",colnames(venn))])
plot(venneuler(RefLogical))
mutALL <- grep("mut",colnames(venn))
mut <- grep("mut",colnames(venn))[!(grep("mut",colnames(venn)) %in% grep("mut5k",colnames(venn)))]
mutLogical <- as.matrix(venn[!(venn$CP %in% mutTable$CP[mutTable$data == "LT2_Ref.mut.snp"]),
                             grep("mut",colnames(venn))[!(grep("mut",colnames(venn)) %in% grep("mut5k",colnames(venn)))]])
plot(venneuler(mutLogical))
mut5kLogical <- as.matrix(venn[!(venn$CP %in% mutTable$CP[mutTable$data == "LT2_Ref.mut5k.snp"]),grep("mut5k",colnames(venn))])
plot(venneuler(mut5kLogical))
dev.off()

mut5ksnp <- venn[venn$CP %in% mutTable$CP[mutTable$data == "LT2_Ref.mut5k.snp"],]
mutsnp <- venn[venn$CP %in% mutTable$CP[mutTable$data == "LT2_Ref.mut.snp"],]

pdf("in_silico_muts_TP_FN.pdf")
mut5ksnp <- venn[venn$CP %in% mutTable$CP[mutTable$data == "LT2_Ref.mut5k.snp"],]
mutsnp <- venn[venn$CP %in% mutTable$CP[mutTable$data == "LT2_Ref.mut.snp"],]
mut5ksnpNoITS <- mut5ksnp[,grep("TMAP", colnames(mut5ksnp))]
RefLogical <- as.matrix(mut5ksnpNoITS[,grep("Ref",colnames(mut5ksnpNoITS))])
vennDiagram(vennCounts(RefLogical))
mutsnpNoITS <- mutsnp[,grep("TMAP", colnames(mutsnp))]
mutLogical <- as.matrix(mutsnpNoITS[,grep("mut",colnames(mutsnpNoITS))[!(grep("mut",colnames(mutsnpNoITS)) %in% grep("mut5k",colnames(mutsnpNoITS)))]])
vennDiagram(vennCounts(mutLogical))
mut5kLogical <- as.matrix(mut5ksnpNoITS[,grep("mut5k",colnames(mut5ksnpNoITS))])
vennDiagram(vennCounts(mut5kLogical))
mutsnpITS <- mutsnp[,grep("its", colnames(mutsnp))]
mutLogical <- as.matrix(mutsnp[,grep("its", colnames(mutsnp))])
vennDiagram(vennCounts(mutLogical))
mut5ksnpITS <- mut5ksnp[,grep("its", colnames(mut5ksnp))]
mut5kLogical <- as.matrix(mut5ksnp[,grep("its", colnames(mut5ksnp))])
vennDiagram(vennCounts(mut5kLogical))
dev.off()






ggplot(autoTable[grep("coral",autoTable$data), ]) + 
  geom_point(aes(x = POS, y = data, color = data))

coral_snp <- sort(autoTable$POS[autoTable$data == "LT2_coral.fastq.bwasw.vcf"])
coral_mut <- sort(autoTable$POS[autoTable$data == "mut.LT2_coral.fastq.vcf"])
snps <- sort(snp$V2)


##dataset for heat map analysis
RefSet <- autoTable[autoTable$reads == "Ref",]
non_coral <- unique(RefSet$POS[RefSet$ref != "coral.vcf" & RefSet$CHROM == "gi|16763390|ref|NC_003197.1|"])
RefSet$POS <- as.factor(RefSet$POS)
ggplot(RefSet[RefSet$CHROM == "gi|16763390|ref|NC_003197.1|" & RefSet$POS %in% non_coral,]) + 
  geom_raster(aes(x = ref, y = POS, fill = QUAL)) + facet_wrap(~mapper, scale = "free_x", as.table = T)

ggplot(autoTable) + geom_bar(aes(x = reads, fill = ref, color = mapper), position = "dodge")


#analysis for true positive and negatives for mutation snp calls in relation to snp reference calls 
#for the same read set
#raw TMAP analysis
T_raw_Ref <- autoTable[autoTable$data == "TMAP_Ref_raw.vcf",]
T_coral_Ref <- autoTable[autoTable$data == "TMAP_coral_raw.vcf",]
T_prinseq_Ref <- autoTable[autoTable$data == "TMAP_prinseq_raw.vcf",]

insert <-mutTable[grep(unique(mutTable$REF)[5],mutTable$REF),]
insert$INDEL <- "insert"

del <-mutTable[grep(unique(mutTable$MUT)[5],mutTable$MUT),]
del$INDEL <- "del"
insert_mut <- insert[insert$data == "LT2_Ref.mut.snp",]
insert_mut5k <- insert[insert$data == "LT2_Ref.mut5k.snp",]
del_mut <- del[del$data == "LT2_Ref.mut.snp",]
del_mut5k <- del[del$data == "LT2_Ref.mut5k.snp",]

indel_mut <- rbind(insert_mut, del_mut)
indel_mut <- arrange(indel_mut, desc(POS))

#still need to sort
indel_mut5k <- rbind(insert_mut5k, del_mut)



TrR <- T_raw_Ref
#mut genome raw TMAP
for(i in 1:nrow(indel_mut)){
  T_raw_Ref$POS[T_raw_Ref$POS >= indel_mut[i,2] & T_raw_Ref$CHROM == "gi|16763390|ref|NC_003197.1|" ] <- 
    T_raw_Ref$POS[T_raw_Ref$POS >= indel_mut[i,2] & T_raw_Ref$CHROM == "gi|16763390|ref|NC_003197.1|" ]- 1
}
T_raw_Ref_mut <- rbind(T_raw_Ref, autoTable[autoTable$data == "TMAP_mut_raw.vcf",])
CPunique <- data.frame(mut = unique(T_raw_Ref_mut$CP) %in% autoTable$CP[autoTable$data == "TMAP_mut_raw.vcf"], 
                       ref = unique(T_raw_Ref_mut$CP) %in% T_raw_Ref$CP)
vennDiagram(vennCounts(as.matrix(CPunique)))

$mut <- (venn$CP %in% autoTable$CP[autoTable$data == i])


