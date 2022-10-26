library(reshape2)

args = commandArgs(trailingOnly=TRUE)
options(scipen=999)

pc1 = read.table(args[1], skip=1, sep="\t")
rownames(pc1) = paste(pc1$V2, pc1$V3-1, sep="-")
pc1$V7 = as.numeric(pc1$V7)


hicMat = matrix(0, ncol=3, nrow=0)
for(f in list.files(path=args[2], pattern=args[3], full.name=T)){
	ia = read.table(f, skip=1, sep="\t", row.names=1)
	ia = ia[,-1]
	colnames(ia) = rownames(ia)
	ia[upper.tri(ia)] <- NA
	diag(ia) <- NA
	ia = melt(as.matrix(ia))
	ia = na.omit(ia)
	ia[,3] = as.numeric(ia[,3])
	#print(head(ia))
	hicMat = rbind(hicMat, ia)
}
#hicMat = hicMat[hicMat$V1 != hicMat$V2, ]


hicMat = data.frame(hicMat, stringsAsFactors=F)
hicMat$pos1 = as.numeric(unlist(strsplit(as.character(hicMat$Var1), "-"))[seq(2,2*nrow(hicMat),2)])
hicMat$pos2 = as.numeric(unlist(strsplit(as.character(hicMat$Var2), "-"))[seq(2,2*nrow(hicMat),2)])
hicMat$distance = abs(hicMat$pos1-hicMat$pos2)

library(dplyr)
freqDis = hicMat %>%
	group_by(distance) %>%
	summarise(freq = sum(value))

total_IA = sum(freqDis$freq)
freqDis$freq = freqDis$freq / total_IA
freqDis$sample = args[3]
write.table(freqDis, paste0(args[3],".FreqDis.txt"), quote=F, sep="\t")

# A compartment
freqDis_A = hicMat %>%
	filter(Var1 %in% rownames(pc1[pc1$V7>0,])) %>%
	filter(Var2 %in% rownames(pc1[pc1$V7>0,])) %>%
        group_by(distance) %>%
        summarise(freq = sum(value))
freqDis_A$freq = freqDis_A$freq / total_IA
freqDis_A$sample = args[3]
write.table(freqDis_A, paste0(args[3],".A.FreqDis.txt"), quote=F, sep="\t")

# B compartment
freqDis_B = hicMat %>%
        filter(Var1 %in% rownames(pc1[pc1$V7<0,])) %>%
        filter(Var2 %in% rownames(pc1[pc1$V7<0,])) %>%
        group_by(distance) %>%
        summarise(freq = sum(value))
freqDis_B$freq = freqDis_B$freq / total_IA
freqDis_B$sample = args[3]
write.table(freqDis_B, paste0(args[3],".B.FreqDis.txt"), quote=F, sep="\t")



