#####
# compartmentalization 
# saddlePlot
#


library(reshape2)


options(scipen=999)

pc1="../analysis/ABcompart/CLP1.PC1.txt"
mat_dir="../analysis/matrix/homer_500k_distNorm/"
sample="CLP1"

saddlePlot <- function(pc1, mat_dir, sample){

  pc1 = read.table(pc1, skip=1, sep="\t")
  #pc1 = read.table("../analysis/ABcompart/CLP1.PC1.txt",skip=1, sep="\t")
  rownames(pc1) = paste(pc1$V2, pc1$V3-1, sep="-")
  pc1$V7 = as.numeric(pc1$V7)
  pc1 = pc1[!is.na(pc1$V7),]
  pc1 = pc1[order(pc1$V7),]
  
  
  hicMat = matrix(0, ncol=3, nrow=0)
  for(f in list.files(path=mat_dir, pattern=sample, full.name=T)){
  	ia = read.table(f, skip=1, sep="\t", row.names=1)
    #ia = read.table("../analysis/matrix/homer_500k_distNorm/CLP1.chr1.500k.txt", skip=1, sep="\t", row.names=1)
  	ia = ia[,-1]
  	colnames(ia) = rownames(ia)
  	ia[upper.tri(ia)] <- NA
  	diag(ia) <- NA
  	ia = melt(as.matrix(ia))
  	ia = na.omit(ia)
  	ia[,3] = as.numeric(ia[,3])
  	hicMat = ia
  	#print(head(ia))
  	hicMat = rbind(hicMat, ia)
  }
  hicMat$Var1 = as.character(hicMat$Var1)
  hicMat$Var2 = as.character(hicMat$Var2)
  #
  #hicMat = hicMat[hicMat$value !=0, ]
  

  bin_order = rownames(pc1)[rownames(pc1) %in% hicMat$Var1]
  image(log(as.matrix(ia[bin_order,bin_order])))

  windows = as.numeric(seq(1,length(bin_order),length.out=31))
  res = matrix(0, ncol=30, nrow=30)
  for(i in 1:30){
    for(j in 1:30){
      regions1 = bin_order[windows[i]:windows[i+1]]
      regions2 = bin_order[windows[j]:windows[j+1]]
      res[i,j] = mean(hicMat[(hicMat$Var1 %in% regions1) & (hicMat$Var2 %in% regions2), 3])
    }
  }
  
  
  png(paste0("../figure/Compartmentalization/", sample, ".saddleplot.png"), 430, 400)
  #layout(matrix(c(1,4,2,3), 2, 2, byrow = TRUE), heights=c(1,5), widths=c(5,1))
  #par(mar=c(0.5,1,1,0))
  #barplot(pc1[bin_order, 7], col =colorRampPalette(c("darkblue", "white","darkred"))(length(bin_order)), border = NA,axes =F, space=0, xaxs='i')
  #par(mar=c(1,1,0,0))
  par(mar=c(1,1,1,1))
  image(res,col = colorRampPalette(c("white", "red", "black"))(20),frame.plot=TRUE,axes=F)
  #par(mar=c(1,0.5,0,1))
  #barplot(pc1[bin_order, 7],  col =colorRampPalette(c("darkblue", "white","darkred"))(length(bin_order)), horiz=TRUE,border = NA,axes=F, space=0, yaxs='i')
  dev.off()
  
  remove(hicMat)
  
}

for (i in c("LT1", "ST1","MPP2", "CLP1","CMP1","GMP1","MEP1","MKP1","MK1","G1")) {
 saddlePlot(paste0("../analysis/ABcompart/",i,".PC1.txt"), "../analysis/matrix/homer_500k_distNorm/", paste0(i,".chr1.500k"))
}



#
saddlePlot1 <- function(pc1, mat){
  
  pc1 = read.table(pc1, skip=1, sep="\t")
  pc1 = read.table("../homer/CLP1.PC1.bedGraph",skip=1)
  rownames(pc1) = paste(pc1$V1, pc1$V2, sep="-")
  pc1$V4 = as.numeric(pc1$V4)
  pc1 = pc1[!is.na(pc1$V4),]
  pc1 = pc1[order(pc1$V4),]
  
    ia = read.table(mat, skip=1, sep="\t", row.names=1)
    #ia = read.table("../analysis/matrix/homer_500k_distNorm/CLP1.chr1.500k.txt", skip=1, sep="\t", row.names=1)
    ia = ia[,-1]
    colnames(ia) = rownames(ia)
    qt = quantile(na.omit(as.vector(as.matrix(ia))), c(0.05, 0.95))
    ia[ia > qt[2]] <- qt[2]
    ia[ia < qt[1]] <- qt[1]
  
  bin_order = rownames(pc1)[rownames(pc1) %in% rownames(ia)]
  ia = ma3x3.matrix(as.matrix(ia[bin_order,bin_order]))

  sample = strsplit(basename(mat), "[.]")[[1]][1]
  png(paste0("../figure/Compartmentalization/", sample, "_chr1.saddleplot.png"), 430, 400)
  layout(matrix(c(1,4,2,3), 2, 2, byrow = TRUE), heights=c(1,6), widths=c(6,1))
  par(mar=c(0.5,1,1,0))
  barplot(pc1[bin_order, 4], col =colorRampPalette(c("darkblue", "white","darkred"))(length(bin_order)), border = NA,axes =F, space=0, xaxs='i')
  par(mar=c(1,1,0,0))
  image(ia,col = colorRampPalette(c("royalblue2","white", "salmon"))(20),frame.plot=TRUE,axes=F)
  par(mar=c(1,0.5,0,1))
  barplot(pc1[bin_order, 4],  col =colorRampPalette(c("darkblue", "white","darkred"))(length(bin_order)), horiz=TRUE,border = NA,axes=F, space=0, yaxs='i')
  dev.off()

  
}

saddlePlot1("../analysis/ABcompart/CLP1.PC1.txt", paste0("../analysis/matrix/homer_500k_distNorm/CLP1.chr1.500k.txt"))

for (i in c("LT1", "ST1","MPP2", "CLP1","CMP1","GMP1","MEP1","MKP1","MK1","G1")) {
  saddlePlot1(paste0("../analysis/ABcompart/",i,".PC1.txt"), paste0("../analysis/matrix/homer_500k_distNorm/", i,".chr1.500k.txt"))
}

#

compartmentStrength <- function(pc1, mat){
  pc1 = read.table(pc1, skip=1, sep="\t")
  #pc1 = read.table("../analysis/ABcompart/CLP1.PC1.txt",skip=1, sep="\t")
  rownames(pc1) = paste(pc1$V2, pc1$V3-1, sep="-")
  
    ia = read.table(mat, skip=1, sep="\t", row.names=1)
    #ia = read.table("../analysis/matrix/homer_500k_balance/CMP1.chr14.500k.txt", skip=1, sep="\t", row.names=1)
    ia = ia[,-1]
    colnames(ia) = rownames(ia)
    ia[upper.tri(ia)] <- NA
    diag(ia) <- NA
    ia = melt(as.matrix(ia))
    ia = na.omit(ia)
    ia[,3] = as.numeric(ia[,3])

  hicMat = ia
  hicMat$Var1 = as.character(hicMat$Var1)
  hicMat$Var2 = as.character(hicMat$Var2)
  #
  region_A = rownames(pc1[pc1$V9 == "A",])
  region_B = rownames(pc1[pc1$V9 == "B",])
  AA = mean(hicMat[hicMat$Var1 %in% region_A & hicMat$Var2 %in% region_A, 3])
  AB = mean(hicMat[(hicMat$Var1 %in% region_A & hicMat$Var2 %in% region_B) | (hicMat$Var1 %in% region_B & hicMat$Var2 %in% region_A), 3])
  BB = mean(hicMat[hicMat$Var1 %in% region_B & hicMat$Var2 %in% region_B, 3])
  return(c(AA, AB, BB))
}

compartmentStrength("../analysis/ABcompart/G1.PC1.txt", "../analysis/matrix/homer_500k_balance/G1.chr1.500k.txt")

compStrength = data.frame(row.names = c("sample", "chr", "AA", "AB", "BB"))
for (i in Sys.glob("../analysis/matrix/homer_500k_balance/*.500k.txt")) {
  bn = strsplit(basename(i), "[.]")[[1]]
  sample = bn[1]
  chr = bn[2]
  pc1 = paste0("../analysis/ABcompart/", sample,".PC1.txt")
  compStrength[[basename(i)]] = c(sample, chr, compartmentStrength(pc1, i))
}

compStrength = data.frame(t(compStrength), stringsAsFactors = F)
compStrength$chr = factor(compStrength$chr, levels = paste0("chr", c(seq(1:19), "X")))

compStrength$cell = stringr::str_sub(compStrength$sample, 1, -2)
compStrength$rep = stringr::str_sub(compStrength$sample, -1, -1)
compStrength$cell = factor(compStrength$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))
compStrength$AA = as.numeric(compStrength$AA)
compStrength$AB = as.numeric(compStrength$AB)
compStrength$BB = as.numeric(compStrength$BB)
compStrength$total = compStrength$AB+compStrength$AA+compStrength$BB
compStrength$CS = (compStrength$AA+compStrength$BB)/compStrength$total
g=ggplot(compStrength[compStrength$chr != "chrX",], aes(cell, CS)) +
  #geom_jitter(width=0.20,aes(col=chr, shape=rep)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(lwd=1, outlier.shape = NA, fill=NA)+ #ylim(0.65,0.8)+
  theme_classic(base_size = 20) +
  scale_colour_manual(values=rainbow(25)[1:20])+xlab("")+ylab("Compartment Strength")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("../figure/Compartmentalization/CompartmentStrength.box.pdf", device = "pdf", width = 5, height = 4)

ggplot(compStrength[compStrength$chr != "chrX",], aes(cell, log2(BB/AA))) +
  geom_jitter(width=0.20,aes(col=chr, shape=rep)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(lwd=1, outlier.shape = NA, fill=NA)+ 
  theme_classic(base_size = 20) +
  scale_colour_manual(values=rainbow(25)[1:20])+xlab("")+ylab("BB/AA")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

