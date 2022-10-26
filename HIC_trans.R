#####
# trans-interaction
#

library(data.table)
library(reshape2)
library(gplots)
library(dplyr)

sum_trans_heatmap <- function(mat){
  LT1_trans = data.frame(fread(mat, sep="\t", drop=1), row.names = 1, stringsAsFactors = F)
  colnames(LT1_trans) = rownames(LT1_trans)
  LT1_trans.m = melt(as.matrix(LT1_trans))
  LT1_trans.m$chr1 = sapply(strsplit(as.character(LT1_trans.m$Var1), "-"), "[", 1)
  LT1_trans.m$chr2 = sapply(strsplit(as.character(LT1_trans.m$Var2), "-"), "[", 1)

  LT1_trans.m1 = LT1_trans.m[,3:5] %>% group_by(chr1, chr2) %>%
    summarise(mean(value))
  
  LT1_trans.m1 = dcast(LT1_trans.m1, chr1 ~ chr2)
  rownames(LT1_trans.m1) = LT1_trans.m1[,1]
  LT1_trans.m1 = LT1_trans.m1[,-1]
  
  chrom_order = paste0("chr", c(seq(1:19)))
  LT1_trans.m1 = LT1_trans.m1[chrom_order, chrom_order]
  diag(LT1_trans.m1) = NA
  LT1_trans.m1[upper.tri(LT1_trans.m1)] <- NA
  #heatmap.2(as.matrix(LT1_trans.m1),Rowv=FALSE, Colv=FALSE, trace="none", key=FALSE,keysize = 0.4, col = colorRampPalette(c("navy", "white", "firebrick3"))(20))
  
  par(mar=c(1,1.5,1,1))
  LT1_trans.m1[LT1_trans.m1>0.002] = 0.002
  image(t(as.matrix(LT1_trans.m1)), col = colorRampPalette(c("blue", "white", "red"))(n = 20),frame.plot=F,axes=F)
  axis(3, at=seq(0,1, length=19), labels=colnames(LT1_trans.m1), lwd=0, las=2)
  axis(2, at=seq(0,1, length=19), labels=colnames(LT1_trans.m1), lwd=0, las=1)
}



pdf("../figure/transInteraction/trans_chrom_summary.pdf",15,6)
layout(matrix(1:10, 2, 5, byrow = TRUE))
for(i in c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G")){
  sum_trans_heatmap(paste0("../analysis/matrix/homer_merge_2M_balance/",i,".2M.txt"))
}
dev.off()



plot_trans_heatmap <- function(mat, chr){
  if(chr == "all"){
    chr = paste0("chr", c(seq(1:19), "X"))
  }
  
  LT1_trans = data.frame(fread(mat, sep="\t", drop=1), row.names = 1, stringsAsFactors = F)
  #LT1_trans = data.frame(fread("../analysis/matrix/homer_2M_balance/MK1.2M.txt", sep="\t", drop=1), row.names = 1, stringsAsFactors = F)
  colnames(LT1_trans) = rownames(LT1_trans)
  LT1_trans.m = melt(as.matrix(LT1_trans))
  LT1_trans.m$chr1 = sapply(strsplit(as.character(LT1_trans.m$Var1), "-"), "[", 1)
  LT1_trans.m$chr2 = sapply(strsplit(as.character(LT1_trans.m$Var2), "-"), "[", 1)
  
  if(length(chr) == 2){
    LT1_trans.m = LT1_trans.m[LT1_trans.m$chr1 == chr[1] & LT1_trans.m$chr2 == chr[2],]
  } else{
    LT1_trans.m = LT1_trans.m[LT1_trans.m$chr1 %in% chr & LT1_trans.m$chr2 %in% chr,]
  }
  LT1_trans.m[LT1_trans.m$chr1 == LT1_trans.m$chr2, 3] <- NA
  LT1_trans.m1 = dcast(LT1_trans.m, Var1 ~ Var2, value.var="value")
  rownames(LT1_trans.m1) = LT1_trans.m1$Var1
  LT1_trans.m1 = LT1_trans.m1[,-1]
  LT1_trans.m1 = as.matrix(LT1_trans.m1)
  q = quantile(na.omit(as.vector(LT1_trans.m1)), c(0.1, 0.9))
  LT1_trans.m1[LT1_trans.m1 > q[2]] <- q[2]
  LT1_trans.m1[LT1_trans.m1 < q[1]] <- q[1]
  par(mar=c(1,1,1,1))
  #image((LT1_trans.m1), col = colorRampPalette(c("purple", "black", "yellow"))(n = 20),frame.plot=TRUE,axes=F)
  image((LT1_trans.m1[,rev(1:ncol(LT1_trans.m1))]), col = colorRampPalette(c("blue", "white", "red"))(n = 20),frame.plot=TRUE,axes=F)
  
}

plot_trans_heatmap("../analysis/matrix/homer_2M_balance/LT1.2M.txt", c("chr1","chrX"))

plot_trans_heatmap("../analysis/matrix/homer_2M_balance/LT1.2M.txt", paste0("chr", c(seq(1:8))))


for(i in Sys.glob("../analysis/matrix/homer_2M_balance/*2M.txt")){
  sample = strsplit(basename(i), "[.]")[[1]][1]
  #pdf(paste0("../figure/transInteraction/", sample,".2M.pdf"),4.2,4)
  png(paste0("../figure/transInteraction/", sample,".2M.png"), 1100,1000, res=300)
  #png(paste0("../figure/transInteraction/", sample,".2M_chr1-8.png"), 650,600, res=300)
  #par(mar=c(1,1,1,1))
  #plot_trans_heatmap(i, paste0("chr", c(seq(1:8))))
  plot_trans_heatmap(i, "all")
  dev.off()

}

# 2M resulation
plot_trans_2M_heatmap_pc1 <- function(mat, chr, pc1){
  pc1 = read.table(pc1,skip=1)
  #pc1 = read.table("../homer/MK1.PC1.bedGraph",skip=1)
  
  LT1_trans = data.frame(fread(mat, sep="\t", drop=1), row.names = 1, stringsAsFactors = F)
  #LT1_trans = data.frame(fread("../analysis/matrix/homer_2M_balance/MK1.2M.txt", sep="\t", drop=1), row.names = 1, stringsAsFactors = F)
  colnames(LT1_trans) = rownames(LT1_trans)
  LT1_trans.m = melt(as.matrix(LT1_trans))
  LT1_trans.m$chr1 = sapply(strsplit(as.character(LT1_trans.m$Var1), "-"), "[", 1)
  LT1_trans.m$chr2 = sapply(strsplit(as.character(LT1_trans.m$Var2), "-"), "[", 1)
  

  LT1_trans.m = LT1_trans.m[LT1_trans.m$chr1 == chr[1] & LT1_trans.m$chr2 == chr[2],]

  LT1_trans.m[LT1_trans.m$chr1 == LT1_trans.m$chr2, 3] <- NA
  LT1_trans.m1 = dcast(LT1_trans.m, Var1 ~ Var2, value.var="value")
  rownames(LT1_trans.m1) = LT1_trans.m1$Var1
  LT1_trans.m1 = LT1_trans.m1[,-1]
  LT1_trans.m1 = as.matrix(LT1_trans.m1)
  q = quantile(na.omit(as.vector(LT1_trans.m1)), c(0.1, 0.9))
  LT1_trans.m1[LT1_trans.m1 > q[2]] <- q[2]
  LT1_trans.m1[LT1_trans.m1 < q[1]] <- q[1]

  layout(matrix(c(2,1,4,3), 2, 2, byrow = TRUE), heights=c(6,1), widths=c(1,6))
  par(mar=c(0,0,1,1))
  image((LT1_trans.m1[,rev(1:ncol(LT1_trans.m1))]),col = colorRampPalette(c("blue", "white", "red"))(20),frame.plot=TRUE,axes=F)
  par(mar=c(0,0.5,1,0.5))
  pc1_rev = rev(pc1[pc1$V1==chr[2], 4])
  barplot(pc1_rev,  col = ifelse(pc1_rev > 0, "darkred", "darkblue"), horiz=TRUE,border = NA,axes=F, space=0, yaxs='i')
  par(mar=c(0.5,0,0.5,1))
  barplot(pc1[pc1$V1==chr[1], 4], col = ifelse(pc1[pc1$V1==chr[1], 4] > 0, "darkred", "darkblue"), border = NA,axes =F, space=0, xaxs='i')
  
}
plot_trans_2M_heatmap_pc1("../analysis/matrix/homer_2M_balance/MPP2.2M.txt", c("chr1","chr2"), paste0("../homer/","MPP2",".PC1.bedGraph"))
barplot(rep(1,20),col=colorRampPalette(c("blue", "white", "red"))(20), axes=F, border = NA,space=0)

for(i in Sys.glob("../analysis/matrix/homer_2M_balance/*2M.txt")){
  sample = strsplit(basename(i), "[.]")[[1]][1]
  if(sample=="MPP1"){
    next
  }
  png(paste0("../figure/transInteraction/", sample,".2M_chr1chr2.png"), 640,600, res=300)
  par(mar=c(0,0,0,0))
  plot_trans_2M_heatmap_pc1(i, c("chr1","chr2"), paste0("../homer/",sample,".PC1.bedGraph"))
  dev.off()
}



# 500K resulation
plot_trans_500K_heatmap_pc1 <- function(mat, abs_bed, chr, pc1){
  pc1 = read.table(pc1,skip=1)
  
  abs_bed = read.table(abs_bed)
  chr1_abs = abs_bed[abs_bed$V1 == chr[1], 4]
  chr2_abs = abs_bed[abs_bed$V1 == chr[2], 4]

  mat1 = read.table(mat, stringsAsFactors = F)
  mat1 = mat1[mat1$V1 %in% chr1_abs & mat1$V2 %in% chr2_abs,]
  mat1 = dcast(mat1, V1 ~ V2, value.var="V3")
  rownames(mat1) = mat1$V1
  mat1 = mat1[,-1]

  mat1 = as.matrix(mat1)
  mat1[is.na(mat1)] <- 0

  # hic matrix normalized by quantile
  q = quantile(mat1, probs = c(0.15,0.85))
  #q = quantile(mat1, probs = c(0.2,0.8))
  mat1[mat1 > q[2]] <- q[2]
  mat1[mat1 < q[1]] <- q[1]
  #image((mat1),col = colorRampPalette(c("blue", "white", "red"))(20),frame.plot=TRUE,axes=F)
  
  layout(matrix(c(2,1,4,3), 2, 2, byrow = TRUE), heights=c(6,1), widths=c(1,6))
  par(mar=c(0,0,1,1))
  image((mat1[,rev(1:ncol(mat1))]),col = colorRampPalette(c("blue", "white", "red"))(20),frame.plot=TRUE,axes=F)
  par(mar=c(0,0.5,1,0.5))
  pc_rev = rev(pc1[pc1$V1==chr[2], 4])
  barplot(pc_rev,  col = ifelse(pc_rev > 0, "darkred", "darkblue"), horiz=TRUE,border = NA,axes=F, space=0, yaxs='i')
  par(mar=c(0.5,0,0.5,1))
  barplot(pc1[pc1$V1==chr[1], 4], col = ifelse(pc1[pc1$V1==chr[1], 4] > 0, "darkred", "darkblue"), border = NA,axes =F, space=0, xaxs='i')
  
}
plot_trans_500K_heatmap_pc1("../hicpro/G1/iced/500000/G1_500000_iced.matrix","../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", c("chr11","chr12"), paste0("../homer/","G1",".PC1.bedGraph"))
plot_trans_500K_heatmap_pc1("../hicpro/ST1/iced/500000/ST1_500000_iced.matrix","../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", c("chr1","chr2"), paste0("../homer/","ST1",".PC1.bedGraph"))
plot_trans_500K_heatmap_pc1("../hicpro/CLP1/iced/500000/CLP1_500000_iced.matrix","../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", c("chr1","chr12"), paste0("../homer/","CLP1",".PC1.bedGraph"))
plot_trans_500K_heatmap_pc1("../hicpro/LT1/iced/500000/LT1_500000_iced.matrix","../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", c("chr1","chr12"), paste0("../homer/","LT1",".PC1.bedGraph"))

for(i in Sys.glob("../hicpro/*/iced/500000/*_500000_iced.matrix")){
  sample = strsplit(basename(i), "_")[[1]][1]
  if(sample=="MPP1"){
    next
  }
  png(paste0("../figure/transInteraction/", sample,".500K_chr1chr2.png"), 640,600, res=300)
  par(mar=c(0,0,0,0))
  plot_trans_500K_heatmap_pc1(i, "../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", c("chr1","chr2"), paste0("../homer/",sample,".PC1.bedGraph"))
  dev.off()
}

# stat trans-telomere interaction index
library(dplyr)
Tans_Inter_Telo_Centro = data.frame()

for (i in Sys.glob("../analysis/matrix/homer_2M_balance/*2M.txt")) {
  sample = strsplit(basename(i), "[.]")[[1]][1]
  
  LT1_trans = data.frame(fread(i, sep="\t", drop=1), row.names = 1, stringsAsFactors = F)
  #LT1_trans = data.frame(fread("../analysis/matrix/homer_2M_balance/MK1.2M.txt", sep="\t", drop=1), row.names = 1, stringsAsFactors = F)
  colnames(LT1_trans) = rownames(LT1_trans)
  q = quantile(na.omit(as.vector(as.matrix(LT1_trans))), 0.9)
  LT1_trans[LT1_trans > 1] <- q
  
  LT1_trans[lower.tri(LT1_trans)] <- NA
  LT1_trans.m = melt(as.matrix(LT1_trans))
  LT1_trans.m = na.omit(LT1_trans.m)
  LT1_trans.m$chr1 = sapply(strsplit(as.character(LT1_trans.m$Var1), "-"), "[", 1)
  LT1_trans.m$chr2 = sapply(strsplit(as.character(LT1_trans.m$Var2), "-"), "[", 1)
  LT1_trans.m = LT1_trans.m[LT1_trans.m$chr1 != LT1_trans.m$chr2,]
  LT1_trans.m = LT1_trans.m[LT1_trans.m$chr1 %in% paste0("chr", c(seq(1:19), "X")) & LT1_trans.m$chr2 %in% paste0("chr", c(seq(1:19), "X")),]
  LT1_trans.m$pos1 = as.numeric(sapply(strsplit(as.character(LT1_trans.m$Var1), "-"), "[", 2))
  LT1_trans.m$pos2 = as.numeric(sapply(strsplit(as.character(LT1_trans.m$Var2), "-"), "[", 2))
  
  chr_length = data.frame(row.names = "len")
  for(j in paste0("chr", c(seq(1:19), "X"))){
    chr_length[[j]] = max(LT1_trans.m[LT1_trans.m$chr1 == j, "pos1"])
  }
  chr_length = as.data.frame(t(chr_length))
  LT1_trans.m$chr1_len = chr_length[match(LT1_trans.m$chr1, rownames(chr_length)),]
  LT1_trans.m$chr2_len = chr_length[match(LT1_trans.m$chr2, rownames(chr_length)),]
  LT1_trans.m$region1 = "other"
  LT1_trans.m$region2 = "other"
  LT1_trans.m[LT1_trans.m$pos1 < (LT1_trans.m$chr1_len *0.2),"region1"] = "Centromere"
  LT1_trans.m[LT1_trans.m$pos2 < (LT1_trans.m$chr2_len *0.2),"region2"] = "Centromere"
  LT1_trans.m[LT1_trans.m$pos1 > (LT1_trans.m$chr1_len *0.8),"region1"] = "Telomere"
  LT1_trans.m[LT1_trans.m$pos2 > (LT1_trans.m$chr2_len *0.8),"region2"] = "Telomere"
  
  LT1_trans.stat = LT1_trans.m  %>% group_by(chr1, chr2,region1,region2) %>%
    summarise(sum=sum(value))
  
  LT1_trans.total = LT1_trans.m %>% group_by(chr1,chr2) %>%
    summarise(total=sum(value)) 
  
  LT1_trans.stat = merge(LT1_trans.stat, LT1_trans.total, by=c("chr1", "chr2"))
  LT1_trans.stat$ratio = LT1_trans.stat$sum / LT1_trans.stat$total

  LT1_trans.stat = LT1_trans.stat[LT1_trans.stat$region1!="other" & LT1_trans.stat$region2!="other",]
  LT1_trans.stat$sample = sample
  
  Tans_Inter_Telo_Centro = rbind(Tans_Inter_Telo_Centro, LT1_trans.stat)
}



Tans_Inter_Telo_Centro$cell = stringr::str_sub(Tans_Inter_Telo_Centro$sample, 1, -2)
Tans_Inter_Telo_Centro$rep = stringr::str_sub(Tans_Inter_Telo_Centro$sample, -1, -1)
Tans_Inter_Telo_Centro$cell = factor(Tans_Inter_Telo_Centro$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))
# remove MPP1
Tans_Inter_Telo_Centro = Tans_Inter_Telo_Centro[Tans_Inter_Telo_Centro$sample!="MPP1",]

# Rabl index
Rabl_index = matrix(nrow=0,ncol=4)
for(i in unique(Tans_Inter_Telo_Centro$sample)){
  for(c1 in paste0("chr", c(seq(1:19), "X"))){
    for(c2 in paste0("chr", c(seq(1:19), "X"))){
      if(c1 != c2){
        ratio = Tans_Inter_Telo_Centro[Tans_Inter_Telo_Centro$sample == i & Tans_Inter_Telo_Centro$chr1 == c1 & Tans_Inter_Telo_Centro$chr2 == c2,]
        ratio = sum(ratio[ratio$region1 == ratio$region2, "sum"]) / sum(ratio[ratio$region1 != ratio$region2, "sum"])
        Rabl_index = rbind(Rabl_index, c(i, c1, c2, ratio))
      }
    }
  }
}
Rabl_index = data.frame(Rabl_index, stringsAsFactors = F)
colnames(Rabl_index) = c("sample", "chr1", "chr2", "Rabl_index")
Rabl_index$cell = stringr::str_sub(Rabl_index$sample, 1, -2)
Rabl_index$rep = stringr::str_sub(Rabl_index$sample, -1, -1)
Rabl_index$cell = factor(Rabl_index$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))
Rabl_index$Rabl_index = as.numeric(Rabl_index$Rabl_index)
ggplot(na.omit(Rabl_index), aes(cell, Rabl_index))+
  #geom_jitter(width=0.20,aes(shape=rep)) +
  stat_boxplot(geom ='errorbar', width = 0.6, coef=1.5) +
  geom_boxplot(outlier.size=-1)+
  theme_classic(base_size = 20) +xlab("")+ylab("Rabl_index")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))


Rabl_index1 = na.omit(Rabl_index) %>% group_by(sample, cell, rep) %>%
  dplyr::summarise(rabl = mean(Rabl_index)) 
ggplot(Rabl_index1, aes(cell, rabl, fill=rep))+
  geom_point(aes(shape=rep), size=4) +
  theme_minimal(base_size = 20) +xlab("")+ylab("Rabl index")+
  theme(axis.text.x = element_text(color="black",angle = 30, hjust = 1),
        axis.text.y = element_text(color="black"))

ggsave("../figure/transInteraction/Rabl_index.pdf", device = "pdf", width = 6, height = 3)

ggplot(Rabl_index[Rabl_index$chr1=="chr1" & Rabl_index$chr2=="chr2",], aes(cell, Rabl_index, fill=rep))+
  geom_point(aes(shape=rep), size=4) +
  theme_minimal(base_size = 20) +xlab("")+ylab("Rabl index")+
  theme(axis.text.x = element_text(color="black",angle = 30, hjust = 1),
        axis.text.y = element_text(color="black"))


# TCI
TCI = Tans_Inter_Telo_Centro[Tans_Inter_Telo_Centro$region1=="Centromere" & Tans_Inter_Telo_Centro$region2=="Centromere", ]

ggplot(TCI, aes(cell, ratio))+
  #geom_jitter(width=0.20,aes(shape=rep)) +
  stat_boxplot(geom ='errorbar', width = 0.6, coef=1.5) +
  geom_boxplot(outlier.size=-1)+
  theme_classic(base_size = 20) +xlab("")+ylab("TCI")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))


TCI = TCI %>% group_by(sample, cell, rep) %>%
  summarise(sums=sum(sum), totals=sum(total)) 

ggplot(TCI, aes(cell, sums/totals, fill=rep))+
  geom_point(aes(shape=rep), size=4) +
  theme_minimal(base_size = 20) +xlab("")+ylab("TCI")+
  theme(axis.text.x = element_text(color="black",angle = 30, hjust = 1))
ggsave("../figure/transInteraction/TCI.pdf", device = "pdf", width = 6, height = 3)

# TTI
TTI = Tans_Inter_Telo_Centro[Tans_Inter_Telo_Centro$region1=="Telomere" & Tans_Inter_Telo_Centro$region2=="Telomere", ]
ggplot(TTI[TTI$chr1!="chrX" & TTI$chr2!="chrX",], aes(cell, ratio))+
  geom_jitter(width=0.20,aes(shape=rep)) +
  stat_boxplot(geom ='errorbar', width = 0.6, coef=1.5) +
  geom_boxplot(outlier.size=-1)+
  theme_classic(base_size = 20) +xlab("")+ylab("TTI")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

TTI = TTI %>% group_by(sample, cell, rep) %>%
  summarise(sums=sum(sum), totals=sum(total)) 

ggplot(TTI, aes(cell, sums/totals, fill=rep))+
  geom_point(aes(shape=rep), size=4) +
  theme_minimal(base_size = 20) +xlab("")+ylab("TTI")+
  theme(axis.text.x = element_text(color="black",angle = 30, hjust = 1))
ggsave("../figure/transInteraction/TTI.pdf", device = "pdf", width = 6, height = 3)


# TCTI
TCTI = Tans_Inter_Telo_Centro[Tans_Inter_Telo_Centro$region1 != Tans_Inter_Telo_Centro$region2, ]
ggplot(TCTI[TCTI$chr1!="chrX" & TCTI$chr2!="chrX",], aes(cell, ratio))+
  geom_jitter(width=0.20,aes(shape=rep)) +
  stat_boxplot(geom ='errorbar', width = 0.6, coef=1.5) +
  geom_boxplot(outlier.size=-1)+
  theme_classic(base_size = 20) +xlab("")+ylab("TTI")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

TCTI = TCTI %>% group_by(sample, cell, rep) %>%
  summarise(sums=sum(sum), totals=sum(total)) 

ggplot(TCTI, aes(cell, sums/totals, fill=rep))+
  geom_point(aes(shape=rep), size=4) +
  theme_minimal(base_size = 20) +xlab("")+ylab("TCTI")+
  theme(axis.text.x = element_text(color="black",angle = 30, hjust = 1))
ggsave("../figure/transInteraction/TCTI.pdf", device = "pdf", width = 6, height = 3)



# Image
T_counts = read.table("../data/image_Telomeres.count.txt", header = 1, sep="\t")
T_counts = na.omit(melt(T_counts))
ggplot(T_counts, aes(variable, value, fill=variable))+
  geom_boxplot() +
  theme_classic(base_size = 25) + xlab("") + ylab("Telomere numbers") +
  scale_fill_d3()+
  theme(axis.text.x = element_text(color="black",angle = 30, hjust = 1))
ggsave("../figure/transInteraction/Image_Telomere_count.pdf", device = "pdf", width = 6, height = 3)

C_counts = read.table("../data/image_Centromeres.count.txt", header = 1, sep="\t")
C_counts = na.omit(melt(C_counts))
ggplot(C_counts, aes(variable, value, fill=variable))+
  geom_boxplot() +
  theme_classic(base_size = 25) + xlab("") + ylab("Centromere numbers") +
  scale_fill_d3()+
  theme(axis.text.x = element_text(color="black",angle = 30, hjust = 1))
ggsave("../figure/transInteraction/Image_Centromere_count.pdf", device = "pdf", width = 6, height = 3)


# interaction between compartment A/B
# AB compartment

Trans_compartment_interaction = function(mat_path, abs_bed, epig1){
  epig1 = read.table(epig1, skip = 1)
  
  abs_bed = read.table(abs_bed)
  #abs_bed = "../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed"
  
  epig1_abs = merge(epig1, abs_bed, by=c("V1","V2"))
  epig1_abs_A = epig1_abs[epig1_abs$V4.x > 0, 6]
  epig1_abs_B = epig1_abs[epig1_abs$V4.x < 0, 6]
  #mat_path = "../hicpro/MPP2/iced/500000/MPP2_500000_iced.matrix"
  
  mat = read.table(mat_path, stringsAsFactors = F)
  mat$chr1 = abs_bed[mat$V1, 1]
  mat$chr2 = abs_bed[mat$V2, 1]
  # remove diag
  mat = mat[mat$V1 != mat$V2,]
  
  #
  mat_cis = mat[mat$chr1 == mat$chr2,]
  mat_cis_mean = mean(mat_cis$V3)
  mat_trans = mat[mat$chr1 != mat$chr2,]
  mat_trans_mean = mean(mat_trans$V3)
  
  mat_cis_AA = mean(mat_cis[mat_cis$V1 %in% epig1_abs_A & mat_cis$V2 %in% epig1_abs_A, 3])
  mat_cis_AB = mean(mat_cis[(mat_cis$V1 %in% epig1_abs_A & mat_cis$V2 %in% epig1_abs_B)|(mat_cis$V1 %in% epig1_abs_B & mat_cis$V2 %in% epig1_abs_A), 3])
  mat_cis_BB = mean(mat_cis[mat_cis$V1 %in% epig1_abs_B & mat_cis$V2 %in% epig1_abs_B, 3])
  
  mat_trans_AA = mean(mat_trans[mat_trans$V1 %in% epig1_abs_A & mat_trans$V2 %in% epig1_abs_A, 3])
  mat_trans_AB = mean(mat_trans[(mat_trans$V1 %in% epig1_abs_A & mat_trans$V2 %in% epig1_abs_B)|(mat_trans$V1 %in% epig1_abs_B & mat_trans$V2 %in% epig1_abs_A), 3])
  mat_trans_BB = mean(mat_trans[mat_trans$V1 %in% epig1_abs_B & mat_trans$V2 %in% epig1_abs_B, 3])
  
  return(c(mat_cis_AA/mat_cis_mean,mat_cis_AB/mat_cis_mean,mat_cis_BB/mat_cis_mean,
           mat_trans_AA/mat_trans_mean,mat_trans_AB/mat_trans_mean,mat_trans_BB/mat_trans_mean))
  
}

Trans_compartment_interaction("../hicpro/CLP1/iced/500000/CLP1_500000_iced.matrix", "../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", "../homer/CLP1.PC1.bedGraph")

Trans_Inter_Compartment = data.frame(row.names = c("Cis_AA","Cis_AB","Cis_BB","Trans_AA","Trans_AB","Trans_BB"))
for (i in Sys.glob("../hicpro/*/iced/500000/*_500000_iced.matrix")) {
  sample = strsplit(basename(i), "_")[[1]][1]
  if(sample=="MPP1"){
    next
  }
  Trans_Inter_Comp = Trans_compartment_interaction(i, "../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", paste0("../homer/",sample,".PC1.bedGraph"))
  Trans_Inter_Compartment[[sample]] = Trans_Inter_Comp
}
Trans_Inter_Compartment = as.data.frame(t(Trans_Inter_Compartment))
Trans_Inter_Compartment.m = melt(as.matrix(Trans_Inter_Compartment))
Trans_Inter_Compartment.m$cell = stringr::str_sub(Trans_Inter_Compartment.m$Var1, 1, -2)
Trans_Inter_Compartment.m$rep = stringr::str_sub(Trans_Inter_Compartment.m$Var1, -1, -1)
Trans_Inter_Compartment.m$cell = factor(Trans_Inter_Compartment.m$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))

ggplot(Trans_Inter_Compartment.m[Trans_Inter_Compartment.m$Var2 %in% c("Trans_AA","Trans_AB","Trans_BB"),], aes(Var2, value, col=cell))+
  geom_jitter(aes(shape=rep), size=4, width = 0.1) +
  theme_classic(base_size = 20) +xlab("")+ylab("Normalized interaction")+
  scale_color_brewer(palette="Paired") +
  theme(axis.text.x = element_text(color="black", angle = 30, hjust = 1))
ggsave("../figure/transInteraction/TransInteraction_compartment.pdf", device = "pdf", width = 6, height = 5.3)


ggplot(Trans_Inter_Compartment.m[Trans_Inter_Compartment.m$Var2 %in% c("Cis_AA","Cis_AB","Cis_BB"),], aes(Var2, value, col=cell))+
  geom_jitter(aes(shape=rep), size=4, width = 0.1) +
  theme_classic(base_size = 20) +xlab("")+ylab("Normalized interaction")+
  scale_color_brewer(palette="Paired") +
  theme(axis.text.x = element_text(color="black", angle = 30, hjust = 1))
ggsave("../figure/transInteraction/CisInteraction_compartment.pdf", device = "pdf", width = 6, height = 5.3)


ggplot(Trans_Inter_Compartment.m, aes(Var2, value, col=cell))+
  geom_jitter(aes(shape=rep), size=4, width = 0.1) +
  theme_minimal(base_size = 20) +xlab("")+ylab("Normalized interaction")+
  scale_color_brewer(palette="Paired") +
  theme(axis.text.x = element_text(color="black",angle = 30, hjust = 1))
ggsave("../figure/transInteraction/Cis_Trans_Interaction_compartment.pdf", device = "pdf", width = 8, height = 5.3)



Trans_Inter_Compartment$tran_AA_BB = Trans_Inter_Compartment$Trans_AA / Trans_Inter_Compartment$Trans_BB

Trans_Inter_Compartment$cell = stringr::str_sub(rownames(Trans_Inter_Compartment), 1, -2)
Trans_Inter_Compartment$rep = stringr::str_sub(rownames(Trans_Inter_Compartment), -1, -1)
Trans_Inter_Compartment$cell = factor(Trans_Inter_Compartment$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))


ggplot(Trans_Inter_Compartment, aes(cell, Trans_AA, col=cell))+
  geom_point(aes(shape=rep), size=4) +
  theme_classic(base_size = 20) +xlab("")+
  scale_color_brewer(palette="Paired") +
  theme(legend.position = "none",
        axis.text.x = element_text(color="black",angle = 30, hjust = 1))
ggsave("../figure/transInteraction/TransInteraction_AA_.pdf", device = "pdf", width = 6, height = 3)

ggplot(Trans_Inter_Compartment, aes(cell, Trans_BB, col=cell))+
  geom_point(aes(shape=rep), size=4) +
  theme_classic(base_size = 20) +xlab("")+
  scale_color_brewer(palette="Paired") +
  theme(legend.position = "none",
        axis.text.x = element_text(color="black",angle = 30, hjust = 1))
ggsave("../figure/transInteraction/TransInteraction_BB_.pdf", device = "pdf", width = 6, height = 3)





