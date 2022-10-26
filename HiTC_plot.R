#####
# plot heatmaps of hic matrix
#

library(RColorBrewer)

hic_heatmap <- function(mat_path, abs_bed, chrom, epig1){
  
  abs_bed = read.table(abs_bed)
  chr1_abs = abs_bed[abs_bed$V1 == chrom, 4]
  
  mat = read.table(mat_path, stringsAsFactors = F)
  mat_chr1 = mat[mat$V1 %in% chr1_abs & mat$V2 %in% chr1_abs,]
  mat_chr1 = dcast(mat_chr1, V1 ~ V2, value.var="V3")
  rownames(mat_chr1) = mat_chr1$V1
  mat_chr1 = mat_chr1[,-1]
  
  # symmetric matrix
  mat_chr1[lower.tri(mat_chr1)] = t(mat_chr1)[lower.tri(mat_chr1)]
  mat_chr1 = as.matrix(mat_chr1)
  
  mat_chr1[is.na(mat_chr1)] <- 0
  qt = quantile(mat_chr1, probs = 0.95)
  mat_chr1[mat_chr1 > qt] <- qt
  # hic matrix normalized by quantile
  mat_chr1 = mat_chr1/qt
  
  
  epig1 = read.table(epig1, skip = 1)
  epig1 = epig1[epig1$V1 == chrom,]
  
  #plot
  layout(matrix(c(1,2), 2, 1, byrow = TRUE), heights=c(4,1), widths=c(1,1))
  par(mar=c(0,1,1,1))
  image((mat_chr1[,rev(colnames(mat_chr1))]), col = colorRampPalette(c("white", "red"))(n = 20),frame.plot=TRUE,axes=F,cex=1)
  box(lwd = 1)
  par(mar=c(1,1,0,1))
  barplot(epig1$V4, col = ifelse(epig1$V4 > 0, "darkred", "darkblue"), border = NA,axes =F, space=0, xaxs='i') 
}
hic_heatmap("../hicpro/MK2/iced/500000/MK2_500000_iced.matrix", abs_bed, "chr2","../homer/MK2.PC1.bedGraph")

for (i in c("LT", "ST", "CLP","CMP","GMP","MEP","MKP","MK","G")){
  mat_path = paste0("../hicpro/", i,"1/iced/500000/", i,"1_500000_iced.matrix")
  chrom = "chr1"
  abs_bed = "../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed"
  epig1 = paste0("../homer/", i, "1.PC1.bedGraph")
  png(paste0("../figure/heatmap/", i, "_chr1.hp.png"), res=300, width = 900, height = 1000)
  hic_heatmap(mat_path, abs_bed, chrom,epig1)
  dev.off()
}


mat_path = "../hicpro/MPP2/iced/500000/MPP2_500000_iced.matrix"
chrom = "chr1"
abs_bed = "../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed"
epig1 = "../homer/MPP2.PC1.bedGraph"
png(("../figure/heatmap/MPP2_chr1.hp.png"), res=300, width = 900, height = 1000)
hic_heatmap(mat_path, abs_bed, chrom,epig1)
dev.off()





##
mat1_path = "../hicpro/G1/iced/500000/G1_500000_iced.matrix"
mat2_path = "../hicpro/LT1/iced/500000/LT1_500000_iced.matrix"
abs_bed = "../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed"
hic_heatmap_sub <- function(mat1_path, mat2_path, abs_bed, chrom){
  
  abs_bed = read.table(abs_bed)
  chr1_abs = abs_bed[abs_bed$V1 == chrom, 4]
  
  mat1 = read.table(mat1_path, stringsAsFactors = F)
  mat1 = mat1[mat1$V1 %in% chr1_abs & mat1$V2 %in% chr1_abs,]
  mat1 = dcast(mat1, V1 ~ V2, value.var="V3")
  rownames(mat1) = mat1$V1
  mat1 = mat1[,-1]
  
  # symmetric matrix
  mat1[lower.tri(mat1)] = t(mat1)[lower.tri(mat1)]
  mat1 = as.matrix(mat1)
  
  mat1[is.na(mat1)] <- 0
  diag(mat1) <- 0
  
  # hic matrix normalized by quantile
  sm = sum(mat1)
  qt = quantile(mat1, probs = 0.5)
  mat1 = mat1/sm

  
  
  mat2 = read.table(mat2_path, stringsAsFactors = F)
  mat2 = mat2[mat2$V1 %in% chr1_abs & mat2$V2 %in% chr1_abs,]
  mat2 = dcast(mat2, V1 ~ V2, value.var="V3")
  rownames(mat2) = mat2$V1
  mat2 = mat2[,-1]
  
  # symmetric matrix
  mat2[lower.tri(mat2)] = t(mat2)[lower.tri(mat2)]
  mat2 = as.matrix(mat2)
  
  mat2[is.na(mat2)] <- 0
  diag(mat2) <- 0
  
  # hic matrix normalized by quantile
  sm = sum(mat2)
  qt = quantile(mat2, probs = 0.5)
  mat2 =  mat2/sm
  
  # substract
  mat_sub = mat1 - mat2
  qt = quantile(mat_sub, probs = c(0.05,0.95))
  mat_sub[mat_sub > qt[2]] <- qt[2]
  mat_sub[mat_sub < qt[1]] <- qt[1]
  par(mar=c(1,1,1,1))
  image((mat_sub[,rev(colnames(mat_sub))]), col = colorRampPalette(c("blue", "white","red"))(n = 20),frame.plot=TRUE,axes=F)

}

abs_bed = "../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed"

png(("../figure/heatmap/G-LT_chr1.hp.png"), res=300, width = 900, height = 800)
hic_heatmap_sub("../hicpro/G1/iced/500000/G1_500000_iced.matrix", "../hicpro/LT1/iced/500000/LT1_500000_iced.matrix", abs_bed, chrom)
dev.off()
png(("../figure/heatmap/MK-LT_chr1.hp.png"), res=300, width = 900, height = 800)
hic_heatmap_sub("../hicpro/MK1/iced/500000/MK1_500000_iced.matrix", "../hicpro/LT1/iced/500000/LT1_500000_iced.matrix", abs_bed, chrom)
dev.off()
png(("../figure/heatmap/MEP-LT_chr1.hp.png"), res=300, width = 900, height = 800)
hic_heatmap_sub("../hicpro/MEP1/iced/500000/MEP1_500000_iced.matrix", "../hicpro/LT1/iced/500000/LT1_500000_iced.matrix", abs_bed, chrom)
dev.off()
png(("../figure/heatmap/CLP-LT_chr1.hp.png"), res=300, width = 900, height = 800)
hic_heatmap_sub("../hicpro/CLP1/iced/500000/CLP1_500000_iced.matrix", "../hicpro/LT1/iced/500000/LT1_500000_iced.matrix", abs_bed, chrom)
dev.off()
png(("../figure/heatmap/ST-LT_chr1.hp.png"), res=300, width = 900, height = 800)
hic_heatmap_sub("../hicpro/ST1/iced/500000/ST1_500000_iced.matrix", "../hicpro/LT1/iced/500000/LT1_500000_iced.matrix", abs_bed, chrom)
dev.off()

png(("../figure/heatmap/MPP-LT.hp.png"), res=300, width = 900, height = 800)
abs_bed = "../hicpro/CLP1/raw/1000000/CLP1_1000000_abs.bed"
hic_heatmap_sub("../hicpro/MPP2/iced/1000000/MPP2_1000000_iced.matrix", "../hicpro/LT1/iced/1000000/LT1_1000000_iced.matrix", abs_bed, chrom)
dev.off()




######
# Start End Index
library(data.table)
SEI = data.frame(row.names = c("sample","chr", "SEI"))
for (i in Sys.glob("../analysis/matrix/homer_500k_balance/*500k.txt")) {
  mat = data.frame(fread(i, sep="\t", drop=1), row.names = 1)
  diag(mat) <- 0 
  bn = strsplit(basename(i), "[.]")[[1]]
  sample = bn[1]
  chr = bn[2]
  mat_dim = ncol(mat)
  l = as.numeric(mat_dim*0.1)
  SEindex = mat[1:l, (mat_dim-l+1):mat_dim]
  SEindex = sum(SEindex)*2 / sum(mat)
  SEI[[basename(i)]] = c(sample, chr, SEindex)
}

SEI.t = data.frame(t(SEI), stringsAsFactors = F)
SEI.t$chr = factor(SEI.t$chr, levels = paste0("chr", c(seq(1:19), "X")))
SEI.t$SEI = as.numeric(SEI.t$SEI)
SEI.t$cell = stringr::str_sub(SEI.t$sample, 1, -2)
SEI.t$rep = stringr::str_sub(SEI.t$sample, -1, -1)
SEI.t$cell = factor(SEI.t$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))
ggplot(SEI.t[SEI.t$chr %in%  paste0("chr", c(seq(1:19))), ], aes(cell, SEI))+
  geom_jitter(width=0.20,aes(col=chr, shape=rep)) +
  geom_boxplot(outlier.size=-1, fill=NA)+
  theme_minimal(base_size = 20) +
  scale_colour_manual(values=rainbow(25)[1:20])+xlab("")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("../figure/SEindex.pdf", device = "pdf", width = 6, height = 3)






hic_heatmap_region <- function(mat_path, abs_bed, region, epig1){
  abs_bed = read.table(abs_bed)
  abs_bed$V2= as.numeric(abs_bed$V2)
  abs_bed$V3= as.numeric(abs_bed$V3)
  #chr1_abs = abs_bed[abs_bed$V1 == region[1], 4]
  chr1_abs = abs_bed[abs_bed$V1 == region[1] & abs_bed$V2 >= as.numeric(region[2]) & abs_bed$V3 <= as.numeric(region[3]), 4]
  
  mat = read.table(mat_path, stringsAsFactors = F)
  mat_chr1 = mat[mat$V1 %in% chr1_abs & mat$V2 %in% chr1_abs,]
  mat_chr1 = dcast(mat_chr1, V1 ~ V2, value.var="V3")
  rownames(mat_chr1) = mat_chr1$V1
  mat_chr1 = mat_chr1[,-1]
  
  # symmetric matrix
  mat_chr1[lower.tri(mat_chr1)] = t(mat_chr1)[lower.tri(mat_chr1)]
  mat_chr1 = as.matrix(mat_chr1)
  
  mat_chr1[is.na(mat_chr1)] <- 0
  qt = quantile(mat_chr1, probs = 0.95)
  mat_chr1[mat_chr1 > qt] <- qt
  # hic matrix normalized by quantile
  mat_chr1 = mat_chr1/qt
  
  
  epig1 = read.table(epig1, skip = 1)
  epig1 = epig1[epig1$V1 == region[1],]
  epig1 = epig1[epig1$V2 >= as.numeric(region[2]) & epig1$V3 <= as.numeric(region[3]),]
  
  #plot
  layout(matrix(c(1,2), 2, 1, byrow = TRUE), heights=c(4,1), widths=c(1,1))
  par(mar=c(0,1,1,1))
  image(log(mat_chr1[,rev(colnames(mat_chr1))]), col = colorRampPalette(c("blue","white", "red"))(n = 20),frame.plot=TRUE,axes=F,cex=1)
  box(lwd = 1)
  par(mar=c(1,1,0,1))
  barplot(epig1$V4, col = ifelse(epig1$V4 > 0, "darkred", "darkblue"), border = NA,axes =F, space=0, xaxs='i') 
}



# GR specific
abs_bed = "../hicpro/CLP1/raw/100000/CLP1_100000_abs.bed"
png(paste0("../figure/heatmap/G_chr17_Fpr1.hp.png"), res=300, width = 900, height = 1000)
hic_heatmap_region("../hicpro/G2/iced/100000/G2_100000_iced.matrix", abs_bed, c("chr17",15000000,25000000),"../homer/G2.PC1.bedGraph")
dev.off()
png(paste0("../figure/heatmap/GMP_chr17_Fpr1.hp.png"), res=300, width = 900, height = 1000)
hic_heatmap_region("../hicpro/GMP2/iced/100000/GMP2_100000_iced.matrix", abs_bed, c("chr17",15000000,25000000),"../homer/GMP2.PC1.bedGraph")
dev.off()
png(paste0("../figure/heatmap/CMP_chr17_Fpr1.hp.png"), res=300, width = 900, height = 1000)
hic_heatmap_region("../hicpro/CMP1/iced/100000/CMP1_100000_iced.matrix", abs_bed, c("chr17",15000000,25000000),"../homer/CMP1.PC1.bedGraph")
dev.off()



