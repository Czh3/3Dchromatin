setwd("/lustre/user/liclab/zhangc/proj/blood3D_3/script/")

require(HiTC)
library(GenomicRanges)
library(rtracklayer)
library(wesanderson)
library("ggsci")

gtfRangedData<- import.gff("/lustre/user/liclab/publicData/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf")

#### homer
# the pc1 value called by homer maybe inverse
ABcomp_500k = read.table("../homer/LT1.PC1.bedGraph", skip = 1)[,1:3]
rownames(ABcomp_500k) = paste(ABcomp_500k$V1, ABcomp_500k$V2, sep = "_")
for (i in Sys.glob("../homer/*.PC1.bedGraph")){
  ABcomp = read.table(i, skip = 1)
  rownames(ABcomp) = paste(ABcomp$V1, ABcomp$V2, sep = "_")
  sample = strsplit(basename(i), "[.]")[[1]][1]
  ABcomp_500k[[sample]] = ABcomp[rownames(ABcomp_500k), 4]
}

ABcomp_500k = ABcomp_500k[,!colnames(ABcomp_500k) %in% c("MPP1")]
ABcomp_500k = na.omit(ABcomp_500k)

pca = prcomp(t(ABcomp_500k[,-c(1:3)]), scale = F, center = T)
pca_plot = as.data.frame(pca$x)
pca_plot$cell = toupper(gsub("2","", gsub("1", "", rownames(pca_plot))))
pca_plot$cell = factor(pca_plot$cell, level = c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))



pdf("../figure/AB_PC1_PCA_500k.pdf", 5, 4)
ggplot(pca_plot, aes(PC1, PC2, color=cell)) +
  geom_point(size=3) +
  theme_classic(base_size=20) +
  scale_color_d3()
dev.off()


pc1.cor = cor(ABcomp_500k[,-c(1:3)], method = "spearman")
pheatmap(pc1.cor, scale="none", border_color=NA,  color = rev(colorRampPalette(c("red", "white", "blue"))(20)))


pdf("../figure/AB_cor_heatmap_500k.pdf", 5.3, 5)
pheatmap(pc1.cor, scale="none", border_color=NA,  color = rev(colorRampPalette(c("red", "white", "blue"))(20)))
dev.off()


my_line <- function(x,y,...){
  smoothScatter(x, y, nrpoints = 0, add = TRUE)
}
panel.cor <- function(x, y, digits = 1, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  beta = lm(x ~ y)$coefficients[2]
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0("r=", txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1)
}
png("../figure/AB_rep_scatter_500k.png",2800,2800, res=300)
pairs(ABcomp_500k[,-c(1:3)], lower.panel = my_line, upper.panel = panel.cor, xlim=c(-0.1,0.1), ylim=c(-0.1,0.1))
dev.off()

# AB proportion
library(reshape2)
library(ggplot2)

ABcomp_binary = ifelse(ABcomp_500k[,-c(1:3)] > 0, 1, -1 )
rownames(ABcomp_binary) = paste(ABcomp_500k$V1, ABcomp_500k$V2, sep="_")
ABprop = ifelse(ABcomp_binary>0, "A", "B")  
ABprop = apply(ABprop, 2, table)
ABprop.m = melt(ABprop)
ABprop.m$Var2 = factor(ABprop.m$Var2, levels = c("LT1","LT2","ST1","ST2","MPP2","CLP1","CLP2",
                                                 "CMP1","CMP2","MEP1","MEP2","GMP1","GMP2","MKP1","MKP2","MK1","MK2","G1","G2"))
pdf("../figure/AB_proportion_500k.pdf", 6, 3)
ggplot(ABprop.m, aes(Var2, value, fill = Var1)) +
  geom_bar(stat="identity", position = "fill") +
  theme_classic(base_size = 20) + xlab("") + ylab("Proportion") +
  scale_fill_manual( values = c("red","darkblue"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black", angle = 90, hjust = 1,vjust = 0.5))
dev.off()

# AB length
ABcomp_length = matrix(nrow=0,ncol=19)
for (i in 2:nrow(ABcomp_500k)){
  ABcomp_length = rbind(ABcomp_length, ABcomp_500k[i,] * ABcomp_500k[i-1,] < 0)
}
ABcomp_length_stat = colSums(ABcomp_length)
ABcomp_length_stat = nrow(ABcomp_500k)/ABcomp_length_stat * 0.5
ABcomp_length_stat = as.data.frame(ABcomp_length_stat)
ABcomp_length_stat$sample = rownames(ABcomp_length_stat)
ABcomp_length_stat$cell = stringr::str_sub(rownames(ABcomp_length_stat), 1, -2)
ABcomp_length_stat$rep = stringr::str_sub(rownames(ABcomp_length_stat), -1, -1)
ABcomp_length_stat$cell = factor(ABcomp_length_stat$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))
ggplot(ABcomp_length_stat, aes(cell, ABcomp_length_stat, col=cell))+
  geom_point(aes(shape=rep), size=4) +
  theme_classic(base_size = 20) +xlab("")+ylab("compartment length (Mb)")+
  scale_color_d3() +
  theme(legend.position = "none",
        axis.text.x = element_text(color="black",angle = 30, hjust = 1))
ggsave("../figure/compartment_length.pdf", device = "pdf", width = 5, height = 4)


# AB gene expression
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = as.matrix(expr[,c(12,7,19,20,2,17,18,16,25)])
colnames(expr) = c("LT","ST","MPP","CMP","CLP","MEP","GMP","MK","G")

pc_gene_plot = matrix(0, nrow = 0, ncol = 3)
gene_length = "/lustre/user/liclab/zhangc/Taolab/guan/rna-seq/reference/gene_position.protein_coding.bed"
for (i in Sys.glob("../homer/*.PC1.bedGraph")){

  cmd = paste0("grep -v track ",i," | tr  ' ' '\t' |/lustre/user/liclab/biotools/bedtools2/bin/intersectBed -a ",gene_length," -b - -wa -wb -f 0.5 > ../data/temp.pc_gene")
  system(cmd)
  pc_gene = read.table("../data/temp.pc_gene")
  pc_gene$PC = ifelse(pc_gene$V10 >0, "A", "B")
  

  sample = gsub("1|2","", strsplit(basename(i), "[.]")[[1]][1])
  sample_rep = strsplit(basename(i), "[.]")[[1]][1]
  if(sample_rep %in% c("MKP1", "MKP2", "MPP1")){
    next
  }
  pc_gene = cbind(pc_gene, expr[match(pc_gene$V4, rownames(expr)),])
  pc_gene = pc_gene[,c("PC", sample)]
  pc_gene$sample = sample_rep
  colnames(pc_gene) = c("PC", "expr", "sample")
  pc_gene_plot = rbind(pc_gene_plot, pc_gene)
}
system("rm ../data/temp.pc_gene")
pc_gene_plot$sample = factor(pc_gene_plot$sample,
                                levels=c("LT1","LT2","ST1","ST2","MPP1","MPP2","CLP1","CLP2",
                                         "CMP1","CMP2","MEP1","MEP2","GMP1","GMP2","MK1","MK2","G1","G2"))
ggplot(pc_gene_plot, aes(sample, log(expr+1), fill=PC))+
  geom_boxplot(outlier.size=-1, size=1)+ ylim(0,9)+
  theme_classic(base_size = 20) +xlab("")+ylab("log(TPM+1)")+
  theme( axis.text = element_text(color="black"),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("../figure/AB_comp_expr.boxplot.pdf", device = "pdf", width = 12, height = 3.5)


# AB switch
ABcomp_500k = ABcomp_500k[,c("LT1","LT2","ST1","ST2","MPP2","CLP1","CLP2",
                                 "CMP1","CMP2","MEP1","MEP2","GMP1","GMP2","MKP1","MKP2","MK1","MK2","G1","G2")]

ABcomp_binary = ifelse(ABcomp_500k > 0, 1, -1 )

ABcomp_binary_switch = ABcomp_binary[abs(rowSums(ABcomp_binary)) != ncol(ABcomp_binary),]
pheatmap(ABcomp_binary_switch, scale="none", border_color=NA, cluster_cols = F,
         show_rownames = F, color = colorRampPalette(c("purple", "black", "yellow"))(20))

q=quantile(as.matrix(ABcomp_500k), probs=c(0.1,0.9))
ABcomp_500k[ABcomp_500k > q[2]] <- q[2]
ABcomp_500k[ABcomp_500k < q[1]] <- q[1]
png("../figure/AB_PC1_switch_heatmap_500k.png", width = 1200, height = 1000, res=300)
pheatmap(ABcomp_500k[rownames(ABcomp_binary_switch),], scale="none", border_color=NA, cluster_cols = F,
         show_rownames = F, color = colorRampPalette(c("purple", "black", "yellow"))(20))
dev.off()
ph=pheatmap(ABcomp_500k[rownames(ABcomp_binary_switch),], scale="none", border_color=NA, cluster_cols = F,
            show_rownames = F, color = colorRampPalette(c("purple", "black", "yellow"))(20))
write.table(ph$tree_row$labels[ph$tree_row$order], "../data/AB_dynamic.regions.bed", row.names = F, quote = F)

# H3K4me1
AB_H3K4me1 = read.table("../data/AB_dynamic.regions.H3K4me1.bed", stringsAsFactors = F)
rownames(AB_H3K4me1) = paste(AB_H3K4me1$V1, AB_H3K4me1$V2, sep = "_")
AB_H3K4me1 = AB_H3K4me1[,-c(1:3)]
colnames(AB_H3K4me1) = c("LT","ST","MPP","CMP","GMP","MEP","CLP","G")
AB_H3K4me1 = AB_H3K4me1[, c("LT","ST","MPP", "CLP","CMP","MEP","GMP","G")]
AB_H3K4me1 = AB_H3K4me1[ph$tree_row$labels[ph$tree_row$order],]
AB_H3K4me1 = as.matrix(AB_H3K4me1)
class(AB_H3K4me1) <- "numeric"
AB_H3K4me1 = na.omit(AB_H3K4me1)
colm = apply(AB_H3K4me1, 2, median)
AB_H3K4me1 = t(t(AB_H3K4me1)/colm)
AB_H3K4me1 = log10(AB_H3K4me1+0.01)
AB_H3K4me1[AB_H3K4me1 > 1.2] <- 1.2
AB_H3K4me1[AB_H3K4me1 < -1] <- -1
png("../figure/AB_PC1_switch_H3K4me1_heatmap_500k.png", width = 1000, height = 1000, res=300)
pheatmap((AB_H3K4me1), scale="none", border_color=NA, cluster_rows = F,cluster_cols = F,
            show_rownames = F, color = colorRampPalette(c("blue", "white", "red"))(20))
dev.off()
# H3K27ac
AB_H3K27ac = read.table("../data/AB_dynamic.regions.H3K27ac.bed", stringsAsFactors = F)
rownames(AB_H3K27ac) = paste(AB_H3K27ac$V1, AB_H3K27ac$V2, sep = "_")
AB_H3K27ac = AB_H3K27ac[,-c(1:3)]
colnames(AB_H3K27ac) = c("LT","ST","MPP","CMP","GMP","MEP","EryA","EryB","G","Mono","MF","B","CD4","CD8","NK","CLP")
AB_H3K27ac = AB_H3K27ac[, c("LT","ST","MPP", "CLP","CMP","MEP","GMP","G")]
AB_H3K27ac = AB_H3K27ac[ph$tree_row$labels[ph$tree_row$order],]

AB_H3K27ac = as.matrix(AB_H3K27ac)
AB_H3K27ac = na.omit(AB_H3K27ac)
class(AB_H3K27ac) <- "numeric"
AB_H3K27ac = na.omit(AB_H3K27ac)
colm = apply(AB_H3K27ac, 2, median)
AB_H3K27ac = t(t(AB_H3K27ac)/colm)
#AB_H3K27ac = log10(AB_H3K27ac+0.01)
AB_H3K27ac[AB_H3K27ac > 3] <- 3
AB_H3K27ac[AB_H3K27ac < 5] <- 5
#png("../figure/AB_PC1_switch_H3K27ac_heatmap_500k.png", width = 1000, height = 1000, res=300)
pheatmap((AB_H3K27ac), scale="none", border_color=NA, cluster_rows = F,cluster_cols = F,
         show_rownames = F, color = colorRampPalette(c("blue", "white", "red"))(20))
#dev.off()


##### cell type specific A
ABcomp_binary_cellType = (ABcomp_binary[,c(1,3,5,seq(6,19,2))] * ABcomp_binary[,c(1,3,4,seq(6,19,2))+1])
ABcomp_binary_cellType = ABcomp_binary_cellType[rowSums(ABcomp_binary_cellType) == 10,]
ABcomp_binary_cellType = ABcomp_binary[rownames(ABcomp_binary_cellType),c(1,3,5,seq(6,19,2)) ]
ABcomp_binary_cellType = ABcomp_binary_cellType[abs(rowSums(ABcomp_binary_cellType)) != 10,]
ABcomp_binary_cellType = as.data.frame(ABcomp_binary_cellType)
attach(ABcomp_binary_cellType)
ABcomp_binary_cellType = ABcomp_binary_cellType[order(LT1, ST1,MPP2, CLP1,CMP1, MEP1,GMP1,MKP1,MK1,G1, decreasing = T),]
detach(ABcomp_binary_cellType)
pheatmap(ABcomp_500k[rownames(ABcomp_binary_cellType),], scale="none", border_color=NA, cluster_cols = F,cluster_rows = T, show_rownames= F,
         clustering_method = "ward.D2",  color = colorRampPalette(c("purple", "black", "yellow"))(20))

breaksList = seq(-1, 1, by = 0.1)
png("../figure/AB_PC1_switch_heatmap_500k.png", width = 900, height = 800, res=300)
pheatmap(ABcomp_500k[rownames(ABcomp_binary_cellType),], scale="none", border_color=NA, cluster_cols = F,
         clustering_method = "ward.D2", cluster_rows = T, show_rownames= F,
         col=colorRampPalette(c("purple", "black", "yellow"))(length(breaksList)),breaks = breaksList)
dev.off()
ph=pheatmap(ABcomp_500k[rownames(ABcomp_binary_cellType),], scale="none", border_color=NA, cluster_cols = F,
            clustering_method = "ward.D2", cluster_rows = T, show_rownames= F,
            col=colorRampPalette(c("purple", "black", "yellow"))(length(breaksList)),breaks = breaksList)

pdf("../figure/AB_PC1_switch_heatmap_500k.pdf", 5, 4)
ph
dev.off()
write.table(ph$tree_row$labels[ph$tree_row$order], "../data/AB_dynamic.regions.bed", row.names = F, quote = F)



if(FALSE){
ABcomp_binary_cellType = ABcomp_500k[rowSums(ABcomp_binary[,c(1,3,seq(6,19,2))] * ABcomp_binary[,c(1,3,seq(6,19,2))+1]) == 9,]
ABcomp_binary_cellType_rep = (ABcomp_binary_cellType[,c(1,3,5,seq(6,19,2))] +ABcomp_binary_cellType[,c(1,3,5,seq(6,19,2))+1])/2
ABcomp_binary_cellType_region = c()
for (i in 1:ncol(ABcomp_binary_cellType_rep)){
  AB_mean = ABcomp_binary_cellType_rep[,i] - rowMeans(ABcomp_binary_cellType_rep[,-i])
  AB_mean = AB_mean[AB_mean>=0.8]
  AB_mean1 = AB_mean[!names(AB_mean) %in% names(ABcomp_binary_cellType_region)]
  ABcomp_binary_cellType_region = ABcomp_binary_cellType_region[!names(ABcomp_binary_cellType_region)%in%names(AB_mean)]
  ABcomp_binary_cellType_region = c(ABcomp_binary_cellType_region, AB_mean1)
}
png("../figure/AB_PC1_switch_heatmap_500k.png", width = 1200, height = 1000, res=300)
pheatmap(ABcomp_500k[names(ABcomp_binary_cellType_region),], scale="none", border_color=NA, cluster_cols = F,cluster_rows = T,
         clustering_method = "ward.D", show_rownames = F, color = colorRampPalette(c("purple", "black", "yellow"))(20))
dev.off()
ph=pheatmap(ABcomp_500k[names(ABcomp_binary_cellType_region),], scale="none", border_color=NA, cluster_cols = F,cluster_rows = T,
            clustering_method = "ward.D", show_rownames = F, color = colorRampPalette(c("purple", "black", "yellow"))(20))

write.table(ph$tree_row$labels[ph$tree_row$order], "../data/AB_dynamic.regions.bed", row.names = F, quote = F)


pheatmap(ABcomp_binary[names(ABcomp_binary_cellType_region),], scale="none", border_color=NA, cluster_cols = F,cluster_rows = F,
         clustering_method = "ward.D",cutree_rows = 6, show_rownames = F, color = colorRampPalette(c("purple", "black", "yellow"))(20))
}

# H3K4me1
AB_H3K4me1 = read.table("../data/AB_dynamic.regions.H3K4me1.bed", stringsAsFactors = F)
rownames(AB_H3K4me1) = paste(AB_H3K4me1$V1, AB_H3K4me1$V2, sep = "_")
AB_H3K4me1 = AB_H3K4me1[,-c(1:3)]
colnames(AB_H3K4me1) = c("LT","ST","MPP","CMP","GMP","MEP","CLP","G")
AB_H3K4me1 = AB_H3K4me1[, c("LT","ST","MPP", "CLP","CMP","MEP","GMP","G")]
AB_H3K4me1 = AB_H3K4me1[ph$tree_row$labels[ph$tree_row$order],]
AB_H3K4me1 = as.matrix(AB_H3K4me1)
class(AB_H3K4me1) <- "numeric"
AB_H3K4me1 = na.omit(AB_H3K4me1)

colm = apply(AB_H3K4me1, 2, function(x) quantile(x, prob=0.1))
AB_H3K4me1 = t(t(AB_H3K4me1)/colm)


#AB_H3K4me1 = log10(AB_H3K4me1+0.01)
breaksList = seq(0, 3, by = 0.1)
png("../figure/AB_PC1_switch_H3K4me1_heatmap_500k.png", width = 700, height = 800, res=300)
pheatmap(log(AB_H3K4me1), scale="none", border_color=NA, cluster_rows = F,cluster_cols = F,
         show_rownames = F, color = colorRampPalette(c("blue","white", "red"))(length(breaksList)),breaks = breaksList)
dev.off()


# stat cell type specifc
ABcomp_binary_cellType = (ABcomp_binary[,c(1,3,5,seq(6,19,2))] * ABcomp_binary[,c(1,3,4,seq(6,19,2))+1])
ABcomp_binary_cellType = ABcomp_binary_cellType[rowSums(ABcomp_binary_cellType) == 10,]
ABcomp_binary_cellType = ABcomp_binary[rownames(ABcomp_binary_cellType),c(1,3,5,seq(6,19,2)) ]
ABcomp_binary_cellType = as.data.frame(ABcomp_binary_cellType)
ABcomp_binary_cellType_conserve =  table(rowSums(ABcomp_binary_cellType))
ABcomp_binary_cellType_conserve = as.data.frame(ABcomp_binary_cellType_conserve)
ggplot(ABcomp_binary_cellType_conserve, aes(Var1, Freq)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size=20) + ylab("number of 500 Kb bins")+
  theme(legend.position = "none",
      axis.text.x = element_text(color="black",angle = 30, hjust = 1))
ggsave("../figure/ABcompartment_distribution.pdf", device = "pdf", width = 6, height = 4)


for(i in 1:9){
  print(nrow(ABcomp_binary_cellType[ABcomp_binary_cellType[,i] == 1 & ABcomp_binary_cellType[,i+1] == -1, ]))
}




# GR specific
clust5 = cutree(ph$tree_row, k=5)[ph$tree_row[["order"]]]
clust5 = as.data.frame(clust5)
pheatmap(ABcomp_500k[rownames(clust5)[clust5$clust5=="5"],], scale="none", border_color=NA, cluster_cols = F,
         clustering_method = "ward.D2", cluster_rows = T, show_rownames= F,
         col=colorRampPalette(c("purple", "black", "yellow"))(length(breaksList)),breaks = breaksList)

AB_H3K4me1_c5 = AB_H3K4me1[rownames(clust5)[clust5$clust5=="5"],]
AB_H3K4me1_c5.m = melt(AB_H3K4me1_c5)
AB_H3K4me1_c5.m$Var1 = as.character(AB_H3K4me1_c5.m$Var1)
AB_H3K4me1_c5.m$Var2 = factor(AB_H3K4me1_c5.m$Var2, levels = c("LT","ST","MPP", "CLP","CMP","GMP", "MEP","MK","G"))
AB_H3K4me1_c5.m = rbind(AB_H3K4me1_c5.m,c("x","MK",0))
AB_H3K4me1_c5.m$value = as.numeric(AB_H3K4me1_c5.m$value)
library("scales") 
ggplot(AB_H3K4me1_c5.m, aes(Var2, log(value+1), fill=Var2))+
  #stat_boxplot(geom ='errorbar', width = 0.6, coef=1.5) +
  geom_boxplot(outlier.size=-1)+
  theme_classic(base_size = 20) +xlab("")+ylab("log(TPM+1)")+
  scale_fill_manual(values = pal_d3("category10")(10)[-8])+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("../figure/AB_Granulocyte_specifcA_cluster5.H3K4me1.pdf", device = "pdf", width = 5, height = 4)

t.test(AB_H3K4me1_c5.m[AB_H3K4me1_c5.m$Var2 == "G",3], AB_H3K4me1_c5.m[AB_H3K4me1_c5.m$Var2 != "G",3])
# expr
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(2,3,5:7,9,12,14,16:20,25,27,28,31,32)]
expr = expr[rowMeans(expr)>=1,]

AB_switch_expr = read.table("../data/AB_dynamic.regions.RNAexpr.bed", stringsAsFactors = F, sep="\t")
AB_switch_expr = AB_switch_expr[AB_switch_expr$V7 == "protein_coding",]
AB_switch_expr$name = paste(AB_switch_expr$V1, AB_switch_expr$V2, sep = "_")
AB_switch_expr_c5 = AB_switch_expr[AB_switch_expr$name %in% rownames(clust5)[clust5$clust5=="5"],]
AB_switch_expr = na.omit(expr[unique(AB_switch_expr_c5$V8),])

AB_switch_expr.m = melt(AB_switch_expr)
AB_switch_expr.m$cell = toupper(gsub(".2","", gsub(".1", "", AB_switch_expr.m$variable)))
AB_switch_expr.m$cell = factor(AB_switch_expr.m$cell, level = c("LT.HSC", "ST.HSC", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","GR"))
ggplot(AB_switch_expr.m, aes(cell, log(value+1), fill=cell))+
  #stat_boxplot(geom ='errorbar', width = 0.6, coef=1.5) +
  geom_boxplot(outlier.size=-1)+
  theme_classic(base_size = 20) +xlab("")+ylab("log(TPM+1)")+
  scale_fill_manual(values = pal_d3("category10")(10)[-8])+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("../figure/AB_Granulocyte_specifcA_cluster5.expr.pdf", device = "pdf", width = 5, height = 4)

t.test(AB_switch_expr.m[AB_switch_expr.m$cell == "GR",2], AB_switch_expr.m[AB_switch_expr.m$cell != "GR",2])

write.table(AB_switch_expr, "../data/AB_Granu.TMP1.bed", row.names = T, quote = F, sep="\t")

# GO
library(goseq)
library(org.Mm.eg.db)

genes = rep(0, nrow(expr))
names(genes) = rownames(expr)
genes[names(genes) %in% rownames(AB_switch_expr)] = 1
pwf=nullp(genes,"mm10","geneSymbol", plot.fit = FALSE)
GO_G=goseq(pwf,"mm10","geneSymbol")
GO_G = GO_G[GO_G$ontology == "BP",][6:1,]
GO_G$term = factor(GO_G$term, levels = GO_G$term)
ggplot(GO_G, aes(term,-log10(over_represented_pvalue) )) +
  geom_bar(stat = "identity", fill="black") + coord_flip()+
  theme_minimal() + xlab("") + ylab("-log10(p-value)") +
  theme(axis.text = element_text(size=10, color="black"))+
  scale_colour_brewer(palette="Dark2")
ggsave("../figure/AB_Granulocyte_specifcA.GeneOntology.pdf", device = "pdf", width = 5, height = 2)


# MEP MK

ABcomp_binary_cellType_MEP = ABcomp_binary[ABcomp_binary[,10] > 0 & ABcomp_binary[,11] > 0 & ABcomp_binary[,16] <0 & ABcomp_binary[,17]<0,]
ABcomp_binary_cellType_MK = ABcomp_binary[ABcomp_binary[,10] < 0 & ABcomp_binary[,11] < 0 & ABcomp_binary[,16] >0 & ABcomp_binary[,17]>0,]

pheatmap(ABcomp_500k[rownames(ABcomp_binary_cellType_MK),c(10,11,16,17)], scale="none", border_color=NA, cluster_cols = F,cluster_rows = T,
         clustering_method = "ward.D", show_rownames = F, color = colorRampPalette(c("purple", "black", "yellow"))(20))
#write.table(rownames(ABcomp_binary_cellType_MK), "../data/MEP_MK.AB.bed", row.names = F, quote = F)


ABcomp_binary_cellType_MK = read.table("../data/MEP_MK.AB_expr.bed")
ABcomp_binary_cellType_MK_gene = ABcomp_binary_cellType_MK[ABcomp_binary_cellType_MK$V7 == "protein_coding", 8]


# expr
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(2,3,5:7,9,12,14,16:20,25,27,28,31,32)]
expr = expr[rowMeans(expr)>=1,]

AB_switch_expr = read.table("../data/MEP_MK.AB_expr.bed", stringsAsFactors = F, sep="\t")
AB_switch_expr = AB_switch_expr[AB_switch_expr$V7 == "protein_coding",]
AB_switch_expr$name = paste(AB_switch_expr$V1, AB_switch_expr$V2, sep = "_")
AB_switch_expr = na.omit(expr[unique(AB_switch_expr$V8),])

AB_switch_expr.m = melt(AB_switch_expr[,c("mep.1","mep.2","mk.1","mk.2")])

AB_switch_expr.m$cell = toupper(gsub(".2","", gsub(".1", "", AB_switch_expr.m$variable)))
ggplot(AB_switch_expr.m, aes(variable , log(value+1), fill=cell))+
  #stat_boxplot(geom ='errorbar', width = 0.6, coef=1.5) +
  geom_boxplot(outlier.size=-1)+
  theme_classic(base_size = 20) +xlab("")+ylab("log(TPM+1)")+
  scale_fill_manual(values = pal_d3("category10")(10)[-8])+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))



