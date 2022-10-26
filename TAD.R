### TAD

# number
TAD_num = read.table("../analysis/TAD_boundary/TAD_number.txt")
TAD_num$cell = stringr::str_sub(TAD_num$V2, 1, -2)
TAD_num$rep = stringr::str_sub(TAD_num$V2, -1, -1)
TAD_num$cell = factor(TAD_num$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))

ggplot(TAD_num, aes(cell, V1, col=cell))+
  geom_point(aes(shape=rep), size=4) +
  theme_classic(base_size = 20) + xlab("") + ylab("TAD numbers") +
  scale_color_brewer(palette="Paired") +
  theme(legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black", angle = 30, hjust = 1,vjust = 1))
ggsave("../figure/TAD/TAD_number.pdf", device = "pdf", width = 6, height = 4)


# structrue gene
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)

expr = expr[,c(2,3,5:7,9,12,14,16:20,25,27,28,31,32)]
Structure_gene = c("Ctcf","Smc3", "Rad21","Yy1" )
pheatmap(log(expr[Structure_gene,]+1), scale="row", border_color = NA,
         colorRampPalette(rev(brewer.pal(n = 11, name = 'RdBu')))(20))
expr.m = melt(as.matrix(expr[Structure_gene,]))
expr.m$cell = toupper(gsub(".2","", gsub(".1", "", expr.m$Var2 )))
expr.m$cell = factor(expr.m$cell, level = c("LT.HSC", "ST.HSC", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","GR"))
expr.m$rep = gsub("[a-z|.]", "", expr.m$Var2 )
ggplot(expr.m, aes(cell, value, rep = rep, fill=cell)) +
  geom_bar(stat="identity",position=position_dodge(), color="black") + ylab("TPM") +
  facet_grid(Var1~., scales = "free") +
  theme_classic(base_size=20)+
  theme(axis.text.x = element_text(size = 13,angle = 45, vjust = 1, hjust = 1))+
  scale_fill_brewer(palette="Paired")
  

# GMP G
TAD_GMP_G = read.table("../analysis/TAD_loop/Diff_GMP_G.tad.scores.txt")
TAD_GMP_G$ratio = TAD_GMP_G$V11 / TAD_GMP_G$V12
TAD_GMP_G = TAD_GMP_G[order(TAD_GMP_G$ratio, decreasing = T),]

TAD_GMP_G = TAD_GMP_G[TAD_GMP_G$V11 > 1.6 & TAD_GMP_G$V12 < 1.1,]
write.table(TAD_GMP_G[,-1], "../analysis/TAD_loop/Diff_GMP_G.tad.bed", row.names = F, col.names = F, quote = F, sep="\t")
TAD_GMP_G = read.table("../analysis/TAD_loop/Diff_GMP_G.tad.anno.bed")
TAD_GMP_G = TAD_GMP_G[TAD_GMP_G$V16 == "protein_coding",]

options(scipen=999)
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)

#expr = expr[,c(2,3,5:7,9,12,14,16:20,25,27,28,31,32)]


expr = expr[,c(18,31,25,27)]

TAD_GMP_G_expr = expr[as.character(unique(TAD_GMP_G$V17)),]
TAD_GMP_G_expr = TAD_GMP_G_expr[rowMeans(TAD_GMP_G_expr) >= 1 ,]
pheatmap(TAD_GMP_G_expr, scale = "row")
boxplot(log10(TAD_GMP_G_expr+1))
boxplot(log10(expr[rowMeans(expr)>=1,]+1))


# tad score
tads = read.table("../analysis/TAD_loop/merge.tad.2D.tad.scores.txt")
rownames(tads) = tads$V1
tads = tads[,11:30]
#colnames(tads) = c("LT","ST","MPP","CMP","CLP","MEP","GMP","MKP","MK","G")
colnames(tads) = c("CLP1",   "G1",     "LT1",    "MK1",    "MPP1",  "CLP2",   "G2",     "LT2",    "MK2",    "MPP2",  "CMP1",   "GMP1",   "MEP1",   "MKP1",   "ST1",  "CMP2",   "GMP2",   "MEP2",   "MKP2",   "ST2")
tads = tads[,-5] # remove MPP1
boxplot(log(tads+1))
tads = tads[rowMeans(tads)!=0,]

breaksList = seq(-2, 2, by = 0.1)
png("../figure/TADscore_Homer.png",1000,1000, res=300)
pheatmap(tads, scale = "row", show_rownames = F, cluster_cols = T,  clustering_method="ward.D",border_color = NA,
         col=colorRampPalette(rev(brewer.pal(n = 11, name = 'RdYlBu')))(length(breaksList)),breaks = breaksList)
dev.off()

tads = tads[,c("LT1","LT2","ST1","ST2","MPP2","CLP1","CLP2",
                               "CMP1","CMP2","MEP1","MEP2","GMP1","GMP2","MKP1","MKP2","MK1","MK2","G1","G2")]


group = factor(c("LT","LT","ST","ST","MPP","CLP","CLP",
                 "CMP","CMP","MEP","MEP","GMP","GM2","MKP","MKP","MK","MK","G","G"))

results =sapply(1:nrow(tads) ,function(x){ summary(aov(as.numeric(tads[x,])~group))[[1]][["Pr(>F)"]][1]} )
tads_sig = tads[results<0.05, ]
tads_sig = tads_sig[apply(tads_sig,1,max) > 2 & apply(tads_sig,1,min) < 1.2, ]
tads_sig[tads_sig>2.5] <- 2.5
pheatmap(tads_sig, scale="none",color = colorRampPalette(c("blue", "white", "red"))(20),
         clustering_method = "ward.D", border_color = NA, show_rownames = F)



write.table(as.data.frame(gene_clust), "../analysis/TAD_loop/merge.tad.2D.tad.scores.bed", quote = F)

expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(2,3,5:7,9,12,14,16:20,25,27,28,31,32)]

tads_expr = read.table("../analysis/TAD_loop/merge.tad.2D.tad.scores.anno.bed",stringsAsFactors = F)
expr = expr[rowMeans(expr)>=1,]
tads_expr_1 = na.omit(expr[unique(tads_expr[tads_expr$V4 == 2,8]),])
tads_expr_1[tads_expr_1>10] <- 10
pheatmap(tads_expr_1, scale = "none")
boxplot(log(na.omit(expr[unique(tads_expr[tads_expr$V9 == 1,5]),]+1)))



