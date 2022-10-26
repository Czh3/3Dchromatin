# loops are called by JUICER, and merged by Homer


# loop length
Loops = read.table("../analysis/loops/merged.loops.loop.scores.txt", header = T, comment.char="!", sep="\t")
Loops_length = Loops$start2 - Loops$start1
densityplot(log10(Loops_length))
plot(density(log10(Loops_length)))
# loop enrich score

loopEnrich = function(loop, hic){
  chrom = gsub("chr","",loop$chr1)
  loop1s = loop$start1 -200000
  loop1e = loop$end1 +200000
  loop2s = loop$start2 -200000
  loop2e = loop$end2 +200000
  pos1 = paste(chrom,loop1s,loop1e, sep=":")
  pos2 = paste(chrom,loop2s,loop2e, sep=":")

  
  cmd = paste0("./straw KR ",hic," ",pos1," ",pos2 ," BP 10000 > temp.straw ")
  #cmd = paste0("/lustre/user/liclab/zhangc/tools/juicer-master/CPU/common/juicer_tools dump oe KR ../Juicer/merge/ST.hic ",pos," ",pos ,"  BP 5000 > temp.straw")
  system(cmd)
  
  mat = tryCatch(read.table("temp.straw"), error=function(e){return(NA)}) 

  mat = mat[mat$V2+mat$V1 >= loop$start1*2 & mat$V2+mat$V1 <= loop$end2*2, ]
  if(nrow(mat) < 50){
    return(NA)
  }
  mat[is.na(mat)]<- 0 
  extend = 20000
  loop$start1 = loop$start1 - extend
  loop$end1 = loop$end1 + extend
  loop$start2 = loop$start2 - extend
  loop$end2 = loop$end2 - extend
  mat_l = mean(mat[mat$V1 >= loop$start1 & mat$V1 <= loop$end1 & mat$V2 >= loop$start2 & mat$V2 <= loop$end2, 3])
  mat_b = mean(mat[!(mat$V1 >= loop$start1 & mat$V1 <= loop$end1 & mat$V2 >= loop$start2 & mat$V2 <= loop$end2), 3])
  
  if(mat_b==0){
    return(0)
  }else{
    return(mat_l/mat_b)
  }
}

Loops = read.table("../analysis/loops/merged.loops.loop.scores.txt", header = T, comment.char="!", sep="\t")
rownames(Loops) = Loops[,1]

loopsER = as.data.frame(matrix(0, nrow=nrow(Loops),ncol=0))
rownames(loopsER) = rownames(Loops)
for(h in c("ST","CLP","MEP","GMP","MKP","G")){
  hic = paste0("../Juicer/merge/",h,".hic")
  loopsE = c()
  for(i in 1:nrow(Loops)){
    loopsE = c(loopsE, loopEnrich(Loops[i,], hic))
  }
  loopsER[[h]] = loopsE
}

loopsER = na.omit(loopsER)
boxplot(loopsER)

loopsER1 = (loopsER)

# normalize by colsum
col_Sum = apply(loopsER1, 2, sum)/nrow(loopsER1)
col_Sum = col_Sum/col_Sum[1]
loopsER1 = as.data.frame(t(apply(loopsER1, 1, function(x) x/col_Sum)))
#boxplot((loopsER1+1))
#myc
loopsER1["chr15:61975000-63750000",]

loopsER_sig = na.omit(loopsER1[apply(loopsER1, 1, function(x) max(x) >= 1.5) & apply(loopsER1, 1, function(x) min(x) < 1),])

breaksList = seq(0.5, 3, by = 0.05)
ph = pheatmap((loopsER_sig), scale = "none", show_rownames = F, clustering_method="ward.D2",cluster_cols = F, border_color = NA,
              annotation_row = clust3, annotation_colors = ann_colors,annotation_legend =T, #cutree_rows = 7,
              col=colorRampPalette(c("gray","white","red"))(length(breaksList)),breaks = breaksList)

breaksList = seq(-2, 2, by = 0.1)
ph = pheatmap((loopsER_sig), scale = "row", show_rownames = F, clustering_method="ward.D2",cluster_cols = F, border_color = NA, cutree_rows = 4,
         col=colorRampPalette(c("gray","white","white","red"))(length(breaksList)),breaks = breaksList)
clust3 = cutree(ph$tree_row, k=3)[ph$tree_row[["order"]]]
clust3 = as.data.frame(clust3)
clust3$clust3 = paste("cluster",clust3$clust3)
a = brewer.pal(3, "Dark2")
names(a) = paste("cluster",1:3)
ann_colors = list(
  clust3 = a
)
breaksList = seq(-2, 2, by = 0.05)
ph = pheatmap((loopsER_sig), scale = "row", show_rownames = F, clustering_method="ward.D2",cluster_cols = F, border_color = NA,
              annotation_row = clust3, annotation_colors = ann_colors,annotation_legend =T, #cutree_rows = 7,
              col=colorRampPalette(c("gray","white","white","red"))(length(breaksList)),breaks = breaksList)

pdf("../figure/loops/diffLoops_306.heatmap.pdf",5,4)
ph
dev.off()
png("../figure/loops/diffLoops_306.heatmap.png",1100,801, res=300)
ph
dev.off()

Loops = read.table("../analysis/loops/merged.loops.loop.scores.txt", header = T, comment.char="!", sep="\t")
rownames(Loops) = Loops[,1]
Loops1 = cbind(Loops[rownames(clust3), 2:7] ,clust3)
colnames(Loops1) = NA
clust3_pos = rbind(Loops1[, c(1:3,7)], Loops1[, c(4:6,7)])


write.table(Loops1[,-c(3:5)], paste0("../analysis/loops/Cluster3/Cluster1-3.loops"), quote = F, row.names = F, col.names = F, sep="\t")


for (i in 1:7){
  write.table(clust7_pos[clust7_pos[,4] == paste("cluster",i), 1:3], paste0("../analysis/loops/Cluster7/Cluster",i,".loops"), quote = F, row.names = F, col.names = F, sep="\t")
}

write.table(Loops1[,-c(3:5)], paste0("../analysis/loops/Cluster7/Cluster1-7.loops"), quote = F, row.names = F, col.names = F, sep="\t")

for (i in 1:7){
  write.table(Loops1[clust7_pos[,4] == paste("cluster",i), 1:6], paste0("../analysis/loops/Cluster7/Cluster",i,".loops.bedpe"), quote = F, row.names = F, col.names = F, sep="\t")
}




# gene expr
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(2,3,5:7,9,12,14,16:18,20,25,27,31,32)]
colnames(expr) = toupper(gsub(".1|.2","",colnames(expr)))
expr = expr[rowMeans(expr)!=0,]

  
c_gene = read.table("../analysis/loops/Cluster3/Cluster1-3.loops.gene", sep="\t", stringsAsFactors = F)
c_gene = c_gene[,c(4,8)]
colnames(c_gene) = c("cluster", "gene")
expr$gene = rownames(expr)

c_geneExpr = merge(expr,c_gene,by="gene")
c_geneExpr.m = melt(c_geneExpr)
c_geneExpr.m = unique(c_geneExpr.m)
c_geneExpr.m$variable = as.character(c_geneExpr.m$variable)

ggplot(c_geneExpr.m, aes(variable, log(value+1), fill=variable)) +
  geom_boxplot(size=1)+
  facet_wrap(.~cluster)



# C1 G
c_geneExpr.m1 = c_geneExpr.m[c_geneExpr.m$cluster == "cluster 1",]
c_geneExpr.m1[c_geneExpr.m1$variable != "GR", 3] <-  "other"
ggplot(c_geneExpr.m1, aes(variable, log(value+1), fill=variable)) +
  geom_boxplot(size=1)+
  theme_classic(base_size = 20) + xlab("Cluster 1") + ylab("log(TPM+1)") +
  scale_fill_manual( values = c(brewer.pal(7, "Dark2")[1],"gray"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        legend.position = "none",
        axis.text = element_text(colour = "black"))
ggsave("../figure/loops/MEP_specific.geneExpr.pdf", device = "pdf", width = 4, height = 4)

t.test(c_geneExpr.m1[c_geneExpr.m1$variable=="other",4], c_geneExpr.m1[c_geneExpr.m1$variable!="other",4])
# C2 
c_geneExpr.m1 = c_geneExpr.m[c_geneExpr.m$cluster == "cluster 2",]
c_geneExpr.m1[c_geneExpr.m1$variable != "ST.HSC", 3] <-  "other"
ggplot(c_geneExpr.m1, aes(variable, log(value+1), fill=variable)) +
  geom_boxplot(size=1)+
  theme_classic(base_size = 20) + xlab("Cluster 2") + ylab("log(TPM+1)") +
  scale_fill_manual( values = c(brewer.pal(7, "Dark2")[2],"gray"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        legend.position = "none",
        axis.text = element_text(colour = "black"))
ggsave("../figure/loops/G_specific.geneExpr.pdf", device = "pdf", width = 4, height = 4)
t.test(c_geneExpr.m1[c_geneExpr.m1$variable=="other",4], c_geneExpr.m1[c_geneExpr.m1$variable!="other",4])


# binary
loopsER_binary = ifelse(loopsER1 >= 1.5, 1, 0)
pheatmap((loopsER_binary), scale = "none", show_rownames = F, clustering_method="ward.D2",cluster_cols = F, border_color = NA,
              annotation_row = clust7, annotation_colors = ann_colors,annotation_legend =T, #cutree_rows = 7,
              col=colorRampPalette(c("gray","red"))(2))





### homer
loops_homer = read.table("../analysis/loops/merged.loops.loop.scores.txt",sep="\t", comment.char = "!", header = T)
rownames(loops_homer) = loops_homer[,1]
loops_homer = loops_homer[,11:20]
colnames(loops_homer) = gsub("...homer_merge.|.Score", "", colnames(loops_homer) )
loops_homer = loops_homer[,-c(1,3,9)]

col_Sum = apply(loops_homer, 2, sum)/nrow(loops_homer)
col_Sum = col_Sum/col_Sum[1]
loops_homer = as.data.frame(t(apply(loops_homer, 1, function(x) x/col_Sum)))
boxplot(log(loops_homer+1))

#

#loops_homer["chr15:61975000-63750000",]


loops_homer_sig = na.omit(loops_homer[apply(loops_homer, 1, function(x) max(x)/mean(x) >= 1.5),])
breaksList = seq(-2, 2, by = 0.1)
ph = pheatmap((loops_homer_sig), scale = "row", show_rownames = F, clustering_method="ward.D",cluster_cols = F, border_color = NA,
              #annotation_row = clust7, annotation_colors = ann_colors,
              col=colorRampPalette(c("gray","white","red"))(length(breaksList)),breaks = breaksList)

clust7 = cutree(ph$tree_row, k=6)[ph$tree_row[["order"]]]
clust7 = as.data.frame(clust7)
clust7$clust7 = paste("cluster",clust7$clust7)
a = brewer.pal(6, "Dark2")
names(a) = paste("cluster",1:6)
ann_colors = list(
  clust7 = a
)
breaksList = seq(-2, 2, by = 0.1)
ph = pheatmap((loops_homer_sig), scale = "row", show_rownames = F, clustering_method="ward.D",cluster_cols = F, border_color = NA,
         annotation_row = clust7, annotation_colors = ann_colors,
         col=colorRampPalette(c("gray","white","white","red"))(length(breaksList)),breaks = breaksList)

pdf("../figure/loops/diffLoops_911.heatmap.pdf",5,4)
ph
dev.off()
png("../figure/loops/diffLoops_911.heatmap.png",1100,801, res=300)
ph
dev.off()




Loops = read.table("../analysis/loops/merged.loops.loop.scores.txt", header = T, comment.char="!", sep="\t")
rownames(Loops) = Loops[,1]
Loops1 = cbind(Loops[rownames(clust7), 2:7] ,clust7)
colnames(Loops1) = NA
clust7_pos = rbind(Loops1[, c(1:3,7)], Loops1[, c(4:6,7)])


write.table(Loops1[,-c(3:5)], paste0("../analysis/loops/Cluster6/Cluster1-6.loops"), quote = F, row.names = F, col.names = F, sep="\t")

for (i in 1:6){
  write.table(clust7_pos[clust7_pos[,4] == paste("cluster",i), 1:3], paste0("../analysis/loops/Cluster6/Cluster",i,".loops"), quote = F, row.names = F, col.names = F, sep="\t")
}


# gene expr
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(2,3,5:7,9,12,14,16:18,20,25,27,31,32)]
colnames(expr) = toupper(gsub(".1|.2","",colnames(expr)))
expr = expr[rowMeans(expr)!=0,]


c_gene = read.table("../analysis/loops/Cluster6/Cluster1-6.loops.gene", sep="\t", stringsAsFactors = F)
c_gene = c_gene[,c(4,8)]
colnames(c_gene) = c("cluster", "gene")
expr$gene = rownames(expr)

c_geneExpr = merge(expr,c_gene,by="gene")
c_geneExpr = unique(c_geneExpr)
c_geneExpr.m = melt(c_geneExpr)
c_geneExpr.m = unique(c_geneExpr.m)
c_geneExpr.m$variable = as.character(c_geneExpr.m$variable)

ggplot(c_geneExpr.m, aes(variable, log(value+1), fill=variable)) +
  geom_boxplot(size=1)+
  facet_wrap(.~cluster)

a = c_geneExpr[c_geneExpr$cluster == "cluster 2",]
a = a[a[,14] / apply(a[,-c(1,18)], 1, mean) > 5,]
a[a$GR > 20,]

# C3
c_geneExpr.m1 = c_geneExpr.m[c_geneExpr.m$cluster == "cluster 3",]
c_geneExpr.m1[c_geneExpr.m1$variable != "CLP", 3] <-  "other"
ggplot(c_geneExpr.m1, aes(variable, log(value+1), fill=variable)) +
  geom_boxplot(size=1)+
  theme_classic(base_size = 20) + xlab("Cluster 3") + ylab("log(TPM+1)") +
  scale_fill_manual( values = c(brewer.pal(6, "Dark2")[3],"gray"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        legend.position = "none",
        axis.text = element_text(colour = "black"))
ggsave("../figure/loops/CLP_specific.geneExpr.pdf", device = "pdf", width = 4, height = 4)

#t.test(log(c_geneExpr.m1[c_geneExpr.m1$variable=="other",4]+1), log(c_geneExpr.m1[c_geneExpr.m1$variable!="other",4]+1))
wilcox.test(log(c_geneExpr.m1[c_geneExpr.m1$variable=="other",4]+1), log(c_geneExpr.m1[c_geneExpr.m1$variable!="other",4]+1))

# C2 
c_geneExpr.m1 = c_geneExpr.m[c_geneExpr.m$cluster == "cluster 2",]
c_geneExpr.m1[c_geneExpr.m1$variable != "GR", 3] <-  "other"
ggplot(c_geneExpr.m1, aes(variable, log(value+1), fill=variable)) +
  geom_boxplot(size=1)+
  theme_classic(base_size = 20) + xlab("Cluster 2") + ylab("log(TPM+1)") +
  scale_fill_manual( values = c(brewer.pal(6, "Dark2")[2],"gray"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        legend.position = "none",
        axis.text = element_text(colour = "black"))
ggsave("../figure/loops/G_specific.geneExpr.pdf", device = "pdf", width = 4, height = 4)
#t.test(log(c_geneExpr.m1[c_geneExpr.m1$variable=="other",4]+1), log(c_geneExpr.m1[c_geneExpr.m1$variable!="other",4]+1))
wilcox.test(log(c_geneExpr.m1[c_geneExpr.m1$variable=="other",4]+1), log(c_geneExpr.m1[c_geneExpr.m1$variable!="other",4]+1))


# C4
c_geneExpr.m1 = c_geneExpr.m[c_geneExpr.m$cluster == "cluster 4",]
c_geneExpr.m1[c_geneExpr.m1$variable != "MEP", 3] <-  "other"
ggplot(c_geneExpr.m1, aes(variable, log(value+1), fill=variable)) +
  geom_boxplot(size=1)+
  theme_classic(base_size = 20) + xlab("Cluster 4") + ylab("log(TPM+1)") +
  scale_fill_manual( values = c(brewer.pal(6, "Dark2")[4],"gray"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        legend.position = "none",
        axis.text = element_text(colour = "black"))
ggsave("../figure/loops/MEP_specific.geneExpr.pdf", device = "pdf", width = 4, height = 4)
#t.test(log(c_geneExpr.m1[c_geneExpr.m1$variable=="other",4]+1), log(c_geneExpr.m1[c_geneExpr.m1$variable!="other",4]+1))
wilcox.test(log(c_geneExpr.m1[c_geneExpr.m1$variable=="other",4]+1), log(c_geneExpr.m1[c_geneExpr.m1$variable!="other",4]+1))








# JUICER
Loops = read.table("/lustre/user/liclab/zhangc/proj/blood3D_3/analysis/loops/merge.Juicer.FDR.bed", sep="\t", header = F)
Loops = na.omit(Loops)
rownames(Loops) = gsub(" ", "", apply(Loops[,1:6], 1, function(x) paste(x, collapse = ":")))
Loops = Loops[,-c(1:6)]
colnames(Loops) = c("MEP","LT","CMP","MPP","CLP","MKP","G","MK","GMP","ST")
Loops = Loops[,c("ST","CLP","CMP","MEP","GMP","MKP","MK","G")]
Loops[Loops==0] <- 1e-50
Loops[Loops>1] <- 1
Loops = -log10(Loops)


Loops_sig = Loops[apply(Loops, 1, max) > 5 & apply(Loops, 1, mean) < 2, ]
breaksList = seq(0, 10, by = 0.1)
ph = pheatmap((Loops_sig), scale = "none", show_rownames = F, clustering_method="ward.D",cluster_cols = F, border_color = NA,cutree_rows = 5,
         col=colorRampPalette(c("gray","white","yellow","orange","red"))(length(breaksList)),breaks = breaksList)

clust5 = cutree(ph$tree_row, k=5)[ph$tree_row[["order"]]]
clust5 = as.data.frame(clust5)
clust5$clust5 = paste("cluster",clust5$clust5)
a = brewer.pal(5, "Dark2")
names(a) = paste("cluster",1:5)
ann_colors = list(
  clust5 = a
)
breaksList = seq(0, 8, by = 0.1)
pheatmap((Loops_sig), scale = "none", show_rownames = F, clustering_method="ward.D",cluster_cols = F, border_color = NA,
         annotation_row = clust5, annotation_colors = ann_colors,
         col=colorRampPalette(c("white","white","red"))(length(breaksList)),breaks = breaksList)



write.table(clust5, paste0("../analysis/loops/Cluster5/Cluster1-5.loops"), quote = F, row.names = T, col.names = F, sep="\t")


# gene expr
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(2,3,5:7,9,12,14,16:18,20,25,27,31,32)]
colnames(expr) = toupper(gsub(".1|.2","",colnames(expr)))
expr = expr[rowMeans(expr)!=0,]


c_gene = read.table("../analysis/loops/Cluster5/Cluster1-5.loops.gene", sep="\t", stringsAsFactors = F)
c_gene = c_gene[,c(4,8)]
colnames(c_gene) = c("cluster", "gene")
expr$gene = rownames(expr)

c_geneExpr = merge(expr,c_gene,by="gene")
c_geneExpr.m = melt(c_geneExpr)
c_geneExpr.m = unique(c_geneExpr.m)
c_geneExpr.m$variable = as.character(c_geneExpr.m$variable)

ggplot(c_geneExpr.m, aes(variable, log(value+1), fill=variable)) +
  geom_boxplot(size=1)+
  facet_wrap(.~cluster)




# gene at loops anchors
# gene expr
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(2,3,5:7,9,12,14,16:18,20,25,27,31,32)]
colnames(expr) = toupper(colnames(expr))
expr = expr[rowMeans(expr)!=0,]

anchor_genes = read.table("../analysis/loops/merge.loopsAnchor.geneTSS.bed")
anchor_genes_expr = expr[unique(anchor_genes$V7),]

breaksList = seq(-3, 3, by = 0.1)
pheatmap(anchor_genes_expr, scale="row",show_rownames = F, clustering_method="ward.D",
         col=colorRampPalette(c("blue","white","red"))(length(breaksList)),breaks = breaksList)




# motif
motif = read.table("../analysis/loops/merge.loopsAnchor.Motif/knownResults.txt", sep = "\t", comment.char = "!",header=T, stringsAsFactors = F)
motif = motif[1:5,]
motif$Motif = sapply(strsplit(motif$Motif.Name, "/"), "[[", 1)
motif$Motif = factor(motif$Motif, levels = rev(motif$Motif))
ggplot(motif, aes(Motif, -log10(P.value))) +
  geom_bar(stat="identity")+
  coord_flip() +
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"))

ggsave("../figure/loops/merge.motif.pdf", device = "pdf", width = 4, height = 4)


