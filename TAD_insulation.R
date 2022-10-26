####
# insulation score at boundaries
#


setwd("/lustre/user/liclab/zhangc/proj/blood3D_3/script/")
# LT
IS_LT1 = read.table("../data/InsulationScore_at_LT1_boundary.tab",sep="\t",skip=2,row.names = 1,stringsAsFactors = F)
colnames(IS_LT1) = c("gene",1:(ncol(IS_LT1)-1))
IS_LT1 = IS_LT1[,-1]
IS_LT1 = IS_LT1[rownames(IS_LT1) != "MPP1",]

IS_LT1 = melt(as.matrix(IS_LT1))

IS_LT1$cell = stringr::str_sub(IS_LT1$Var1, 1, -2)
IS_LT1$rep = stringr::str_sub(IS_LT1$Var1, -1, -1)
IS_LT1$cell = factor(IS_LT1$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))


ggplot(IS_LT1, aes(Var2, value, group= Var1, color=cell))+
  #geom_smooth(method="loess", se=F, fullrange=T) +
  geom_line(size=1) +
  theme_classic(base_size=20)  + xlab("") + ylab("")+
  scale_color_brewer(palette="Paired") +
  scale_x_continuous(breaks=c(0,16,30),
                     labels=c("-600K", "boundaries", "+600K"))
ggsave("../figure/InsulationScore_at_LT1_boundary.pdf", device = "pdf", width = 6, height = 4)


# G
IS_G1 = read.table("../data/InsulationScore_at_G1_boundary.tab",sep="\t",skip=2,row.names = 1,stringsAsFactors = F)
colnames(IS_G1) = c("gene",1:(ncol(IS_G1)-1))
IS_G1 = IS_G1[,-1]
IS_G1 = IS_G1[rownames(IS_G1) != "MPP1",]

IS_G1 = melt(as.matrix(IS_G1))

IS_G1$cell = stringr::str_sub(IS_G1$Var1, 1, -2)
IS_G1$rep = stringr::str_sub(IS_G1$Var1, -1, -1)
IS_G1$cell = factor(IS_G1$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))


ggplot(IS_G1, aes(Var2, value, group= Var1, color=cell))+
  #geom_smooth(method="loess", se=F, fullrange=T) +
  geom_line() +
  theme_classic(base_size=20)  + xlab("") + ylab("")+
  scale_color_brewer(palette="Paired") +
  scale_x_continuous(breaks=c(0,16,30),
                     labels=c("-600K", "boundaries", "+600K"))
ggsave("../figure/InsulationScore_at_G1_boundary.pdf", device = "pdf", width = 6, height = 4)


#
IS_LT1 = read.table(gzfile("../data/InsulationScore_at_LT1_boundary.gz"),sep="\t",skip=1,stringsAsFactors = F)
rownames(IS_LT1) = IS_LT1$V4
IS_LT1 = IS_LT1[,-c(1:6)]
IS_LT1 = as.matrix(IS_LT1)
# LT1
pheatmap(na.omit(IS_LT1[,271:300]), cluster_rows = T, cluster_cols = F, show_colnames = F,show_rownames = F, color = colorRampPalette(c("yellow", "black", "purple"))(30), scale="row")
# G1
pheatmap(na.omit(IS_LT1[,151:180]), cluster_rows = T, cluster_cols = F, show_colnames = F,show_rownames = F,  color = colorRampPalette(c("yellow", "black", "purple"))(30), scale="row")
# MK1
pheatmap(na.omit(IS_LT1[,391:420]), cluster_rows = T, cluster_cols = F, show_colnames = F,show_rownames = F,  color = colorRampPalette(c("yellow", "black", "purple"))(30), scale="row")
# CLP1
pheatmap(na.omit(IS_LT1[,1:30]), cluster_rows = T, cluster_cols = F, show_colnames = F,show_rownames = F,  color = colorRampPalette(c("yellow", "black", "purple"))(30), scale="row")


# G1 LT1
LT1_G1 = rbind( IS_LT1[,271:300],IS_LT1[,151:180])
LT1_G1 = na.omit(LT1_G1)
#LT1_G1.1 = t(apply(LT1_G1, 1, function(x){x/max(x)}))

pheatmap(LT1_G1, cluster_rows = F, cluster_cols = F, show_rownames = F, color = rainbow(30)[1:20], scale="row")



###  TAD broundaries loss 
# LT1
TAD_boundary = read.table("../analysis/TAD_boundary/LT1.merge.bed")
TAD_boundary = TAD_boundary[TAD_boundary$V1 != "chrY",]
TAD_boundary = TAD_boundary[TAD_boundary$V1 != "chrX",]
for(i in 1:nrow(TAD_boundary)){
  position_b = (floor((TAD_boundary[i, 3] - TAD_boundary[i, 2])/40000 /2))*40000 + TAD_boundary[i, 2]
  rownames(TAD_boundary)[i] = paste(TAD_boundary[i, 1], position_b, sep="-")
  TAD_boundary[i,2] = position_b - 200000
  TAD_boundary[i,3] = position_b + 200000
}
for(i in  Sys.glob("../analysis/TAD/*/allChr.insulation.bedGraph") ){
  IS_G = read.table(i, sep="\t")
  rownames(IS_G) = paste(IS_G$V1, IS_G$V2, sep="-")
  IS_ratio = c()
  for(j in 1:nrow(TAD_boundary)){
    b = seq(TAD_boundary[j,2],TAD_boundary[j,3], by=40000)
    b = paste(TAD_boundary[j,1], b, sep="-")
    #ratio = min(IS_G[b,4]) -  mean(IS_G[b,4])
    ratio = min(IS_G[b,4])
    IS_ratio = c(IS_ratio, ratio)
  }
  TAD_boundary = cbind(TAD_boundary, IS_ratio)
  
  #TAD_boundary = cbind(TAD_boundary, IS_G[rownames(TAD_boundary),4])
}
TAD_boundary = TAD_boundary[,-c(1:3)]
colnames(TAD_boundary) = sapply(strsplit(Sys.glob("../analysis/TAD/*/allChr.insulation.bedGraph") , "/"), "[", 4)
TAD_boundary = na.omit(TAD_boundary)
TAD_boundary = TAD_boundary[,c("LT1","LT2","ST1","ST2","MPP2","CLP1","CLP2",
                 "CMP1","CMP2","MEP1","MEP2","GMP1","GMP2","MKP1","MKP2","MK1","MK2","G1","G2")]

TAD_boundary_sig = sort(apply(TAD_boundary, 1 , sd))

pdf("../figure/TAD/InsulationScore_StandardDeviation.pdf",4,3.5)
plot(density(TAD_boundary_sig), lwd=1, main="Standard Deviation of Insulation Score")
abline(v = quantile(TAD_boundary_sig, prob=.95), col="red", lwd=1.5, lty=2)
dev.off()


TAD_boundary_sig1 = TAD_boundary_sig[TAD_boundary_sig >= quantile(TAD_boundary_sig, prob=.95)]
TAD_boundary_sig1 = TAD_boundary[apply(TAD_boundary,1,min) < -0.6 & apply(TAD_boundary,1,max)  > 0,]

pdf("../figure/TAD/InsulationScore_heatmap.pdf",4,5)
pheatmap(TAD_boundary[names(TAD_boundary_sig1),],scale="row",color = colorRampPalette(c("gray", "white", "red"))(20), border_color = NA, show_rownames = F)
dev.off()

pdf("../figure/TAD/LossTAD_InsulationScore_heatmap.pdf",4,4)
TAD_boundary_sig1[TAD_boundary_sig1< -2] <- -2
pheatmap(TAD_boundary_sig1,scale="none",color = colorRampPalette(c("black","gray", "white", "red"))(20),
         clustering_method = "ward.D", border_color = NA, show_rownames = F)
dev.off()


# anova
group = factor(c("LT","LT","ST","ST","MPP","CLP","CLP",
                 "CMP","CMP","MEP","MEP","GMP","GM2","MKP","MKP","MK","MK","G","G"))

results =sapply(1:nrow(TAD_boundary) ,function(x){ summary(aov(as.numeric(TAD_boundary[x,])~group))[[1]][["Pr(>F)"]][1]} )
results = p.adjust(results, method="bonferroni")
TAD_boundary_sig = TAD_boundary[results<0.05, ]
TAD_boundary_sig[TAD_boundary_sig< -1.5] <- -1.5
pheatmap(TAD_boundary_sig,scale="none",color = colorRampPalette(c("blue", "white", "red"))(20),
         clustering_method = "ward.D", border_color = NA, show_rownames = F)



# plot hic matrix
library(reshape2)
library(gridGraphics)
plotMatrix = function(hic, pos){
  cmd = paste0("./straw KR ",hic ," ",pos," ",pos ," BP 50000 > temp.straw ")
  system(cmd)
  
  mat = tryCatch(read.table("temp.straw"), error=function(e){return(matrix(0, nrow = 30, ncol = 30))}) 
  mat = mat[mat$V2 - mat$V1 < 3000000, ]
  mat = dcast(mat, V1 ~ V2, value.var="V3")
  mat[is.na(mat)] <- 0
  rownames(mat) = mat$V1
  mat = mat[,-1]
  mat = as.matrix(mat)
  nx = ncol(mat)
  ny = nrow(mat)
  mat1 = matrix(0, nrow = 60, ncol=nx+2)
  for(i in 1:60){
    mat1[i,] = c(rep(0,floor(i/2)+1), as.numeric(diag(mat[,i:nx])), rep(0,floor((i+1)/2)))
  }

  #
  q = quantile(mat, probs=0.97)

  mat1[mat1>q] <- q
  
  par(mar=c(0,0,0,0))
  image((t(mat1)), col = colorRampPalette(c( "white", "red"))(n = 20),frame.plot=F,axes=F)

}

png(paste0("../figure/TAD/heatmap_chr5_38m-50m.png"),res=300, width = 500, height = 700)
par(mfrow=c(9,1)) 
for (i in c("LT", "ST", "CLP","CMP","GMP","MEP","MKP","MK","G")){
  pos="5:38000000:50000000"
  mat_path = paste0("../Juicer/merge/",i,".hic")
  par(mar=c(0,0,0,0))
  plotMatrix(mat_path, pos)
}
dev.off()
