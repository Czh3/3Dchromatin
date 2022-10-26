#####
# population-based 3D modeling of hic data
#


library(reshape2)
library(plotly)
library(igraph)

mat_path = "../hicpro/G1/iced/500000/G1_500000_iced.matrix"
chrom = "chr1"
abs_bed = "../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed"

hic_model <- function(mat_path, abs_bed, chrom){
  
abs_bed = read.table(abs_bed)
chr1_abs = abs_bed[abs_bed$V1 == chrom, 4]

mat = read.table(mat_path, stringsAsFactors = F)
mat_chr1 = mat[mat$V1 %in% chr1_abs & mat$V2 %in% chr1_abs,]
mat_chr1 = dcast(mat_chr1, V1 ~ V2, value.var="V3")
rownames(mat_chr1) = mat_chr1$V1
mat_chr1 = mat_chr1[,-1]

# symmetric matrix
bin_select = intersect(colnames(mat_chr1), rownames(mat_chr1))
mat_chr1 = mat_chr1[bin_select, bin_select]
mat_chr1[lower.tri(mat_chr1)] = t(mat_chr1)[lower.tri(mat_chr1)]
mat_chr1 = as.matrix(mat_chr1)

mat_chr1[is.na(mat_chr1)] <- 0

# distance matrix
mat_min = min(mat_chr1[mat_chr1!=0]) * 0.1

mat_chr1[mat_chr1 == 0] <- mat_min
C = 1 / mat_chr1
#C = apply(mat_chr1, 1, function(x) ifelse(x==0, mat_max, x** -1))

# shortest path

D = graph.adjacency(C, mode="undirected", weighted = T, diag=F)
D = shortest.paths(D, algorithm = "dijkstra")

#image(log(mat_chr1))
#image(log(C))
#image(D)


fit = cmdscale(D, eig=T, k=3)
x = fit$points[,1]
y = fit$points[,2]
z = fit$points[,3]
pos = as.data.frame(fit$points)


pos$cor = 1:nrow(pos)

ax <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)
scene = list(
  xaxis = ax,
  yaxis = ax,
  zaxis = ax)

plot_ly(pos, x = ~V1, y = ~V2, z = ~V3, color = ~cor)  %>% 
  layout(scene=scene)

return(pos)

}

# 500k
model_LT =  hic_model("../hicpro/LT1/iced/500000/LT1_500000_iced.matrix","../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", "chr1")
hic_model("../hicpro/MEP1/iced/500000/MEP1_500000_iced.matrix","../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", "chr1")
#hic_model("../hicpro/GMP1/iced/500000/GMP1_500000_iced.matrix","../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", "chr1")
#hic_model("../hicpro/MK1/iced/500000/MK1_500000_iced.matrix","../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", "chr1")
model_G = hic_model("../hicpro/G1/iced/500000/G1_500000_iced.matrix","../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", "chr1")

pc1_G = read.table("../homer/G1.PC1.bedGraph", skip = 1)
pc1_G = pc1_G[pc1_G$V1 == "chr1",]
model_G$cor = c(-1,pc1_G$V4)
model_G[model_G$cor > 1, 4] <- 1
model_G[model_G$cor < -1, 4] <- -1
plot_ly(model_G, x = ~V1, y = ~V2, z = ~V3, color = ~cor, colors = c("blue","red"))  %>% 
  layout(scene=scene)


set.seed(100)
model_G$cor = sample(c(-1,pc1_G$V4))
model_G[model_G$cor > 1, 4] <- 1
model_G[model_G$cor < -1, 4] <- -1
plot_ly(model_G, x = ~V1, y = ~V2, z = ~V3, color = ~cor, colors = c("blue","red"))  %>% 
  layout(scene=scene)

H3k4_G = read.table("../data/H3K4me1_G_chr1_500K.bed")
model_G$cor = c(2541.5,H3k4_G$V5)
model_G$cor = log10(model_G$cor+1)
plot_ly(model_G, x = ~V1, y = ~V2, z = ~V3, color = ~cor, colors = c("black","yellow"))  %>% 
  layout(scene=scene)


pc1_LT = read.table("../homer/LT1.PC1.bedGraph", skip = 1)
pc1_LT = pc1_LT[pc1_LT$V1 == "chr1",]
model_LT$cor = c(-1,pc1_LT$V4)
model_LT[model_LT$cor > 1, 4] <- 1
model_LT[model_LT$cor < -1, 4] <- -1
plot_ly(model_LT, x = ~V1, y = ~V2, z = ~V3, color = ~cor, colors = c("blue","red"))  %>% 
  layout(scene=scene)



#100k
#
model_LT =  hic_model("../hicpro/LT1/iced/100000/LT1_100000_iced.matrix","../hicpro/CLP1/raw/100000/CLP1_100000_abs.bed", "chr1")
model_LT$pos = floor(as.numeric(rownames(model_LT))/5)*5
p=plot_ly(model_LT, x = ~V1, y = ~V2, z = ~V3, color = ~pos)  %>% 
  layout(scene=scene)


pc1_G = read.table("../homer/LT1.PC1.bedGraph", skip = 1)
pc1_G = pc1_G[pc1_G$V1 == "chr1",]
model_LT$pos = floor(as.numeric(rownames(model_LT))/5)*5
model_LT$cor = pc1_G[match(model_LT$pos,pc1_G$V2/100000), 4]
model_LT[model_LT$cor > 1, 4] <- 1
model_LT[model_LT$cor < -1, 4] <- -1
plot_ly(model_LT, x = ~V1, y = ~V2, z = ~V3, color = ~cor, colors = c("blue","red"))  %>% 
  layout(scene=scene)

#
model_G =  hic_model("../hicpro/G1/iced/100000/G1_100000_iced.matrix","../hicpro/CLP1/raw/100000/CLP1_100000_abs.bed", "chr1")
model_G$pos = floor(as.numeric(rownames(model_G))/5)*5
plot_ly(model_G, x = ~V1, y = ~V2, z = ~V3, color = ~pos)  %>% 
  layout(scene=scene)

pc1_G = read.table("../homer/G1.PC1.bedGraph", skip = 1)
pc1_G = pc1_G[pc1_G$V1 == "chr1",]

model_G$cor = pc1_G[match(model_G$pos,pc1_G$V2/100000), 4]
model_G = na.omit(model_G)
model_G[model_G$cor > 1, 4] <- 1
model_G[model_G$cor < -1, 4] <- -1
plot_ly(model_G, x = ~V1, y = ~V2, z = ~V3, color = ~cor, colors = c("blue","red"))  %>% 
  layout(scene=scene)



#
model_G =  hic_model("../hicpro/MEP1/iced/100000/MEP1_100000_iced.matrix","../hicpro/CLP1/raw/100000/CLP1_100000_abs.bed", "chr1")
model_G$pos = floor(as.numeric(rownames(model_G))/5)*5
plot_ly(model_G, x = ~V1, y = ~V2, z = ~V3, color = ~pos)  %>% 
  layout(scene=scene)

pc1_G = read.table("../homer/MEP1.PC1.bedGraph", skip = 1)
pc1_G = pc1_G[pc1_G$V1 == "chr1",]

model_G$cor = pc1_G[match(model_G$pos,pc1_G$V2/100000), 4]
model_G = na.omit(model_G)
model_G[model_G$cor > 1, 4] <- 1
model_G[model_G$cor < -1, 4] <- -1
plot_ly(model_G, x = ~V1, y = ~V2, z = ~V3, color = ~cor, colors = c("blue","red"))  %>% 
  layout(scene=scene)


#
model_G =  hic_model("../hicpro/CLP1/iced/100000/CLP1_100000_iced.matrix","../hicpro/CLP1/raw/100000/CLP1_100000_abs.bed", "chr1")
model_G$pos = floor(as.numeric(rownames(model_G))/5)*5
plot_ly(model_G, x = ~V1, y = ~V2, z = ~V3, color = ~pos)  %>% 
  layout(scene=scene)

pc1_G = read.table("../homer/CLP1.PC1.bedGraph", skip = 1)
pc1_G = pc1_G[pc1_G$V1 == "chr1",]

model_G$cor = pc1_G[match(model_G$pos,pc1_G$V2/100000), 4]
model_G = na.omit(model_G)
model_G[model_G$cor > 1, 4] <- 1
model_G[model_G$cor < -1, 4] <- -1
plot_ly(model_G, x = ~V1, y = ~V2, z = ~V3, color = ~cor, colors = c("blue","red"))  %>% 
  layout(scene=scene)

#
model_G =  hic_model("../hicpro/CMP1/iced/100000/CMP1_100000_iced.matrix","../hicpro/CLP1/raw/100000/CLP1_100000_abs.bed", "chr1")
model_G$pos = floor(as.numeric(rownames(model_G))/5)*5
plot_ly(model_G, x = ~V1, y = ~V2, z = ~-V3, color = ~pos)  %>% 
  layout(scene=scene)

pc1_G = read.table("../homer/CMP1.PC1.bedGraph", skip = 1)
pc1_G = pc1_G[pc1_G$V1 == "chr1",]

model_G$cor = pc1_G[match(model_G$pos,pc1_G$V2/100000), 4]
model_G = na.omit(model_G)
model_G[model_G$cor > 1, 4] <- 1
model_G[model_G$cor < -1, 4] <- -1
plot_ly(model_G, x = ~V1, y = ~V2, z = ~-V3, color = ~cor, colors = c("blue","red"))  %>% 
  layout(scene=scene)

#
model_G =  hic_model("../hicpro/MK1/iced/100000/MK1_100000_iced.matrix","../hicpro/CLP1/raw/100000/CLP1_100000_abs.bed", "chr1")
model_G$pos = floor(as.numeric(rownames(model_G))/5)*5
plot_ly(model_G, x = ~V1, y = ~V2, z = ~-V3, color = ~pos)  %>% 
  layout(scene=scene)

pc1_G = read.table("../homer/MK1.PC1.bedGraph", skip = 1)
pc1_G = pc1_G[pc1_G$V1 == "chr1",]

model_G$cor = pc1_G[match(model_G$pos,pc1_G$V2/100000), 4]
model_G = na.omit(model_G)
model_G[model_G$cor > 1, 4] <- 1
model_G[model_G$cor < -1, 4] <- -1
plot_ly(model_G, x = ~V1, y = ~V2, z = ~-V3, color = ~cor, colors = c("blue","red"))  %>% 
  layout(scene=scene)

#
model_G =  hic_model("../hicpro/MPP2/iced/100000/MPP2_100000_iced.matrix","../hicpro/CLP1/raw/100000/CLP1_100000_abs.bed", "chr1")
model_G$pos = floor(as.numeric(rownames(model_G))/5)*5
plot_ly(model_G, x = ~V1, y = ~V2, z = ~-V3, color = ~pos)  %>% 
  layout(scene=scene)

pc1_G = read.table("../homer/MPP2.PC1.bedGraph", skip = 1)
pc1_G = pc1_G[pc1_G$V1 == "chr1",]

model_G$cor = pc1_G[match(model_G$pos,pc1_G$V2/100000), 4]
model_G = na.omit(model_G)
model_G[model_G$cor > 1, 4] <- 1
model_G[model_G$cor < -1, 4] <- -1
plot_ly(model_G, x = ~V1, y = ~V2, z = ~-V3, color = ~cor, colors = c("blue","red"))  %>% 
  layout(scene=scene)



# stat
stat_pos = data.frame(row.names = c("sample","chr", "a", "b"))
n=0
for (i in c("LT1", "ST1","MPP2", "CLP1","CMP1","GMP1","MEP1","MKP1","MK1","G1")) {
  print(i)
  for(j in paste0("chr", c(seq(1:19), "X"))){
    n=n+1
    print(j)
    pos = hic_model(paste0("../hicpro/",i,"/iced/500000/",i,"_500000_iced.matrix"),"../hicpro/CLP1/raw/500000/CLP1_500000_abs.bed", j)
    stat_pos[[n]] = c(i, j, max(pos$V1)-min(pos$V1),  max(pos$V2)-min(pos$V2))
   }
}

stat_pos.t = data.frame(t(stat_pos), stringsAsFactors = F)
stat_pos.t$ratio = as.numeric(stat_pos.t$a) / as.numeric(stat_pos.t$b)
stat_pos.t$sample = factor(stat_pos.t$sample, levels = c("LT1", "ST1","MPP2", "CLP1","CMP1","GMP1","MEP1","MKP1","MK1","G1"))
stat_pos.t$chr = factor(stat_pos.t$chr, levels = paste0("chr", c(seq(1:19), "X")))


saveRDS(stat_pos.t, file="Modeling_stat_pos.rds")
stat_pos.t = readRDS("Modeling_stat_pos.rds")

ggplot(stat_pos.t[stat_pos.t$chr %in%  paste0("chr", c(seq(1:19))), ], aes(sample, ratio))+
  geom_boxplot(outlier.size=0)+
  geom_jitter(width=0.20,aes(col=chr)) +
  theme_minimal(base_size = 20) +
  scale_colour_manual(values=rainbow(25)[1:20])+xlab("")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
  #scale_colour_manual(values=colorRampPalette(brewer.pal(n = 9, "Paired"))(20))
ggsave("../figure/Modeling_abRation.pdf", device = "pdf", width = 6, height = 3)

barplot(rep(1,20), col=rainbow(25)[1:20],border = NA,axes =F, space=0)


library(vegan)
#### all chromosome (2Mb)
  mat = read.table("../analysis/matrix/homer_2M_balance/LT1.2M.txt", stringsAsFactors = F,sep="\t", header = T, row.names = 1)
  mat = mat[,-1]

  #mat[is.na(mat)] <- 0
  chroms =  sapply(strsplit(rownames(mat), "-"), "[", 1)
  mat = mat[chroms%in%(paste0("chr", c(seq(1:19)))), chroms%in%(paste0("chr", c(seq(1:19))))]
  chroms_first =  sapply(strsplit(rownames(mat), "-"), "[", 2)
  mat = mat[chroms_first != "0", chroms_first != "0",]
  # distance matrix
  mat_min = min(mat[mat!=0])
  chroms =  sapply(strsplit(rownames(mat), "-"), "[", 1)
  mat[mat == 0] <- mat_min
  colnames(mat) = rownames(mat)
  
  q = quantile(as.matrix(mat),prob=0.98)
  mat1 = mat
  mat1[mat1>q] <- q
  image(as.matrix((mat1)),col = colorRampPalette(c("blue","white", "red"))(n = 20))
  

  C = (mat1)
  chroms =  sapply(strsplit(rownames(C), "-"), "[", 1)
  
  cLAD = read.table("../data/cLAD.GSE17051.mm10.2mb.bed", stringsAsFactors = F)
  cLAD = unique(paste(cLAD$V1,cLAD$V2,sep = "-"))
  
  cLAD = read.table("../homer/LT1.PC1.bedGraph", skip = 1)
  rownames(cLAD) = paste(cLAD$V1, cLAD$V2, sep="-")
  for(i in 1:nrow(C)){
    for(j in 1:ncol(C)){
      if(chroms[i] != chroms[j]){
        # cLAD
        if(rownames(C)[i] %in% cLAD && colnames(C)[j] %in% cLAD){
          C[i,j] = C[i,j]/10
        }
      }
    }
  }

  #diag(C) <- max(C)
  image(as.matrix(log(C[1:90,1:90])),col = colorRampPalette(c("white","orange", "red"))(n = 20))
  image(as.matrix(log10(1/(C[1:345,1:345]))),col = colorRampPalette(c("blue","white", "red"))(n = 20))
  
  #C = 1/mat
  #fit = cmdscale((log(1/(C))), eig=T, k=3)
  fit = isoMDS(log(1/as.matrix(C)), k=3)
  pos = as.data.frame(fit$points)
  #pos$cor = 1:nrow(pos)
  pos$cor = sapply(strsplit(rownames(pos), "[-]"), "[", 1)
  plot_ly(pos, x = ~V1, y = ~V2, z = ~V3, color = ~cor,marker = list(size = 10))  %>% 
    layout(scene=scene)
  
  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  scene = list(
    xaxis = ax,
    yaxis = ax,
    zaxis = ax)
  
  plot_ly(pos, x = ~V1, y = ~V2, z = ~V3, color = ~cor)  %>% 
    layout(scene=scene)
  
  
  pc1_LT = read.table("../homer/LT1.PC1.bedGraph", skip = 1)
  rownames(pc1_LT) = paste(pc1_LT$V1, pc1_LT$V2, sep="-")
  pos$AB = pc1_LT[rownames(pos),4]
  plot_ly(pos, x = ~V1, y = ~V2, z = ~V3, color = ~AB,  colors = c("blue","orange","red"), marker = list(size = 10))  %>% 
 layout(scene=scene)

  
  pos$cor1 = "chromosome"
  pos[1,"cor1"] = "centromere"
  for(i in 2:nrow(pos)){
    if(pos[i,"cor"] != pos[i-1,"cor"]){
      pos[i,"cor1"] = "telomere"
      pos[i-1,"cor1"] = "centromere"
    }
  }
  plot_ly(pos, x = ~V1, y = ~V2, z = ~V3, color = ~cor1,colors = c("green",rgb(0,0,0,alpha=0.2),"orange"),  type = 'scatter3d')  %>% 
    layout(scene=scene)
  

  