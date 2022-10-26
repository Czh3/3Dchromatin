# hic-rep 1.11.1


#remotes::install_github("qunhualilab/hicrep")
library(hicrep)
#remotes::install_github("aidenlab/straw/R")
library("strawr")
library(reshape2)

GM_1k_rep1 <- hic2mat("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/GM12878_1k_rep1.hic/GM12878_11.allValidPairs.hic", chromosome1 = "1", chromosome2 = "1", resol = 500000, method = "KR") 
GM_1k_rep2 <- hic2mat("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/GM12878_1k_rep2.hic/GM12878_12.allValidPairs.hic", chromosome1 = "1", chromosome2 = "1", resol = 500000, method = "KR") 
GM_1m <- hic2mat("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/GM12878_bulk/GM12878.allValidPairs.hic", chromosome1 = "1", chromosome2 = "1", resol =   500000, method = "NONE") 
get.scc(GM_1k_rep1, GM_1k_rep2, resol = 500000, h = 2, lbr = 0, ubr = 5000000)
get.scc(GM_1k_rep1, GM_1m, resol = 500000, h = 2, lbr = 0, ubr = 5000000)


hic_files = list.files("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/", "*hic", recursive=T)

scc = matrix(nrow=0,ncol = 3)
for(i in hic_files){
  for(j in hic_files){
    f1 = hic2mat(paste0("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/",i), chromosome1 = "1", chromosome2 = "1", resol = 500000, method = "KR") 
    f2 = hic2mat(paste0("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/",j), chromosome1 = "1", chromosome2 = "1", resol = 500000, method = "KR") 
    scc.out = get.scc(f1, f2, resol = 500000, h = 3, lbr = 0, ubr = 5000000)
    scc = rbind(scc, c(i,j,scc.out[[3]]))
  }
}

scc = data.frame(scc, stringsAsFactors = F)
scc$X1 = basename(as.character(scc$X1))
scc$X2 = basename(as.character(scc$X2))
scc$X3 = as.numeric(as.character(scc$X3))


scc.matrix = reshape2::dcast(scc, X1~X2)
rownames(scc.matrix) = scc.matrix$X1
scc.matrix = scc.matrix[,-1]
diag(scc.matrix) = NA

pheatmap(as.matrix(scc.matrix[-c(1),-c(1)]),color = colorRampPalette(c("navy", "white", "firebrick3"))(50),display_numbers=T)







