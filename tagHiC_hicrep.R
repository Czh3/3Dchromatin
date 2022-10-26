# hic-rep 1.11.1


remotes::install_github("qunhualilab/hicrep")
library(hicrep)
remotes::install_github("aidenlab/straw/R")
library("strawr")
install.packages("stat")


GM_1k_rep1 <- hic2mat("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/GM12878_1k_rep1.hic/GM12878_11.allValidPairs.hic", chromosome1 = "1", chromosome2 = "1", resol = 40000, method = "NONE") 
GM_1k_rep2 <- hic2mat("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/GM12878_1k_rep2.hic/GM12878_12.allValidPairs.hic", chromosome1 = "1", chromosome2 = "1", resol =   40000, method = "NONE") 
GM_1m <- hic2mat("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/GM12878_bulk_lumeng/bulk_allValidPairs.hic", chromosome1 = "1", chromosome2 = "1", resol =   40000, method = "NONE") 

hic_files = list.files("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/", "*hic", recursive=T)

scc = matrix(nrow=0,ncol = 3)
for(i in hic_files){
  for(j in hic_files){
    f1 = hic2mat(paste0("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/",i), chromosome1 = "1", chromosome2 = "1", resol = 40000, method = "NONE") 
    f2 = hic2mat(paste0("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/",j), chromosome1 = "1", chromosome2 = "1", resol = 40000, method = "NONE") 
    scc.out = get.scc(f1, f2, resol = 40000, h = 6, lbr = 100000, ubr = 50000000)
    scc = rbind(scc, c(i,j,scc.out[[3]]))
  }
}


