#plot_intergrate_gene_body_interaction
library(reshape2)
# gene length
gene_length = read.table("/lustre/user/liclab/zhangc/Taolab/guan/rna-seq/reference/gene_position.bed")

expr_level = read.csv("../data/ST_gene_expr_interation.tab", row.names = 1)


getMatrix = function(pos){
  #cmd = paste0("./straw KR ../Juicer/merge/ST.hic ",pos," ",pos ," BP 10000 > temp.straw ")
  cmd = paste0("/lustre/user/liclab/zhangc/tools/juicer-master/CPU/common/juicer_tools dump oe KR ../Juicer/merge/ST.hic ",pos," ",pos ,"  BP 5000 > temp.straw")
  system(cmd)
  
  mat = tryCatch(read.table("temp.straw"), error=function(e){return(matrix(0, nrow = 30, ncol = 30))}) 
  
  if(nrow(mat)<=30){
    return(matrix(0, nrow = 30, ncol = 30))
  }
     
  mat = dcast(mat, V1 ~ V2, value.var="V3")
  mat[is.na(mat)] <- 0
  rownames(mat) = mat$V1
  mat = mat[,-1]
  #
  rl = intersect(rownames(mat), colnames(mat))
  mat = mat[rl, rl]
  #
  mat1 = t(mat)
  diag(mat1) = 0
  mat = mat + mat1
  
  l = nrow(mat)
  if(l>=30){
    #s = as.integer(c( seq(1,50,length.out=10), seq(51, l-51, length.out=10), seq(l-50,l,length.out=10) ))
    s = as.integer(seq(1,l,length.out = 30))
    mat = mat[s, s]
    return(as.matrix(mat))
  }else{
    return(matrix(0, nrow = 30, ncol = 30))
  }
}



for(j in levels(expr_level$expr)){
  aggr_mat = matrix(0, nrow = 30, ncol = 30)
  gene_select = row.names(expr_level[expr_level$expr == j, ])
  for (i in gene_select){
    g = gene_length[gene_length$V5 == i,][1,]
    g_len = g$V3-g$V2
    #g_len = 500000
    pos = paste(gsub("chr","",g$V1), g$V2-g_len, g$V3+g_len, sep = ':')
    m = getMatrix(pos)
    aggr_mat = aggr_mat + m
  }
  
  #image(log10(aggr_mat), col=rev(brewer.pal(n = 11, name = 'PiYG')))
  aggr_mat = aggr_mat/length(gene_select)
  breaksList = seq(0.5, 2.5, by = 0.1)
  pdf(paste0("../figure/intra-gene_interaction_expr_level_",j,".pdf"),width=5,height=4.7)
  pheatmap(aggr_mat, show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F, 
           border_color = NA,col=colorRampPalette(rev(brewer.pal(n = 11, name = 'PiYG')))(length(breaksList)),
           breaks = breaksList)
  dev.off()
}


breaksList = seq(0.5, 2.5, by = 0.1)

pheatmap(aggr_mat_80, show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F, 
         border_color = NA,col=colorRampPalette(rev(brewer.pal(n = 11, name = 'PiYG')))(length(breaksList)),
         breaks = breaksList)


aggr_mat_0_1 = aggr_mat/length(gene_select)
pheatmap(aggr_mat_0_1, show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F, 
         border_color = NA,col=colorRampPalette(rev(brewer.pal(n = 11, name = 'PiYG')))(length(breaksList)),
         breaks = breaksList)




# GAD

getMatrix = function(pos){
  cmd = paste0("./straw KR ../Juicer/merge/ST.hic ",pos," ",pos ," BP 5000 > temp.straw ")
  #cmd = paste0("/lustre/user/liclab/zhangc/tools/juicer-master/CPU/common/juicer_tools dump oe KR ../Juicer/merge/ST.hic ",pos," ",pos ,"  BP 5000 > temp.straw")
  system(cmd)
  
  mat = tryCatch(read.table("temp.straw"), error=function(e){return(matrix(0, nrow = 60, ncol = 60))}) 
  
  if(nrow(mat)<=60){
    return(matrix(0, nrow = 60, ncol = 60))
  }
  
  mat = dcast(mat, V1 ~ V2, value.var="V3")
  mat[is.na(mat)] <- 0
  rownames(mat) = mat$V1
  mat = mat[,-1]
  #
  rl = intersect(rownames(mat), colnames(mat))
  mat = mat[rl, rl]
  #
  mat1 = t(mat)
  diag(mat1) = 0
  mat = mat + mat1
  
  l = nrow(mat)
  if(l>=60){
    #s = as.integer(c( seq(1,50,length.out=10), seq(51, l-51, length.out=10), seq(l-50,l,length.out=10) ))
    s = as.integer(seq(1,l,length.out = 60))
    mat = mat[s, s]
    return(as.matrix(mat))
  }else{
    return(matrix(0, nrow = 60, ncol = 60))
  }
}

ST_GAD = read.csv("../data/ST_GAD_expr.tab", row.names = 1)

for(j in levels(ST_GAD$GAD)){
  aggr_mat = matrix(0, nrow = 60, ncol = 60)
  gene_select = row.names(ST_GAD[ST_GAD$GAD == j, ])
  for (i in gene_select){
    g = gene_length[gene_length$V5 == i,][1,]
    g_len = g$V3-g$V2
    pos = paste(gsub("chr","",g$V1), g$V2-g_len, g$V3+g_len, sep = ':')
    m = getMatrix(pos)
    aggr_mat = aggr_mat + m
  }
  
  #image(log10(aggr_mat), col=rev(brewer.pal(n = 11, name = 'PiYG')))
  aggr_mat = aggr_mat/length(gene_select)
  breaksList = seq(-1, 2.5, by = 0.1)
  #pdf(paste0("../figure/GAD/ST_GAD_",j,".intergrate.pdf"),width=5,height=4.7)
  png(paste0("../figure/GAD/ST_GAD_",j,".intergrate.png"),width=1000,height=900,res=300)
  
  pheatmap(log(aggr_mat), show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F, 
           border_color = NA,col=colorRampPalette(c("blue", "white","red","black"))(length(breaksList)),
           breaks = breaksList)
  dev.off()
}

head(gene_length[rownames(ST_GAD[order(ST_GAD$GAD),])])


