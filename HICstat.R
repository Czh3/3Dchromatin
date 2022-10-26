setwd("/lustre/user/liclab/zhangc/proj/blood3D_3/script")

library(gridExtra)
library(reshape2)

stat = read.table("/lustre/user/liclab/zhangc/proj/blood3D_3/analysis/stat/stat.txt", header = T)
stat = stat[stat$sample != "MPP1", ]
stat_ct = stat[,c(1,5,7,8)]
stat_ct = melt(stat_ct)
stat_ct$sample = factor(stat_ct$sample, levels = c("LT1","LT2","ST1","ST2","MPP2","CLP1","CLP2",
                                                   "CMP1","CMP2","MEP1","MEP2","GMP1","GMP2","MKP1","MKP2","MK1","MK2","G1","G2"))

pdf("../figure/cis_trans_interactions_stat.pdf", 9, 5)
ggplot(stat_ct, aes(sample, value/1e6, fill = variable)) +
  geom_bar(stat="identity") +
  theme_classic(base_size = 25) + xlab("") + ylab("Interactions (Million)") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black", angle = 90, hjust = 1,vjust = 0.5))
dev.off()

pdf("../figure/cis_trans_interaction_proportion.pdf", 10, 4)
ggplot(stat_ct, aes(sample, value, fill = variable)) +
  geom_bar(stat="identity", position = "fill") +
  theme_classic(base_size = 25) + xlab("") + ylab("interactions") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black", angle = 90, hjust = 1,vjust = 0.5))
dev.off()

stat = read.table("/lustre/user/liclab/zhangc/proj/blood3D_3/analysis/stat/stat.txt", header = T)
stat = stat[stat$sample != "MPP1", ]
stat_ct = stat[,c(1,5,6)]
stat_ct$tran_cis = stat_ct$trans_interaction / stat_ct$cis_interaction
stat_ct = melt(stat_ct[,-c(2,3)])
stat_ct$sample = factor(stat_ct$sample, levels = c("LT1","LT2","ST1","ST2","MPP2","CLP1","CLP2",
                                                   "CMP1","CMP2","MEP1","MEP2","GMP1","GMP2","MKP1","MKP2","MK1","MK2","G1","G2"))

stat_ct$cell = stringr::str_sub(stat_ct$sample, 1, -2)
stat_ct$rep = stringr::str_sub(stat_ct$sample, -1, -1)
stat_ct$cell = factor(stat_ct$cell, levels=c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G"))
# remove MPP1

ggplot(stat_ct, aes(cell, value, color = cell)) +
  geom_point(aes(shape=rep), size=4) +
  theme_classic(base_size = 20) + xlab("") + ylab("Trans/Cis") +
  scale_color_brewer(palette="Paired") +
  theme(legend.position = "none",
      axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black", angle = 90, hjust = 1,vjust = 0.5))
ggsave("../figure/trans_cis_ratio.pdf", device = "pdf", width = 5, height = 4)



# # # # # Freq along distance

#
FreqD = read.table("../analysis/FreqDistance/FreDistribution.txt", header = 1, sep = "\t", stringsAsFactors = F)
colnames(FreqD) = c("Distance", "CLP1","CLP2","CMP1","CMP2", "G1","G2","GMP1","GMP2","LT1","LT2","MEP1","MEP2","MK1", "MK2", "MKP1","MKP2", "MPP1","MPP2","ST1","ST2")
FreqD = FreqD[,colnames(FreqD) != "MPP1"]
FreqD = FreqD[11:300000,] # remove diag
FreqD$Distance = as.numeric(FreqD$Distance)
FreqD$Dist = log10(FreqD$Distance)
FreqD$Dist = as.numeric(sprintf("%.1f", FreqD$Dist))
#FreqD_colS = colSums(FreqD[,3:21])
#FreqD = t(apply(FreqD, 1, function(x) x/c(1,1,FreqD_colS)))
library(dplyr)
FreqD = FreqD %>% group_by(Dist) %>%
  summarise_each(mean)
FreqD.m = melt(as.data.frame(FreqD[-1,-2]), id="Dist")
FreqD.m$sample = gsub("[1|2]","",FreqD.m$variable)
FreqD.m$sample = factor(FreqD.m$sample, level=c("LT", "ST", "MPP", "CLP","CMP", "GMP","MEP","MKP","MK","G"))


ggplot(FreqD.m, aes(Dist, log10(value), group = variable, color = sample)) +
  geom_line()+ylim(-8,-3)+ xlim(5,8.5)+
  theme_classic(base_size = 25) + xlab("log10(Genomic Distance)") + ylab("Contact Probability") +
  scale_color_d3() +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black"))
ggsave("../figure/Freq_dist1.pdf", device = "pdf", width = 6, height = 4)


ggplot(FreqD.m, aes(Dist, log10(value), group = variable, color = sample)) +
  geom_line( alpha=0.8)+ ylim(-7,-4.8)+ xlim(6.6,8.2)+
  theme_void(base_size = 25) + ylab("")+xlab("")+
  scale_color_d3()+
  theme(legend.position = "none")
ggsave("../figure/Freq_dist_zoomIN.pdf", device = "pdf", width = 3.5, height = 3)


mean(FreqD.m[FreqD.m$Dist==7.2 & FreqD.m$sample %in% c("MEP", "MK", "G"),3])/mean(FreqD.m[FreqD.m$Dist==7.2 & !FreqD.m$sample %in% c("MEP", "MK", "G"),3])

# by chromosome
FreqDis_chr = function(homer_background){
  FreqD = read.table(homer_background,row.names = 1, head=T,skip=54, nrow=1923,sep = "\t", stringsAsFactors = F)
  FreqD = FreqD[-1,c(1,seq(4,43,by=2))]
  
  FreqD.m = reshape2::melt(as.data.frame(FreqD), id="Distance")
  FreqD.m$Dist = as.numeric(sprintf("%.1f", log10(FreqD.m$Distance)))
  
  FreqD.m = FreqD.m %>% group_by(Dist,variable) %>%
    summarise_each(mean)
  FreqD.m$value = FreqD.m$value / sum(FreqD.m$value)
  return(FreqD.m)
}

FreqD_chr = matrix(nrow=0, ncol=5)
for(i in c("LT", "ST", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","G")){
  FD = FreqDis_chr(paste0("../homer_merge/",i,"/HiCbackground_100000_bp.txt"))
  FD = as.data.frame(FD)
  FD$sample = i
  FreqD_chr = rbind(FreqD_chr,FD)
}
FreqD_chr$sample = factor(FreqD_chr$sample, level=c("LT", "ST", "MPP", "CLP","CMP", "MEP", "GMP","MKP","MK","G"))

ggplot(FreqD_chr[FreqD_chr$value != 0,], aes(Dist, log10(value), group = sample, color = sample)) +
  geom_line( alpha=0.8)+
  #geom_smooth(method = 'loess', se=F)+
  theme_classic() + xlab("log10(Genomic Distance)") + ylab("Contact Probability") +
  facet_wrap( ~ variable, ncol=5) +
  scale_color_d3() +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))
ggsave("../figure/Freq_dist_byChr.pdf", device = "pdf", width = 7, height = 6)



### AB compartment
LT_A = read.table("../analysis/FreqDistance/G1.FreqDis.txt")
ggplot(LT_A, aes(log10(distance), log10(freq), group = sample, color = sample)) +
  geom_line( alpha=0.8, size=1)+ ylim(-5,-1)+
  theme_classic(base_size = 25) + xlab("log10(Genomic Distance)") + ylab("Contact Probability") +
  scale_color_brewer(palette="Paired") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black"))


samples = c("LT1", "ST1","MPP2", "CLP1","CMP1", "MEP1","GMP1","MKP1","MK1","G1")

mats = data.frame()
for(i in samples){
  mat = read.table(paste0("../analysis/FreqDistance/",i,".FreqDis.txt"))
  mats = rbind(mats, mat)
}
ggplot(mats, aes(log10(distance), log10(freq), group = sample, color=sample))+
  #geom_smooth(method="loess", se=F, fullrange=T) + ylim(-6,-1)
  geom_line(size=1) + ylim(-5,-1) +
  theme_classic(base_size=20) + ggtitle("All")+
  scale_color_brewer(palette="Paired")


mats = data.frame()
for(i in samples){
  mat = read.table(paste0("../analysis/FreqDistance/",i,".A.FreqDis.txt"))
  mats = rbind(mats, mat)
}
mats_A = mats
mats = data.frame()
for(i in samples){
  mat = read.table(paste0("../analysis/FreqDistance/",i,".B.FreqDis.txt"))
  mats = rbind(mats, mat)
}
mats_B = mats
mats_A$com = "A"
mats_B$com = "B"
mats_AB = rbind(mats_A,mats_B)
mats_AB$com = as.factor(mats_AB$com)
ggplot(mats_AB[mats_AB$sample %in% c("LT1","G1","MK1"),], aes(log10(distance), log10(freq),color=sample))+
  #geom_smooth(method="loess", se=F, fullrange=T) + ylim(-6,-1)
  geom_line(aes(linetype=com), size=1) + ylim(-5,-1) +
  theme_classic(base_size=20) + ggtitle("A/B")+
  scale_color_brewer(palette="Paired")
ggplot(mats_AB[], aes(log10(distance), log10(freq),color=sample))+
  #geom_smooth(method="loess", se=F, fullrange=T) + ylim(-6,-1)
  geom_line(aes(linetype=com), size=1) + ylim(-5,-1)+ #xlim(6,8)+
  theme_classic(base_size=20) + ggtitle("A/B")+
  scale_color_brewer(palette="Paired")


mats_B_A = merge(mats_B, mats_A, by=c("distance", "sample"))
ggplot(mats_B_A, aes(log10(distance), log10(freq.y)-log10(freq.x),color=sample))+
  #geom_smooth(method="loess", se=F, fullrange=T) + ylim(-6,-1)
  geom_line( size=1) + xlim(6,8.2)+ 
  theme_classic(base_size=20) + ggtitle("A-B")+
  scale_color_brewer(palette="Paired")


# MNase digestion
MNase = read.table("../data/MNase_digestion.txt", header=T, row.names = 1)
MNase.m = melt(data.matrix(MNase))

ggplot(MNase.m, aes(Var1, value, group=Var2, color=Var2))+
  geom_line(size=1)+
  theme_classic(base_size = 18) + xlab("") + ylab("Gray value")+ 
  theme(axis.text = element_text(colour = "black"))+
  scale_color_manual(values=pal_d3("category10")(10)[5:10])
ggsave("../figure/MNase_digestion_distribution.pdf", device = "pdf", width = 5, height = 4)


MNase = read.table("../data/MNase_digestion.txt", header=T, row.names = 1)
MNase = data.frame(colMeans(MNase[1:80,]))
colnames(MNase) = "Mean gray value"
MNase$sample = factor(rownames(MNase), levels=rownames(MNase))
ggplot(MNase, aes(sample, `Mean gray value`,fill=sample))+
  geom_bar(stat = "identity", size=1)+
  theme_classic(base_size = 18) + xlab("") + #ylab("Gray value")+ 
  theme(axis.text = element_text(colour = "black"))+
  scale_fill_manual(values=pal_d3("category10")(10)[5:10])
ggsave("../figure/MNase_digestion_over2k_mean.pdf", device = "pdf", width = 5, height = 4)





