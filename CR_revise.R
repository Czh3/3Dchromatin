# CR revise


### Tcf3
# RNA-seq
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = as.matrix(expr[,c(12,7,3,2,17,18,32,25)])
colnames(expr) = c("LT","ST","CMP","CLP","MEP","GMP","MK","G")
expr.tcf3 = expr[c("Tcf3","Hmga2","Rcc2","Psat1"),]
expr.tcf3 = melt(expr.tcf3)


ggplot(expr.tcf3, aes(Var2, (value), fill=Var1))+
  geom_bar(stat = "identity") +
  scale_fill_aaas()+
  facet_wrap( ~ Var1, scales = "free_y")+
  theme_classic(base_size = 18) +
  xlab("") + ylab("Expr. (TPM)") + 
  theme(legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black", angle = 30, hjust = 1,vjust = 1)) 
ggsave("../figure/CR/Tcf3_expr.pdf",device = "pdf",width = 9, height = 4)



# Loops motif
LT_loops = read.table("/lustre/user/liclab/zhangc/proj/blood3D_3/analysis/loops/LT.loop/merged_loops.Motif/knownResults.txt",
                      sep = "\t", comment.char = "^", header=T, stringsAsFactors = F)
LT_loops = LT_loops[1:8,]
LT_loops$Motif =  sapply(strsplit(LT_loops$Motif.Name, '/'), "[", 1)
LT_loops$Motif = factor(LT_loops$Motif, levels = rev(LT_loops$Motif))

ggplot(LT_loops, aes(Motif, (-Log.P.value)))+
  geom_bar(stat = "identity") +
  scale_fill_aaas()+
  theme_classic(base_size = 18) +
  xlab("") + ylab("-log(p value)") + coord_flip()+
  theme(legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black")) 

ggsave("../figure/CR/LT_loops_motif.pdf",device = "pdf",width = 5, height = 4)



LT_loops = read.table("/lustre/user/liclab/zhangc/proj/blood3D_3/analysis/loops/ST.loop/merged_loops.Motif/knownResults.txt",
                      sep = "\t", comment.char = "^", header=T, stringsAsFactors = F)
LT_loops = LT_loops[1:8,]
LT_loops$Motif =  sapply(strsplit(LT_loops$Motif.Name, '/'), "[", 1)
LT_loops$Motif = factor(LT_loops$Motif, levels = rev(LT_loops$Motif))

ggplot(LT_loops, aes(Motif, -log10(P.value)))+
  geom_bar(stat = "identity") +
  scale_fill_aaas()+
  theme_classic(base_size = 18) +
  xlab("") + ylab("-log10(p value)") + coord_flip()+
  theme(legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black")) 
ggsave("../figure/CR/ST_loops_motif.pdf",device = "pdf",width = 5, height = 4)


