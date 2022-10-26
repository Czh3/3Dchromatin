# GWAS loci at loop anchor
anchor_gwas = read.table("../analysis/loops_homer_Juicer/intersect.loops.anchor.GWAS.bed")
anchor_gwas = anchor_gwas[,1:9]
anchor_gwas = unique(anchor_gwas)
anchor_gwas_ratio = as.data.frame(table(anchor_gwas$V8) )
anchor_gwas_ratio$T1 = nrow(anchor_gwas)

gwas_database = read.table("../analysis/loops_homer_Juicer/gwas_catalog_v1.0.2-associations_e96_r2019-10-14.mm10.bed", stringsAsFactors = F)
gwas_database_ratio = as.data.frame(table(gwas_database$V8) )
gwas_database_ratio$T2 = nrow(gwas_database)

anchor_gwas_OE = merge(anchor_gwas_ratio,gwas_database_ratio, by="Var1")
anchor_gwas_OE$p = apply(anchor_gwas_OE[,-1],1,function(x){chisq.test(matrix(c(x[1],x[2]-x[1],x[3],x[4]-x[3]),nrow=2), simulate.p.value = TRUE, B = 10000)$p.value})
anchor_gwas_OE$FDR = p.adjust(anchor_gwas_OE$p, method = "fdr")
anchor_gwas_OE = anchor_gwas_OE[order(anchor_gwas_OE$FDR),]


anchor_gwas_OE1 = anchor_gwas_OE[1:15,]
anchor_gwas_OE1$Var1  = gsub("_"," ",as.character(anchor_gwas_OE1$Var1))
anchor_gwas_OE1$Var1 = factor(anchor_gwas_OE1$Var1, levels=rev(anchor_gwas_OE1$Var1))
ggplot(anchor_gwas_OE1[1:10,], aes(Var1,-log10(p) )) +
  geom_bar(stat = "identity", fill="black",width=0.3) + coord_flip()+
  theme_classic() + xlab("") + ylab("-log10(p-value)") +
  theme(axis.text = element_text(size=10, color="black"))+
  scale_colour_brewer(palette="Dark2")
ggsave("../figure/GWAS/loops_anchor_gwas.pdf", device = "pdf", width = 6, height = 4)

ggplot(anchor_gwas_OE1, aes(Var1,-log10(p) )) +
  geom_point(stat = "identity", color="darkred", size=3) + coord_flip()+
  theme_classic() + xlab("") + ylab("-log10(p-value)") +
  theme(axis.text = element_text(size=10, color="black"))+
  scale_colour_brewer(palette="Dark2")
ggsave("../figure/GWAS/loops_anchor_gwas1.pdf", device = "pdf", width = 5, height = 4)




#1 3d Venn
g1 = read.table("../analysis/loops_homer_Juicer/1D.genelist",stringsAsFactors = F)$V1
g3 = read.table("../analysis/loops_homer_Juicer/3D.genelist",stringsAsFactors = F)$V1
require(VennDiagram)
venn.diagram(list("1D" = g1, "3D" = g3),fill = c("red", "green"),height=1500,width=1500,
             alpha = c(0.5, 0.5),lty =1,cex=2,cat.cex=2,cat.pos=0,
             filename = "../figure/GWAS/1d_3d.gene.venn.tiff")



# GO

library(goseq)
library(org.Mm.eg.db)
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
traits = c("Blood_protein_levels", "Red_blood_cell_count", "White_blood_cell_count", "Eosinophil_counts", "Inflammatory_bowel_disease", "Granulocyte_percentage_of_myeloid_white_cells",
           "Mean_corpuscular_hemoglobin", "Systemic_lupus_erythematosus", "Type_1_diabetes", "Body_mass_index")
traits = c("Systemic_lupus_erythematosus", "Type_1_diabetes", "Inflammatory_bowel_disease")
traits_GO_3D = data.frame()
traits_GO_3D_plot = c()
for (i in traits){
  ann_1D_3D = read.table("../analysis/loops_homer_Juicer/intersect.loops.anchor.GWAS.TSS_1D3D.bed", sep="\t", stringsAsFactors = F)
  ann_1D_3D = ann_1D_3D[ann_1D_3D$V8 == i,]
  ann_1D_gene = na.omit(unique(unlist(strsplit(ann_1D_3D$V17, ","))))
  ann_3D_gene = na.omit(unique(unlist(strsplit(ann_1D_3D$V20, ","))))
  
  genes = rep(0, nrow(expr))
  names(genes) = rownames(expr)
  genes[names(genes) %in% ann_1D_gene] = 1
  pwf=nullp(genes,"mm10","geneSymbol", plot.fit = FALSE)
  GO_1D=goseq(pwf,"mm10","geneSymbol")
  GO_1D = GO_1D[GO_1D$ontology == "BP",]
  GO_1D$trait = paste0(i,"_1D")
  traits_GO_3D = rbind(traits_GO_3D, GO_1D)
  traits_GO_3D_plot = c(traits_GO_3D_plot, GO_1D[1:3, "term"])
  
  genes = rep(0, nrow(expr))
  names(genes) = rownames(expr)
  genes[names(genes) %in% ann_3D_gene] = 1
  pwf=nullp(genes,"mm10","geneSymbol", plot.fit = FALSE)
  GO_3D=goseq(pwf,"mm10","geneSymbol")
  GO_3D = GO_3D[GO_3D$ontology == "BP",]
  GO_3D$trait =  paste0(i,"_3D")
  traits_GO_3D = rbind(traits_GO_3D, GO_3D)
  traits_GO_3D_plot = c(traits_GO_3D_plot, GO_3D[1:3, "term"])
}
traits_GO_3D_heatmap = traits_GO_3D[traits_GO_3D$term %in% traits_GO_3D_plot, ]
ggplot(traits_GO_3D_heatmap, aes(trait, term, size=numDEInCat/numInCat, color=-log10(over_represented_pvalue)))+
  geom_point()


genes = rep(0, nrow(expr))
names(genes) = rownames(expr)
genes[names(genes) %in% ann_1D_gene] = 1
pwf=nullp(genes,"mm10","geneSymbol", plot.fit = FALSE)
GO_1D=goseq(pwf,"mm10","geneSymbol")
GO_1D = GO_1D[GO_1D$ontology == "BP",]

genes = rep(0, nrow(expr))
names(genes) = rownames(expr)
genes[names(genes) %in% ann_3D_gene] = 1
pwf=nullp(genes,"mm10","geneSymbol", plot.fit = FALSE)
GO_3D=goseq(pwf,"mm10","geneSymbol")
GO_3D = GO_3D[GO_3D$ontology == "BP",]



# Myc expr
# gene expr
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(2,3,5:7,9,12,14,16:18,20,25,27,31,32)]

expr.myc = melt(expr["Myc",])
expr.myc$sample = toupper(gsub(".1|.2","",expr.myc$variable))
expr.myc$sample = factor(expr.myc$sample, levels=c("LT.HSC", "ST.HSC", "MPP", "CLP","CMP","GMP","MEP","MKP","MK","GR"))
expr.myc$rep = stringr::str_sub(expr.myc$variable, -1, -1)

ggplot(expr.myc, aes(sample, value, fill=rep))+
  geom_bar(stat="identity",position=position_dodge())+
  theme_classic(base_size = 20) + ylab("TPM") + xlab("")+
  scale_fill_manual(values=c("black", "darkgray"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black", angle = 30, hjust = 1,vjust = 1))
ggsave("../figure/GWAS/Myc_expr_barplot.pdf", device = "pdf", width = 5, height = 4)


expr.myc_select = expr.myc[expr.myc$sample %in% c("LT.HSC", "ST.HSC", "CMP", "GMP", "GR"),]
ggplot(expr.myc_select, aes(sample, log(value+1), color=rep, group=rep))+
  geom_line(size=1)+
  geom_point(aes(shape=rep), size=3)+
  theme_classic(base_size = 20) + ylab("log(TPM+1)") + xlab("")+
  scale_color_manual(values=c("black", "gray60"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text( color="black", angle = 30, hjust = 1,vjust = 1))
ggsave("../figure/GWAS/Myc_expr_lineplot.pdf", device = "pdf", width = 5, height = 3)



# GO of specific straits





