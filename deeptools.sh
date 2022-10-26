grep "0-1" ../data/ST_gene_expr_interation.tab | > ../data/expr_0_1.gene.txt
grep "1-10" ../data/ST_gene_expr_interation.tab > ../data/expr_1_10.gene.txt
grep "10-30" ../data/ST_gene_expr_interation.tab > ../data/expr_10_30.gene.txt
grep "30-80" ../data/ST_gene_expr_interation.tab > ../data/expr_30_80.gene.txt
grep ">80" ../data/ST_gene_expr_interation.tab > ../data/expr_80.gene.txt


computeMatrix scale-regions -S /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/GSM1441286_H3K4me1ST_HSC.ucsc.bw /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/GSM1441302_H3K4me2ST_HSC.ucsc.bw /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/GSM1441318_H3K4me3ST_HSC.ucsc.bw /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/GSM1441270_K27Ac_ST_HSC.ucsc.bw -R ../data/ST_gene_expr_interation.tab.deeptools -o ../data/ST_gene_expr_interation.tab.deeptools.mat.gz -p 12 -b 10000 -a 10000 --regionBodyLength 20000

plotProfile -m ../data/ST_gene_expr_interation.tab.deeptools.mat.gz -out ../figure/ST_gene_expr.epigenome.pdf --numPlotsPerRow 1 --yMax 10 30 50 30
#plotHeatmap -m ../data/ST_gene_expr_interation.tab.deeptools.mat.gz -out ../figure/ST_gene_expr.epigenome.heatmap.pdf --colorMap RdBu --whatToShow 'heatmap and colorbar' 


computeMatrix scale-regions -S /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/GSM1441286_H3K4me1ST_HSC.ucsc.bw /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/GSM1441302_H3K4me2ST_HSC.ucsc.bw /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/GSM1441318_H3K4me3ST_HSC.ucsc.bw /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/GSM1441270_K27Ac_ST_HSC.ucsc.bw -R ../data/ST_GAD.tab.deeptools -o ../data/ST_GAD.tab.deeptools.mat.gz -p 12 -b 10000 -a 10000 --regionBodyLength 20000
plotProfile -m ../data/ST_GAD.tab.deeptools.mat.gz -out ../figure/GAD/ST_GAD.deeptools.epigenome.pdf --numPlotsPerRow 1 --yMax 10 30 50 30 --colors darkred darkgray

# MEL cells
computeMatrix scale-regions -S /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF488XKH_MEL_H3K4me1.bigWig /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF515BPN_MEL_CTCF.bigWig /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF198FKD_MEL_RAD21.bigWig /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF906BFI_MEL_POLR2A.bigWig /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF034BBD_MEL_SMC3.bigWig -R ../data/ST_gene_expr_interation.tab.deeptools -o ../data/ST_gene_expr_interation.tab.deeptools.mat.ctcf.gz  -p 12 -b 10000 -a 10000 --regionBodyLength 20000
plotProfile -m ../data/ST_gene_expr_interation.tab.deeptools.mat.ctcf.gz -out ../figure/ST_gene_expr.ctcf.pdf --numPlotsPerRow 1 --yMax 2 2.5 2.5 5 2 --yMin 0 

computeMatrix scale-regions -S /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF488XKH_MEL_H3K4me1.bigWig /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF515BPN_MEL_CTCF.bigWig /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF198FKD_MEL_RAD21.bigWig /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF906BFI_MEL_POLR2A.bigWig /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF034BBD_MEL_SMC3.bigWig -R ../data/ST_GAD.tab.deeptools -o ../data/ST_GAD.tab.deeptools.mat.ctcf.gz -p 12 -b 10000 -a 10000 --regionBodyLength 20000
plotProfile -m  ../data/ST_GAD.tab.deeptools.mat.ctcf.gz -out ../figure/GAD/ST_GAD.deeptools.epigenome.ctcf.pdf --numPlotsPerRow 1 --yMax 2 2.5 2.5 5 2 --yMin 0 --colors darkred darkgray

