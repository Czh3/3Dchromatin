mat=`ls ../hicpro/*/iced/40000/*_40000_iced.matrix`
lab="CLP1 CLP2 CMP1 CMP2 G1 G2 GMP1 GMP2 LT1 LT2 MEP1 MEP2 MK1 MK2 MKP1 MKP2 MPP1 MPP2 ST1 ST2"
abs="/lustre/user/liclab/zhangc/proj/blood3D/HiC/hic_pro_out/hic_results/matrix/cmp/raw/40000/cmp_40000_abs.bed"
#chr1:58-70mb
/lustre/user/liclab/software/anaconda3/bin/python /lustre/user/liclab/zhangc/Taolab/guan/tools/HiCPlotter-master/HiCPlotter.py -tri 1   \
    -f /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/LT1/iced/40000/LT1_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/LT2/iced/40000/LT2_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/G1/iced/40000/G1_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/G2/iced/40000/G2_40000_iced.matrix -o ../figure/HiCplotter/chr1 -n LT1 LT2 G1 G2 -chr chr1 -s 1450 -e 1750 -bed $abs -r 40000 --plotInsulation 1 -mm 6 -w 8

/lustre/user/liclab/software/anaconda3/bin/python /lustre/user/liclab/zhangc/Taolab/guan/tools/HiCPlotter-master/HiCPlotter.py -tri 1   \
    -f /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/LT1/iced/40000/LT1_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/LT2/iced/40000/LT2_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/G1/iced/40000/G1_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/G2/iced/40000/G2_40000_iced.matrix -o ../figure/HiCplotter/chr1 -n LT1 LT2 G1 G2 -chr chr1 -s 1570 -e 1700 -bed $abs -r 40000 -mm 4 -dpi 300 -hist /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/LT1/allChr.insulation.bedGraph /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/LT2/allChr.insulation.bedGraph /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/G1/allChr.insulation.bedGraph /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/G2/allChr.insulation.bedGraph -fhist 1 1 1 1 -hl IS IS IS IS


#chr1:91,769,176-98,628,146
/lustre/user/liclab/software/anaconda3/bin/python /lustre/user/liclab/zhangc/Taolab/guan/tools/HiCPlotter-master/HiCPlotter.py -tri 1   \
    -f /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/LT1/iced/40000/LT1_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/LT2/iced/40000/LT2_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/G1/iced/40000/G1_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/G2/iced/40000/G2_40000_iced.matrix -o ../figure/HiCplotter/chr1 -n LT1 LT2 G1 G2 -chr chr1 -s 2325 -e 2450 -bed $abs -r 40000 -mm 4 -dpi 300 -hist /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/LT1/allChr.insulation.bedGraph /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/LT2/allChr.insulation.bedGraph /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/G1/allChr.insulation.bedGraph /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/G2/allChr.insulation.bedGraph -fhist 1 1 1 1 -hl IS IS IS IS


#chr10:109,416,215-122,030,614
/lustre/user/liclab/software/anaconda3/bin/python /lustre/user/liclab/zhangc/Taolab/guan/tools/HiCPlotter-master/HiCPlotter.py -tri 1   \
    -f /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/LT1/iced/40000/LT1_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/CMP1/iced/40000/CMP1_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/MK1/iced/40000/MK1_40000_iced.matrix /lustre/user/liclab/zhangc/proj/blood3D_3/hicpro/G1/iced/40000/G1_40000_iced.matrix -o ../figure/HiCplotter/chr10 -n LT CMP MK G -chr chr10 -s 2740 -e 3050 -bed $abs -r 40000 -mm 4 -dpi 300 -hist /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/LT1/allChr.insulation.bedGraph /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/CMP1/allChr.insulation.bedGraph /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/MK1/allChr.insulation.bedGraph /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/G1/allChr.insulation.bedGraph -fhist 1 1 1 1 -hl IS IS IS IS -hc 9ACD32 9ACD32 9ACD32 9ACD32


# Hoxa 
abs="/lustre/user/liclab/zhangc/proj/blood3D/HiC/hic_pro_out/hic_results/matrix/cmp/raw/10000/cmp_10000_abs.bed" 
for i in `ls ../hicpro/*/iced/10000/*_10000_iced.matrix|grep -v MK` 
do
{
	lab=`basename $i`
	lab=`echo ${lab/_10000_iced.matrix/}`
	sample=`echo ${lab/[1|2]/}`
	H3K4me1=`ls /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/*_H3K4me1$sample*.ucsc.sorted.bgr`
	H3K4me1=`echo $H3K4me1|cut -f 2 -d " "`
	H3K27ac=`ls /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/*K27Ac_$sample*.ucsc.sorted.bgr` 
	H3K27ac=`echo $H3K27ac|cut -f 2 -d " "` 
	#sample1=`echo $sample | tr '[:upper:]' '[:lower:]'`
	#RNA=`ls /lustre/user/liclab/zhangc/proj/blood3D/RNA_seq/bigwig/$sample1.bgr`
	#RNA=`echo $RNA|cut -f 2 -d " "`  
	/lustre/user/liclab/software/anaconda3/bin/python /lustre/user/liclab/zhangc/Taolab/guan/tools/HiCPlotter-master/HiCPlotter.py -tri 1 	\
		-f $i -o ../figure/HiCplotter/Hoxa.$lab -n $lab 	\
		-s 5170 -e 5270 -chr chr6 -r 10000 	\
		-bed $abs -g  /lustre/user/liclab/zhangc/Taolab/guan/hic_data_rep1/hic_pro_result/plot/mm10.all.sort.gene 	\
		-hist $H3K4me1,$H3K27ac -fhist 1,1 -hl H3K4me1,H3K27ac -hm 30,60
}&
done


