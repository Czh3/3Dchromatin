bedtools sort -i ../data/AB_dynamic.regions.bed > ../data/AB_dynamic.regions.sort.bed

cp ../data/AB_dynamic.regions.sort.bed tmp
for i in GSM1441285_H3K4me1LT_HSC.ucsc.sorted.bgr GSM1441286_H3K4me1ST_HSC.ucsc.sorted.bgr GSM1441287_H3K4me1MPP.ucsc.sorted.bgr GSM1441288_H3K4me1CMP.ucsc.sorted.bgr GSM1441289_H3K4me1GMP.ucsc.sorted.bgr GSM1441290_H3K4me1MEP.ucsc.sorted.bgr GSM1441300_H3K4me1CLP.ucsc.sorted.bgr GSM1441293_H3K4me1GN.ucsc.sorted.bgr
do
{
    bedtools map -a tmp -b /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/$i -c 4 -o sum > tmp1
    mv tmp1 tmp
}
done

mv tmp ../data/AB_dynamic.regions.H3K4me1.bed


cp ../data/AB_dynamic.regions.sort.bed tmp
for i in `ls /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/*K27*sorted.bgr`
do
{
	bedtools map -a tmp -b $i  -c 4 -o sum > tmp1
	mv tmp1 tmp
}
done
mv  tmp ../data/AB_dynamic.regions.H3K27ac.bed

# expr
intersectBed -a ../data/AB_dynamic.regions.bed -b /lustre/user/liclab/zhangc/Taolab/guan/rna-seq/reference/gene_position.all.bed -wa -wb > ../data/AB_dynamic.regions.RNAexpr.bed



# cell 2019
intersectBed -a ../data/AB_dynamic.regions.sort.bed -b /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/cell_2019/ImmGenATAC18_AllOCRsInfo.bed -wa -wb > ../data/AB_dynamic.regions.sort.OCR.bed





## CLP CMP
bedtools sort -i ../data/AB_CLP_CMP.regions.bed > ../data/AB_CLP_CMP.regions.sort.bed
cp ../data/AB_CLP_CMP.regions.sort.bed tmp
for i in GSM1441285_H3K4me1LT_HSC.ucsc.sorted.bgr GSM1441286_H3K4me1ST_HSC.ucsc.sorted.bgr GSM1441287_H3K4me1MPP.ucsc.sorted.bgr GSM1441288_H3K4me1CMP.ucsc.sorted.bgr GSM1441289_H3K4me1GMP.ucsc.sorted.bgr GSM1441290_H3K4me1MEP.ucsc.sorted.bgr GSM1441300_H3K4me1CLP.ucsc.sorted.bgr GSM1441293_H3K4me1GN.ucsc.sorted.bgr
do
{
    bedtools map -a tmp -b /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/$i -c 4 -o sum > tmp1
    mv tmp1 tmp
}
done

mv tmp ../data/AB_CLP_CMP.regions.H3K4me1.bed

intersectBed -a ../data/AB_CLP_CMP.regions.bed -b /lustre/user/liclab/zhangc/Taolab/guan/rna-seq/reference/gene_position.all.bed -wa -wb > ../data/AB_CLP_CMP.regions.RNAexpr.bed

# Granulocytes
bedtools sort -i ../data/AB_Granu.regions.bed > ../data/AB_Granu.regions.sort.bed

cp ../data/AB_Granu.regions.sort.bed tmp
for i in GSM1441285_H3K4me1LT_HSC.ucsc.sorted.bgr GSM1441286_H3K4me1ST_HSC.ucsc.sorted.bgr GSM1441287_H3K4me1MPP.ucsc.sorted.bgr GSM1441288_H3K4me1CMP.ucsc.sorted.bgr GSM1441289_H3K4me1GMP.ucsc.sorted.bgr GSM1441290_H3K4me1MEP.ucsc.sorted.bgr GSM1441300_H3K4me1CLP.ucsc.sorted.bgr GSM1441293_H3K4me1GN.ucsc.sorted.bgr
do
{
    bedtools map -a tmp -b /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/$i -c 4 -o sum > tmp1
    mv tmp1 tmp
}
done
mv tmp ../data/AB_Granu.regions.H3K4me1.bed
intersectBed -a ../data/AB_Granu.regions.sort.bed -b /lustre/user/liclab/zhangc/Taolab/guan/rna-seq/reference/gene_position.all.bed -wa -wb > ../data/AB_Granu.regions.RNAexpr.bed



