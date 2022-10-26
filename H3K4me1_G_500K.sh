awk '{if($1=="chr1")print}' ../homer/G1.PC1.bedGraph > tmp


bedtools map -a tmp -b /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/GSE60103/mm10/GSM1441293_H3K4me1GN.ucsc.sorted.bgr -c 4 -o sum > ../data/H3K4me1_G_chr1_500K.bed

rm tmp
