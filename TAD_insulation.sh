lab=`dirname /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/*/*.bw | cut -f 10 -d "/"`
echo $lab

#computeMatrix reference-point -S /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/*/*.bw -R /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/LT1/LT1.merge.bed -b 600000 -a 600000 -o ../data/InsulationScore_at_LT1_boundary.gz -p 10 -bs 40000 
#plotProfile -m ../data/InsulationScore_at_LT1_boundary.gz -out ../figure/InsulationScore_at_LT1_boundary.png --perGroup --samplesLabel $lab --refPointLabel LT_TAD_boundary --outFileNameData ../data/InsulationScore_at_LT1_boundary.tab


computeMatrix reference-point -S /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/*/*.bw -R /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/TAD/G1/G1.merge.bed -b 600000 -a 600000 -o ../data/InsulationScore_at_G1_boundary.gz -p 10 -bs 40000
plotProfile -m ../data/InsulationScore_at_G1_boundary.gz -out ../figure/InsulationScore_at_G1_boundary.png --perGroup --samplesLabel $lab --refPointLabel G1_TAD_boundary --outFileNameData ../data/InsulationScore_at_G1_boundary.tab
