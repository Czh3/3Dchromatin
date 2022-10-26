for i in `ls ../analysis/loops/*.loop/merged_loops.bedpe`
do
{
	name=`dirname $i`
	name=`basename $name`
	grep -v score $i | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3;print "chr"$4,$5,$6}' > ${i/.bedpe/.bed}
	findMotifsGenome.pl ${i/.bedpe/.bed} mm10 ${i/.bedpe/.Motif} -p 5
}
done



# loop anchor CTCF

computeMatrix reference-point -S /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF515BPN_MEL_CTCF.bigWig /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/ENCFF034BBD_MEL_SMC3.bigWig -R /lustre/user/liclab/zhangc/proj/blood3D_3/analysis/loops/merge.loopsAnchor.deeptools.bed -o ../data/LoopAnchor.deeptools.mat.ctcf.gz  -p 12 --referencePoint center  -b 50000 -a 50000 --binSize 200
plotProfile -m ../data/LoopAnchor.deeptools.mat.ctcf.gz -out ../figure/loops/LoopAnchor.deeptools.ctcf.pdf --numPlotsPerRow 1 --yMax 1.3 1.3 --yMin 0 --colors darkred darkgray --refPointLabel anchor --plotWidth 7


findMotifsGenome.pl ../analysis/loops/merge.loopsAnchor.bed mm10 ../analysis/loops/merge.loopsAnchor.Motif -p 20

intersectBed -a /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/cell_2019/ImmGenATAC18_AllOCRsInfo.3.bed -b ../analysis/loops/merge.loopsAnchor.bed > ../analysis/loops/merge.loopsAnchor.ATAC.bed
findMotifsGenome.pl ../analysis/loops/merge.loopsAnchor.ATAC.bed mm10 ../analysis/loops/merge.loopsAnchor.ATAC.Motif

# cluster 6
# homer
for i in `ls ../analysis/loops/Cluster6/*loops`
do
{
	intersectBed -a /lustre/user/liclab/zhangc/proj/blood3D/epigenomeData/cell_2019/ImmGenATAC18_AllOCRsInfo.3.bed -b $i > ${i/.loops/.ATACpeak}
	findMotifsGenome.pl ${i/.loops/.ATACpeak} mm10 ${i/.loops/.Motif} -p 5
}
done

# meme
for i in `ls ../analysis/loops/Cluster6/*.loops`
do
{
	awk '{print $1"\t"int(($2+$3)/2)-10000"\t"int(($2+$3)/2)+10000}' $i | bedtools getfasta -fi /lustre/user/liclab/publicData/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed - -fo ${i/.loops/.fasta}
	#python segBed.py $i | bedtools getfasta -fi /lustre/user/liclab/publicData/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed - -fo ${i/.loops/.fasta}
	centrimo --ethresh 2000 --o ${i/.loops/.meme} ${i/.loops/.fasta} /home/zhangc/tools/MEME/db/motif_databases/MOUSE/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme &
}
done




