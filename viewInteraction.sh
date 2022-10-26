
# Meis1
#./straw NONE ../Juicer/merge/G.hic 11:18879817:19018985 11:18000000:19000000 BP 5000 | awk '{if($2-$1>=30000)print}'| awk '{if($3>=3)print "chr11\t"$1"\t"$1+100"\tchr11\t"$2"\t"$2+100"\t1"}' 

for i in `ls ../Juicer/merge/*.hic`
do
{
	#./straw NONE $i 11:18879817:19018985 11:16000000:20000000 BP 5000 | awk '{if($2-$1>=20000)print}'| awk '{if($3>=3)print "chr11\t"$1"\t"$2+100"\t.\t"$3"\t"$3"\tR\t#7A67EE\tchr11\t"$1"\t"$1+100"\t.\t.\tchr11\t"$2"\t"$2+100"\t.\t+"}' | sort -k1,1 -k2,2n  > G.meis1.interaction

	name=`basename $i|cut -f 1 -d "."`
	bedToBigBed -as=interact.as -type=bed5+13  G.meis1.interaction /lustre/user/liclab/publicData/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/chrom_mm10.sizes ../data/interaction/$name.meis1.interaction.bb 
	rm G.meis1.interaction
}
done







