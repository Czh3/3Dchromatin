for i in `ls ../Juicer/merge/*.hic`
do
{
	# Myc
	name=`basename $i|cut -f 1 -d "."`
	./straw KR $i 15:61985000:61985500 15 BP 5000 | awk 'BEGIN{OFS="\t"}{if($1=="61985000"){print "chr15\t"$2-5000"\t"$2"\t"$3}if($2=="61985000"){print "chr15\t"$1-5000"\t"$1"\t"$3}}'>  ../data/interaction/"$name".Myc.4c.bedGraph
	a=`awk '{a=a+$4}END{print a}' ../data/interaction/"$name".Myc.4c.bedGraph`
	awk 'BEGIN{OFS="\t"}{$4=100*$4/'$a';print}' ../data/interaction/"$name".Myc.4c.bedGraph > ../data/interaction/"$name".Myc.4c.bedGraph1
	mv ../data/interaction/"$name".Myc.4c.bedGraph1 ../data/interaction/"$name".Myc.4c.bedGraph 


}
done


