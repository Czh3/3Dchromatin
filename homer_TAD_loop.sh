# GMP G
merge2Dbed.pl ../analysis/TAD_loop/GMP.tad.2D.bed ../analysis/TAD_loop/G.tad.2D.bed -tad > ../analysis/TAD_loop/merge_GMP_G.tad.2D.bed
merge2Dbed.pl ../analysis/TAD_loop/GMP.loop.2D.bed ../analysis/TAD_loop/G.loop.2D.bed -loop > ../analysis/TAD_loop/merge_GMP_G.loop.2D.bed


findTADsAndLoops.pl score -tad ../analysis/TAD_loop/merge_GMP_G.tad.2D.bed -loop ../analysis/TAD_loop/merge_GMP_G.loop.2D.bed -o ../analysis/TAD_loop/Diff_GMP_G -d ../homer_merge/GMP ../homer_merge/G -cpu 10 -res 20000 -window 40000

# tad
merge2Dbed.pl ../analysis/TAD_loop/*.tad.2D.bed -tad -res 80000 > ../analysis/TAD_loop/merge.tad.2D.bed
findTADsAndLoops.pl score -tad ../analysis/TAD_loop/merge.tad.2D.bed -o ../analysis/TAD_loop/merge.tad.2D -d ../homer_merge/LT ../homer_merge/ST ../homer_merge/MPP ../homer_merge/CMP ../homer_merge/CLP ../homer_merge/MEP ../homer_merge/GMP ../homer_merge/MKP ../homer_merge/MK ../homer_merge/G -cpu 10 -res 20000 -window 40000


findTADsAndLoops.pl score -tad ../analysis/TAD_loop/merge.tad.2D.bed -o ../analysis/TAD_loop/merge.tad.2D.rep -d ../homer/CLP1  ../homer/G1    ../homer/LT1   ../homer/MK1   ../homer/MPP1 ../homer/CLP2  ../homer/G2    ../homer/LT2   ../homer/MK2   ../homer/MPP2 ../homer/CMP1  ../homer/GMP1  ../homer/MEP1  ../homer/MKP1  ../homer/ST1 ../homer/CMP2  ../homer/GMP2  ../homer/MEP2  ../homer/MKP2  ../homer/ST2 -cpu 15 


# loop

# loops from JUICER

merge2Dbed.pl ../analysis/loops/*.loop/merged_loops.bedpe  -loop -res 25000  > ../analysis/loops/merge.loops.bedpe
cut -f 1,2,3,4,5,6,11,12,13 ../analysis/loops/merge.loops.bedpe | awk 'BEGIN{OFS="\t"}{$1="chr"$1;$4="chr"$4;print}' > ../analysis/loops/merged.loops.bedpe

grep -v "#" ../analysis/loops/merge.loops.bedpe |  awk -F "\t" 'BEGIN{OFS="\t"}{if($1!=""){$1="chr"$1;$4="chr"$4;print}}' > ../analysis/loops/merged.loops.bedpe
~/bin/homer/bin/findTADsAndLoops.pl score -loop ../analysis/loops/merged.loops.bedpe -o ../analysis/loops/merged.loops -d ../homer/CLP1  ../homer/G1    ../homer/LT1   ../homer/MK1   ../homer/MPP1 ../homer/CLP2  ../homer/G2    ../homer/LT2   ../homer/MK2   ../homer/MPP2 ../homer/CMP1  ../homer/GMP1  ../homer/MEP1  ../homer/MKP1  ../homer/ST1 ../homer/CMP2  ../homer/GMP2  ../homer/MEP2  ../homer/MKP2  ../homer/ST2 -cpu 15 -res 10000 -window 30000


~/bin/homer/bin/findTADsAndLoops.pl score -loop ../analysis/loops/merged.loops.bedpe -o ../analysis/loops/merged.loops -d ../homer_merge/LT ../homer_merge/ST ../homer_merge/MPP ../homer_merge/CMP ../homer_merge/CLP ../homer_merge/MEP ../homer_merge/GMP ../homer_merge/MKP ../homer_merge/MK ../homer_merge/G -cpu 14 -res 10000 -window 30000

