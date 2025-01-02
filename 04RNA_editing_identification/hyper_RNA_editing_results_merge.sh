#!/usr/bin/bash
for i in `cat $1`;
do
echo "Fiter Reads for $i"
awk '{if($NF<=($(NF-1)+$NF)*0.2 && $(NF-1)>=3)print}' "$i"_1.UE.list> "$i"_1_SE.UE.filter.list
awk '{if($NF<=($(NF-1)+$NF)*0.2 && $(NF-1)>=3)print}' "$i"_2.UE.list> "$i"_2_SE.UE.filter.list
#wc -l "$i"_1_SE.UE.filter.list>> Hyper_editing.reads.number.txt
#wc -l "$i"_2_SE.UE.filter.list>> Hyper_editing.reads.number.txt
grep 'Aligns to:' "$i"_1.UE.Details|awk  '{print $NF"\t"$(NF-3)"\t"$(NF-6)}'|sed 's/,//g'|awk -F '[()]' '{print $0"\t"$2}'|awk '{if($4=="-")$4="0";else if($4="+")$4="1"; print}' > "$i"_1.position.txt
grep 'Edit indexes:' "$i"_1.UE.Details|awk '{print $3}'|sed 's/,//g'> "$i"_1.ps.txt
paste "$i"_1.position.txt "$i"_1.ps.txt> "$i"_1.pos.txt
perl extract_Reads.pl "$i"_1_SE.UE.filter.list "$i"_1.pos.txt >"$i"_1.cluster.txt
perl Editor_Position.pl "$i"_1.cluster.txt> "$i"_1.site.txt
rm "$i"_1.pos.txt
rm "$i"_1.ps.txt
rm "$i"_1.position.txt
grep 'Aligns to:' "$i"_2.UE.Details|awk  '{print $NF"\t"$(NF-3)"\t"$(NF-6)}'|sed 's/,//g'|awk -F '[()]' '{print $0"\t"$2}'|awk '{if($4=="-")$4="0";else if($4="+")$4="1"; print}' >"$i"_2.position.txt
grep 'Edit indexes:' "$i"_2.UE.Details|awk '{print $3}'|sed 's/,//g'> "$i"_2.ps.txt
paste "$i"_2.position.txt "$i"_2.ps.txt> "$i"_2.pos.txt
perl extract_Reads.pl "$i"_2_SE.UE.filter.list "$i"_2.pos.txt >"$i"_2.cluster.txt
perl Editor_Position.pl "$i"_2.cluster.txt> "$i"_2.site.txt
rm "$i"_2.pos.txt
rm "$i"_2.ps.txt
rm "$i"_2.position.txt
###merge two reads
echo "Merge Pair Reads for $i"
awk '{if($4==0&&$3=="A2G")$3="T2C";else if($4==0&&$3=="A2T")$3="T2A";else if($4==0&&$3=="A2C")$3="T2G";else if($4==0&&$3=="T2G")$3="A2C";else if($4==0&&$3=="T2C")$3="A2G"; else if($4==0&&$3=="T2A")$3="A2T";else if($4==0&&$3=="C2A")$3="G2T";else if($4==0&&$3=="C2G")$3="G2C";else if($4==0&&$3=="C2T")$3="G2A";else if($4==0&&$3=="G2A")$3="C2T";else if($4==0&&$3=="G2T")$3="C2A";else if($4==0&&$3=="G2C")$3="C2G"; print}' "$i"_2.site.txt> "$i"_2.position
awk '{if($4==1&&$3=="A2G")$3="T2C";else if($4==1&&$3=="A2T")$3="T2A";else if($4==1&&$3=="A2C")$3="T2G";else if($4==1&&$3=="T2G")$3="A2C";else if($4==1&&$3=="T2C")$3="A2G"; else if($4==0&&$3=="T2A")$3="A2T";else if($4==1&&$3=="C2A")$3="G2T";else if($4==1&&$3=="C2G")$3="G2C";else if($4==1&&$3=="C2T")$3="G2A";else if($4==1&&$3=="G2A")$3="C2T";else if($4==1&&$3=="G2T")$3="C2A";else if($4==1&&$3=="G2C")$3="C2G"; print}' "$i"_1.site.txt> "$i"_1.position
awk '{if($3==0&&$2=="A2G")$2="T2C";else if($3==0&&$2=="A2T")$2="T2A";else if($3==0&&$2=="A2C")$2="T2G";else if($3==0&&$2=="T2G")$2="A2C";else if($3==0&&$2=="T2C")$2="A2G"; else if($3==0&&$2=="T2A")$2="A2T";else if($3==0&&$2=="C2A")$2="G2T";else if($3==0&&$2=="C2G")$2="G2C";else if($3==0&&$2=="C2T")$2="G2A";else if($3==0&&$2=="G2A")$2="C2T";else if($3==0&&$2=="G2T")$2="C2A";else if($3==0&&$2=="G2C")$2="C2G"; print}' "$i"_2.cluster.txt> "$i"_2.reads.position
awk '{if($3==1&&$2=="A2G")$2="T2C";else if($3==1&&$2=="A2T")$2="T2A";else if($3==1&&$2=="A2C")$2="T2G";else if($3==1&&$2=="T2G")$2="A2C";else if($3==1&&$2=="T2C")$2="A2G"; else if($3==0&&$2=="T2A")$2="A2T";else if($3==1&&$2=="C2A")$2="G2T";else if($3==1&&$2=="C2G")$2="G2C";else if($3==1&&$2=="C2T")$2="G2A";else if($3==1&&$2=="G2A")$2="C2T";else if($3==1&&$2=="G2T")$2="C2A";else if($3==1&&$2=="G2C")$2="C2G"; print}' "$i"_1.cluster.txt> "$i"_1.reads.position
cat "$i"_2.position "$i"_1.position>"$i".event.site
awk '{print $1"\t"$2"\t"$3}' "$i".event.site|sort -k1,1 -k2,2n|uniq>"$i".uniq.position
cat "$i"_1.reads.position "$i"_2.reads.position>"$i".reads.position
rm "$i"_1.reads.position
rm "$i"_2.reads.position
rm "$i"_2.position "$i"_1.position
#rm "$i"_1.uniq.position
for j in A2G A2C A2T C2G C2A C2T G2C G2A G2T T2A T2G T2C
do
perl Cluster150.pl "$i".uniq.position "$j">> "$i".cluster.tem.position
done
grep -v '^$' "$i".cluster.tem.position>"$i".cluster.position
rm "$i".cluster.tem.position
awk -F '[\t,]' '{print $1"\t"$3"\t"$(NF-1)"\t"$2}' "$i".cluster.position > "$i".cluster.bed
done
#####
for j in A2G A2C A2T C2G C2A C2T G2C G2A G2T T2A T2G T2C
do
echo "Merge cluster for $j"
for i in `cat $1`
do
grep ''"$j"'' "$i".cluster.bed|awk '{print $1"\t"$2"\t"$3"\t"'"$i"'}' > "$i"."$j".bed
done
cat *."$j".bed|sort -k1,1 -k2,2n > All."$j".bed
bedtools merge -i All."$j".bed -c 4 -o collapse >All.merge."$j".bed
rm *[0-9]."$j".bed
rm All."$j".bed
done
#### Count Cluster type
echo "Count type"
echo -e -n "Type" >>Cluster_type.txt
echo -e -n "Type" >>Uniq_sites_type.txt
echo -e -n "Type" >>Editing_reads_type.txt
echo -e -n "Type" >>Event_reads_type.txt
for j in `cat $1`
do
echo -e -n "\t$j"  >>Cluster_type.txt
echo -e -n "\t$j"  >>Uniq_sites_type.txt
echo -e -n "\t$j"  >>Editing_reads_type.txt
echo -e -n "\t$j"  >>Event_reads_type.txt
done
echo ""  >>Cluster_type.txt
echo ""  >>Uniq_sites_type.txt
echo ""  >>Editing_reads_type.txt
echo ""  >>Event_reads_type.txt
for i in A2G A2C A2T C2G C2A C2T G2C G2A G2T T2A T2G T2C
do
echo -n "$i"  >>Cluster_type.txt
echo -n "$i"  >>Uniq_sites_type.txt
echo -n "$i"  >>Editing_reads_type.txt
echo -n "$i"  >>Event_reads_type.txt
for j in `cat $1`
do
echo -e -n "\t"  >>Cluster_type.txt
echo -e -n "\t"  >>Uniq_sites_type.txt
echo -e -n "\t"  >>Editing_reads_type.txt
echo -e -n "\t" >>Event_reads_type.txt
grep ''"$i"'' "$j".cluster.position|wc -l|tr -d '\n'  >>Cluster_type.txt
grep ''"$i"'' "$j".uniq.position|wc -l|tr -d '\n'  >>Uniq_sites_type.txt
grep ''"$i"'' "$j".reads.position|wc -l|tr -d '\n'  >>Editing_reads_type.txt
grep ''"$i"'' "$j".event.site|wc -l|tr -d '\n'  >>Event_reads_type.txt
done
echo "" >>Cluster_type.txt
echo "" >>Uniq_sites_type.txt
echo "" >>Editing_reads_type.txt
echo "" >>Event_reads_type.txt
done

echo "All Finish...Enjoy!"

