########################## 1 Download sra_data ############################
for i in `cat miRNA_SRA_list.txt`
do
prefetch "$i" -O ./
fastq-dump --split-3 --gzip "$i".sra
done


########################## 2 Quality control################################
for i in `cat miRNA_SRA_list.txt`
do
java -jar ../Trimmomatic-0.38/trimmomatic-0.38.jar SE \
-phred33 "$i".fastq.gz "$i"_1.clean.fq.gz -threads 10 \
ILLUMINACLIP:../Trimmomatic-0.38/adapters/TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

########################## 3 Mapping reads ################################
#Please install MiRDeep2 software (contains collapse_reads_md.pl, mapper.pl, miRDeep2.pl)
for i in `cat miRNA_SRA_list.txt`
do
#merge reads 
collapse_reads_md.pl "$i"_1.rm_adaptor.fa bta>"$i".merge.fa
#mapping to genome
mapper.pl "$i".merge.fa -n -c -p ~/Genome_bowtie/Bosgenome -l 17 -t "$i".reads_collapsed_vs_genome.arf -s "$i".new.merged.fa >"$i".log
#Identify miRNAs
miRDeep2.pl "$i".new.merged.fa ~/Genome_bowtie/Bosgenome.fa "$i".reads_collapsed_vs_genome.arf bta.mature.fa other.mature.fa Bta.hairpined.fa -t Cow -c -d -v 2>"$i".report.log
quantifier.pl -p Bta.hairpined.fa -m bta.mature.fa -r "$i".new.merged.fa -y "$i"&
awk '{print $1"|"$3"\t"$2}' miRNAs_expressed_all_samples_"$i".csv|sed '1d'|sort -k1,1>"$i".count
echo -e "ID\t"$i""|cat - "$i".count >"$i".count.txt
rm "$i".count
done
paste *.count.txt|awk '{printf("%s\t",$1);for(i=2;i<=NF;i+=2){printf("%s\t",$i)};print ""}'>All.miRNA.expression.count.txt

######################### 4. miRNA target predictions######################
#Download cattle_gene_3UTR.fa from Ensembl BioMart website (https://www.ensembl.org/biomart/martview/). Choosing Sequences and 3'UTR in Attributes, then download and saved as cattle_gene_3UTR.fa.
miranda DE.miRNA.fa cattle_gene_3UTR.fa >Tagets_results.txt
