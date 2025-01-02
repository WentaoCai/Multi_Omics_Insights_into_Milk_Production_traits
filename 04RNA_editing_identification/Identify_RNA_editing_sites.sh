
####################  1. Build BWA index for six converted genome fasta ########################
#Separately transformed all A to G, A to C, A to T, G to C, T to C, T to G in genome. Building the index for these six converted genomes, separately.
perl BWA_genome_index.pl /genome_path/Bosgenome.fa

####################    2. Mapping reads using   ########################
#Separately transformed all A to G, A to C, A to T, G to C, T to C, T to G in genome. Building the index for these six converted genomes, separately.

#Mapping Reads1
for i in `cat SRA_list.txt`
do
STAR --runThreadN 10 \
--genomeDir ./Bovine_genome \
--sjdbGTFfile ./Bovine_genome/Bos_taurus.ARS-UCD1.2.102.chr.gtf \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--readFilesCommand zcat \
--outFilterMismatchNmax 3 \
--readFilesIn "$i"_1.clean.fq.gz \
--outFileNamePrefix "$i"_1-STAR
done

#Mapping Reads2
for i in `cat SRA_list.txt`
do
STAR --runThreadN 10 \
--genomeDir ./Bovine_genome \
--sjdbGTFfile ./Bovine_genome/Bos_taurus.ARS-UCD1.2.102.chr.gtf \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--readFilesCommand zcat \
--outFilterMismatchNmax 3 \
--readFilesIn "$i"_2.clean.fq.gz \
--outFileNamePrefix "$i"_2-STAR
done


################ 3. Identify hyper RNA editing sites from unmapped reads.  ########################
#First extracted all unmapped reads from the initial alignment. Then transformed all As to Gs in the unmapped reads and realigned the transformed RNA reads to the transformed A to G reference genome using BWA. 
#The process of handling unmapped reads is detailed in our script 'hyper_editing_for_bwa.pl'. The unique simple repeats were downloaded from UCSC website.
perl hyper_editing_for_bwa.pl "$i"_1-STARAligned.sortedByCoord.out.bam /path/to/index/Bosgenome /path/to/fasta/Bosgenome.fa unique_simple_repeats.txt "$i"_1 &
perl hyper_editing_for_bwa.pl "$i"_2-STARAligned.sortedByCoord.out.bam /path/to/index/Bosgenome /path/to/fasta/Bosgenome.fa unique_simple_repeats.txt "$i"_2 &

#Merge RNA editing results fA2Grom two reads for each sample
sh RNA_editing_merge.sh SRA_list.txt
get
cat *uniq.position|grep -E 'A2G|T2C'|awk '{print $1"_"$2}' >All.A2G.editing.id

for i in `cat SRA_list.txt`
do
#Get editing reads count of each hyper editing site for each sample
awk '{print $1"_"$2}' "$i".event.site|sort|uniq -c|sort -k1,1nr|awk '{print $2"\t"$1}'>"$i".hyper_editing.count.txt
#Get A2G editing sites for each sample
grep -E 'A2G|T2C' "$i".uniq.position|awk '{print $1"_"$2}' >"$i".A2G.editing.id
done
#Get all hyper editing sites across all mammary tissues
cat *.A2G.editing.id|sort|uniq> All.A2G.editing.id

################ 4. Calculate RNA editing levels.  ########################

for i in `cat SRA_list.txt`
do
#improve bam quality
java -jar -Xmx50g ~/bin/picard.jar MarkDuplicates I="$i"_STARAligned.sortedByCoord.out.bam M="$i"_dedup_metrics.txt O="$i"_sorted_dedup_reads.bam
java -jar -Xmx50g ~/bin/picard.jar AddOrReplaceReadGroups I="$i"_sorted_dedup_reads.bam O="$i".prepared.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="$i"
gatk --java-options -Xmx50G SplitNCigarReads -R /path/to/fasta/Bosgenome.fa -I "$i".prepared.bam -O  "$i".SplitNCigarReads.bam
#Get mutation information from mapped bam file 
python /path/to/reditools2.0/src/cineca/reditools.py -f "$i".SplitNCigarReads.bam -r /path/to/fasta/Bosgenome.fa -o "$i".table
#Extract the mapped editing sites and levels 
perl get_mapped_editing_count.pl "$i".table All.A2G.editing.id >"$i".mapped_count.txt
#Recalculate the editing level from both mapped editing sites and unmapped hyper editing sites
perl get_editing_level.pl "$i".hyper_editing.count.txt "$i".mapped_count.txt>"$i".editing_level.txt
done
