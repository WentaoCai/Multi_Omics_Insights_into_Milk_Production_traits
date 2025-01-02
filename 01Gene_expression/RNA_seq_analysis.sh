########################## 1.1 Download sra_data ############################
for i in `cat SRA_list.txt`
do
prefetch "$i" -O ./
fastq-dump --split-3 --gzip "$i".sra
done


########################## 1.2 Quality control################################
#Pair-end quality control 
for i in `cat Paired.SRA_list.txt`
do
java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-phred33 "$i"_1.fastq.gz "$i"_2.fastq.gz "$i"_1.clean.fq.gz "$i"_1_unpaired.fastq.gz "$i"_2.clean.fq.gz "$i"_2_unpaired.fastq.gz \
-threads 10 ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

#Single-end quality control 
for i in `cat Single.SRA_list.txt`
do
java -jar ../Trimmomatic-0.38/trimmomatic-0.38.jar SE \
-phred33 "$i".fastq.gz "$i"_1.clean.fq.gz -threads 10 \
ILLUMINACLIP:../Trimmomatic-0.38/adapters/TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

########################## 1.3 Mapping Reads ###############################
#Build index
STAR --runMode genomeGenerate \
--genomeDir ./Bovine_genome \
--genomeFastaFiles ./Bovine_genome/Bosgenome.fa \
--sjdbGTFfile ./Bovine_genome/Bos_taurus.ARS-UCD1.2.102.chr.gtf \
--runThreadN 10 \

#Pair-end mapping
for i in `cat Paired.SRA_list.txt`
do
STAR --runThreadN 10 \
--genomeDir ./Bovine_genome \
--sjdbGTFfile ./Bovine_genome/Bos_taurus.ARS-UCD1.2.102.chr.gtf \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--readFilesCommand zcat \
--outFilterMismatchNmax 3 \
--readFilesIn "$i"_1.clean.fq.gz "$i"_1.clean.fq.gz \
--outFileNamePrefix "$i"-STAR
done

#Single-end mapping
for i in `cat Single.SRA_list.txt`
do
STAR --runThreadN 10 \
--genomeDir ./Bovine_genome \
--sjdbGTFfile ./Bovine_genome/Bos_taurus.ARS-UCD1.2.102.chr.gtf \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--readFilesCommand zcat \
--outFilterMismatchNmax 3 \
--readFilesIn "$i".clean.fq.gz \
--outFileNamePrefix "$i"-STAR
done

########################## 1.4 Gene expression quantification ###############################
# FPKM quantification
for i in `cat SRA_list.txt`
do
stringtie -p 10 -e -B -G ./Bovine_genome/Bos_taurus.ARS-UCD1.2.102.chr.gtf \
-o ./"$i"-STAR_STARgenome/"$i".gtf -A ./Expression/"$i".tsv "$i"-STARAligned.sortedByCoord.out.bam
done

cd Expression
for i in `cat ../SRA_list.txt`
do
awk '{print $1"\t"$8}' "$i".tsv|sed '1d'|sort -k1,1>"$i".fpkm
echo -e "ID\t"$i""|cat - "$i".fpkm >"$i".fpkm.txt
rm "$i".fpkm
done
paste *.fpkm.txt|awk '{printf("%s\t",$1);for(i=2;i<=NF;i+=2){printf("%s\t",$i)};print ""}'>All.gene.expression.fpkm.txt
cd ..

# Count quantification
for i in `cat SRA_list.txt`
do
featureCounts -T 10 -p -t exon -g gene_id -a ./Bovine_genome/Bos_taurus.ARS-UCD1.2.102.chr.gtf -o ./Expression/"$i".featureCounts.txt "$i"-STARAligned.sortedByCoord.out.bam
done

cd Expression
for i in `cat ../SRA_list.txt`
do
sed '/^#/d' "$i".featureCounts.txt|awk '{print $1"\t"$NF}'>"$i".count.txt
done
paste *.count.txt|awk '{printf("%s\t",$1);for(i=2;i<=NF;i+=2){printf("%s\t",$i)};print ""}'>All.gene.expression.count.txt