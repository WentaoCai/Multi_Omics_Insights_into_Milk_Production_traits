

###########################Step1 Building index###############################################
###Bowtie2 has prepared in default environment path.
bismark_genome_preparation --verbose /path/to/Bisulfite_Genome
#Note: genome fasta file should put in Bisulfite_Genome before excuting the above command.


###########################Step2 Mapping reads###############################################
for i in mammary_L mammary_N brain_L brain_N WBC_L WBC_N
do
#Mapping clean reads to genome using bowtie2
bismark --bowtie2 /path/to/Bisulfite_Genome/ -1 "$i"_R1.fq.gz -2 "$i"_R2.fq.gz -p 10 -N 1 -D 30
#Remove duplicates for mapped bam file
deduplicate_bismark -bam -p "$i"_R1_bismark_bt2_pe.bam --output_dir 02duplicate
#Get methylation level
bismark_methylation_extractor --comprehensive --gzip --no_overlap --remove_spaces -p --parallel 6 --bedGraph \
--counts --cytosine_report --report --buffer_size 30G --genome_folder /path/to/Bisulfite_Genome \
02duplicate/"$i"_R1_bismark_bt2_pe.deduplicated.bam -o ./03report >log 2>&1
done

###########################Step3 Identify differential methylation regions(DMRs)##############
###metilene should be prepared in 
for i in mammary brain WBC
do
#Prepare input file for metilene
perl merge_methylation_level.pl "$i"_L_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz "$i"_N_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz >"$i".common.methylation.level.txt
echo "chr\tpos\tg1_sample1\tg2_sample2"|cat - "$i".common.methylation.level.txt >"$i".common.methylation.level.withheader.txt
#Run metilene
metilene -a g1 -b g2 "$i".common.methylation.level.withheader.txt >"$i"_DMR.txt
#Filter significant DMR regions using threshold P-value<0.00001
awk '{if($6<=1e-5)print}' "$i"_DMR.txt >"$i"_DMR.1e5.bed
done