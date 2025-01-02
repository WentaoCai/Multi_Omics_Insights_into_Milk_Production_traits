
#######################1. Basic filtering of transcripts
#Annotate lncRNA transcrits using class code 
gffcompare -r Bos_taurus.ARS-UCD1.2.102.chr.gtf -G -o merged_from_list -i mergelist.txt
#Extract class code 'i' and 'u', which represented novel intronic and intergenic transcripts
grep 'class_code "i"' merged_from_list.combined.gtf|awk -F "\"" '{print $2}'>i.transcript.id
grep 'class_code "u"' merged_from_list.combined.gtf|awk -F "\"" '{print $2}'>u.transcript.id
cat u.transcript.id i.transcript.id>iu.transcript.id
awk '{for(i=5;i<=NF;i++) if($i=="-") sum++} {if(sum<5)print $1} {sum=0}'  merged_from_list.tracking.txt >identified_at_least_five.transcript.id
cat iu.transcript.id identified_at_least_five.transcript.id|sort|uniq -d >iu.identified_at_least_five.transcript.id
perl extract_gtf.pl iu.identified_at_least_five.transcript.id merged_from_list.combined.gtf>iu.dentified_at_least_five.transcript.gtf
#obtaining trabnscripts with exon >=2 and length>=200
perl summary_gtf.pl -i iu.dentified_at_least_five.transcript.gtf -o length_exon.txt
perl extract_gtf.pl length_exon.txt iu.gtf> length_exon.gtf
gffread -w length_exon.fa -g Bosgenome.fa length_exon.gtf

#######################Step 2. ORF filtering
wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.5.0.tar.gz
tar -xzvf TransDecoder-v5.5.0.tar.gz
./TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t length_exon.fa
grep '>' longest_orfs.pep|awk '{print $3"\t"$5}'|awk -F '[:\t]' '{if($2<121)print $3}'|sort|uniq > ../orf.id
perl extract_gtf.pl orf.id length_exon.gtf > orf.gtf
gffread -w orf.fa -g Bosgenome.fa orf.gtf

#######################Step 3. Coding protential 
#Install CPC, CNCI, PLEK before the following command
#CPC
python CPC2-beta/bin/CPC2.py -i orf.fa -o CPC_results
grep 'noncoding' CPC_results.txt|awk '{print $1}'>CPC.noncoding.id

#CNCI
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
chmod +x faToTwoBit
./faToTwoBit Bosgenome.fa Bosgenome.2bit
python CNCI.py -f ../orf.gtf -g -o CNCI_results -m ve -p 16 -d ../Bosgenome.2bit
grep 'noncoding' ./CNCI/CNCI_results/CNCI.index|awk '{print $1}'>CNCI.noncoding.id

#PLEK
python ./PLEK.1.2/PLEK.py -fasta orf.fa -out PLEK_output -thread 10
grep 'Non-coding' PLEK_output|awk '{print $3}'|sed 's/^>//g' >PLEK.noncoding.id

# The ovarlapped results across three methods
cat PLEK.noncoding.id CPC.noncoding.id CNCI.noncoding.id|sort|uniq -c|awk '{if($1==3)print $2}'>lncRNA.common.id
perl extract_gtf.pl lncRNA.common.id orf.gtf > lncRNA.common.gtf
gffread -w lncRNA.common.fa -g Bosgenome.fa lncRNA.common.gtf

#Pfam 
#Translate candidate lncRNA sequence to six types of protein sequence. 
#Forward translation had three types start from 0,1,2
perl DNA2protein.pl lncRNA.common.fa 0 lncRNA.pro.f0.fa
perl DNA2protein.pl lncRNA.common.fa 1 lncRNA.pro.f1.fa
perl DNA2protein.pl lncRNA.common.fa 2 lncRNA.pro.f2.fa

#Reverse candidate lncRNA sequence
perl reverse_sequence.pl lncRNA.common.fa lncRNA.rev.fa

#Reverse translation had three types start from 0,1,2
perl DNA2protein.pl lncRNA.rev.fa 0 lncRNA.pro.r0.fa
perl DNA2protein.pl lncRNA.rev.fa 1 lncRNA.pro.r1.fa
perl DNA2protein.pl lncRNA.rev.fa 2 lncRNA.pro.r2.fa

#Install PfamScan and hmmer3 before the following command
#Download protein sequencing from EBI database(https://www.ebi.ac.uk/interpro/download/pfam/)
mkdir Pfam_data
cd Pfam_data
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
gunzip Pfam-A.hmm.gz Pfam-A.hmm.dat.gz
#Build library
hmmpress Pfam-A.hmm
cd ..
#lncRNA tanslated peptide sequences mapped to pfam library
pfam_scan.pl lncRNA.pro.f0.fa -fasta -dir Pfam_data -out Results.f0.txt -cpu 8
pfam_scan.pl lncRNA.pro.f1.fa -fasta -dir Pfam_data -out Results.f1.txt -cpu 8
pfam_scan.pl lncRNA.pro.f2.fa -fasta -dir Pfam_data -out Results.f2.txt -cpu 8
pfam_scan.pl lncRNA.pro.r0.fa -fasta -dir Pfam_data -out Results.r0.txt -cpu 8
pfam_scan.pl lncRNA.pro.r1.fa -fasta -dir Pfam_data -out Results.r1.txt -cpu 8
pfam_scan.pl lncRNA.pro.r2.fa -fasta -dir Pfam_data -out Results.r2.txt -cpu 8

#Extract and remove the Pfam mapped lncRNA transcripts
cat Results.f*.txt Results.r*.txt |grep -v "^#"|awk '{print $1}'|sort|uniq > pfam.coding.id
cat pfam.coding.id lncRNA.common.id|sort|uniq -u > lncRNA.final.id
perl extract_gtf.pl lncRNA.final.id lncRNA.common.gtf> lncRNA.final.gtf
