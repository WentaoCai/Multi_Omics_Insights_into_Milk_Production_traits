#####################################Step1 Mappping reads using bwa #################################
#Three cattle
for i in 2181 6819 Daisy 
do
#Six tissues
for j in MG Lung Heart Liver Spleen Kidney
do
#Five types of histone modification
for k in H3K4Me3_input H3K4Me1_input H3K27ac_input H3K27Me3_input CTCF_input H3K4Me3_ChIP H3K4Me1_ChIP H3K27ac_ChIP H3K27Me3_ChIP CTCF_ChIP
do
bwa mem -t 10 /path/to/Bosgenome "$i"_1.clean.fq.gz "$i"_2.clean.fq.gz|samtools sort >"$i".bam
samtools view -h "$i".bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > "$i".filtered.bam
samtools rmdup "$i".filtered.bam "$i".rmdup.bam
done


#####################################Step2 Call peaks using MACS2 #################################
#Three cattle
for i in 2181 6819 Daisy 
do
#Six tissues
for j in MG Lung Heart Liver Spleen Kidney
do
#Five types of histone modification
for k in H3K4Me3 H3K4Me1 H3K27ac H3K27Me3 CTCF
do
macs2 callpeak -c "$i"_"$j"_"$k"_input.rmdup.bam -t "$i"_"$j"_"$k"_ChIP.rmdup.bam -q 0.05 -f BAM -g 2.5e9 -n "$i"_"$j"_"$k".macs2 2>"$i"_"$j"_"$k".macs2.log&
done
done
done

#####################################Step3 Merge Peaks using ChIP-R #################################
#Six tissues
for j in MG Lung Heart Liver Spleen Kidney
do
#Five types of histone modification
for k in H3K4Me3 H3K4Me1 H3K27ac H3K27Me3 CTCF
do
chipr -i 2181_"$j"_"$k".macs2_peaks.narrowPeak 6819_"$j"_"$k".macs2_peaks.narrowPeak Daisy_"$j"_"$k".macs2_peaks.narrowPeak -m 2 -o "$j"-"$k"
done
done
