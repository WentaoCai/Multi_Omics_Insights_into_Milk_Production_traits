
####################>>>>>>>> SNP heritability using one genomic relationship matrix (GRM) model <<<<<<<<<<####################
#Build GRM using genotype of 23,566 Holstein bulls
gcta64 --autosome-num 33 --bfile genotype23566.bta9.1_29 --make-grm --out genotype23566.bta9.1_29 --thread-num 20

#Calculate the SNP heritability for 12 traits
for j in Milk Protein Fat Fat_per Pro_per SCS DCE DSB AFC Heifer_Conc_Rate Dtr_Preg_Rate Dtr_Preg_Rate Cow_Conc_Rate
do
gcta64 --reml --grm genotype23566.bta9.1_29 --autosome-num 29 --pheno "$j".gcta.pheno --reml-pred-rand --out "$j"_gcta --thread-num 20
done


####################>>>>>>>> Variance component analysis of 29 chromosomes using multiple GRMs model <<<<<<<<<<####################
#Split variants according to chromosome
for i in {1..29}
do
awk '{if($1=="'$i'") print $2}' genotype23566.bta9.1_29_QC.bim >chr"$i".SNPlist
done

#Make genomic relationship matrix 
for i in {1..29}
do
gcta64  --bfile genotype23566.bta9.1_29 --autosome-num 29 --extract chr"$i".SNPlist --make-grm chr"$i" --out chr"$i" --thread-num 10
echo 'chr'"$i"'' >>multi_chr.list
done

#Calculate the SNP heritability for 12 traits with multiple GRMs simultaneously
for j in Milk Protein Fat Fat_per Pro_per SCS DCE DSB AFC Heifer_Conc_Rate Dtr_Preg_Rate Dtr_Preg_Rate Cow_Conc_Rate
do
gcta64 --reml --mgrm multi_chr.list --autosome-num 29 --pheno "$j".gcta.pheno --reml-pred-rand --out "$j"_gcta_chr --thread-num 20
done


####################>>>>>>>>  Variance component analysis of genomic location using multiple GRMs model <<<<<<<<<<####################

### Variants were annotated into nine groups based on their genomic regions(Synonymous, Splicing, 3'-UTR, Missense, Downstream, 5'-UTR, Intron, Intergenic, Upstream) using VEP online tools(https://www.ensembl.org/vep)
#Make genomic relationship matrix for each functional class
for i in Synonymous Splicing 3UTR Missense Downstream 5UTR Intron Intergenic Upstream
do
gcta64 --bfile genotype23566.bta9.1_29 --autosome-num 29 --extract "$i".SNPlist --make-grm Region_"$i" --out Region_"$i" --thread-num 10
echo 'Region_'"$i"'' >>multi_location.list
done

#Calculate the SNP heritability for 12 traits with multiple GRMs simultaneously
for j in Milk Protein Fat Fat_per Pro_per SCS DCE DSB AFC Heifer_Conc_Rate Dtr_Preg_Rate Dtr_Preg_Rate Cow_Conc_Rate
do
gcta64 --reml --mgrm multi_location.list --autosome-num 29 --pheno "$j".gcta.pheno --reml-pred-rand --out "$j"_gcta_regions --thread-num 20
done


####################>>>>>>>> Variance component analysis of function class using multiple GRMs model <<<<<<<<<<####################

#Generate genotype position input file for bedtools
awk '{print $1"\t"$4"\t"$4"\t"$2}' genotype23566.bta9.1_29.bim >genotype23566.bta9.1_29.variants_pos.bed

#Variants were categorized into different groups based on their locations within functional regions using bedtools.
for i in class1 class2 ...
do
bedtools window -b "$i".gene.bed -a genotype23566.bta9.1_29.variants_pos.bed -w 5000|awk '{print $4}'|sort|uniq >"$i".SNPlist
#Note: -w set to 5000 for Specific gene, DEGs DE lncRNAs, DE miRNAs, -w set to 100000 for RNA editing, while -w set to 0 for DMRs and enhancers
#For eQTLs of sQTLs, skip the above command. 
#Make genomic relationship matrix
gcta64 --bfile genotype23566.bta9.1_29 --autosome-num 29 --extract "$i".SNPlist --make-grm Class_"$i" --out Class_"$i" --thread-num 10
echo 'Class_'"$i"'' >>multi_class.list
done

#Calculate variance for permutated variants with multiple GRMs simultaneously
for j in Milk Protein Fat Fat_per Pro_per SCS DCE DSB AFC Heifer_Conc_Rate Dtr_Preg_Rate Dtr_Preg_Rate Cow_Conc_Rate
do
gcta64 --reml --mgrm multi_class.list --autosome-num 29 --pheno "$j".gcta.pheno --reml-pred-rand --out "$j"_gcta_class --thread-num 20
done



####################>>>>>>>> Variance component analysis of permutated class using multiple GRMs model <<<<<<<<<<####################

for i in class1 class2 ...
do
#shuffle sites in genome with similar patterns within genotype.
perl extract_postion.pl "$i".SNPlist genotype23566.bta9.1_29.variants_pos.bed > "$i".SNPlist.pos.bed
#Permtation test for 10 times
for t in {1..10}
do
#Shuffle new position 
perl shuffle_with_genotype.pl genotype23566.bta9.1_29.bim "$i".SNPlist.pos.bed > "$i"."$t".SNPlist

gcta64 --bfile genotype23566.bta9.1_29 --autosome-num 29 --extract "$i"."$t".SNPlist --make-grm Class_"$i" --out Class_"$i"."$t" --thread-num 10
echo 'Class_'"$i"'.'"$t"'' >>multi_class."$t".list
done
done


#Permtation test for 10 times
for t in {1..10}
do
#Calculate variance for permutated variants with multiple GRMs simultaneously
for j in Milk Protein Fat Fat_per Pro_per SCS DCE DSB AFC Heifer_Conc_Rate Dtr_Preg_Rate Dtr_Preg_Rate Cow_Conc_Rate
do
gcta64 --reml --mgrm multi_class."$t".list --autosome-num 29 --pheno "$j".gcta.pheno --reml-pred-rand --out "$j"_gcta_class."$t" --thread-num 20
done
done
