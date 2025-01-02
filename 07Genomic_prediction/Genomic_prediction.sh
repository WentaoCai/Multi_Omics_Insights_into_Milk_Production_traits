############# Build GRM for all genotype variants ###########################
ldak5.2.linux --calc-kins-direct Kinship_all --bfile Reference19575.bta9.QC --ignore-weights YES --power -1 --max-threads 20

############## Run GBLUP ###########################
for i in Milk Protein Fat Fat_per Pro_per
do
ldak5.2.linux --reml "$i"_all --pheno Ref."$i".phen --grm Kinship_all --max-threads 50
ldak5.2.linux --calc-blups "$i"_all --remlfile "$i"_all.reml --grm Kinship_all --bfile Reference19575.bta9.QC
##Predict the phenotype of validation population based on their genotype 
plink --cow --bfile Validation3991.bta9.QC --out "$i"_all --score "$i"_all.blup.full 1 2 5 sum
##The phenotype were located in the third column of Val.phe.csv, the predictive reliability were calculated using Pearson correlation
Rscript GBLUP_correlation.R Val.phe.csv 3 "$i"_all.profile >"$i"_all.reliability.txt
done

############# Build GRM for variants of two class ###########################
# Variants located in functional class 
plink --chr-set 29 --bfile Reference19575.bta9.QC --extract functional.list --make-bed --out Functional
ldak5.2.linux --calc-kins-direct Functional --bfile Functional --ignore-weights YES --power -1 --max-threads 20
echo 'Functional' >twocomp.list
plink --chr-set 29 --bfile Reference19575.bta9.QC --exclude functional.list --make-bed --out other
ldak5.2.linux --calc-kins-direct Other --bfile Other --ignore-weights YES --power -1 --max-threads 20
echo 'Other' >>twocomp.list


############## Run MultiGBLUP ###########################
for i in Milk Protein Fat Fat_per Pro_per
do
ldak5.2.linux --reml "$i"_2comp --pheno Ref."$i".phen --mgrm twocomp.list --bfile Reference19575.bta9.QC --ignore-weights YES --power -1 --max-threads 20
ldak5.2.linux --calc-blups "$i"_2comp --mgrm merge.list --remlfile "$i"_2comp.reml --bfile ../../Reference19575.bta9.QC --max-threads 50
##The phenotype were located in the third column of Val.phe.csv, the predictive reliability were calculated using Pearson correlation
plink --cow --bfile Validation3991.bta9.QC --out "$i"_2comp --score "$i"_2comp.blup 1 2 5 sum
##The phenotype were located in the third column of Val.phe.csv
Rscript GBLUP_correlation.R Val.phe.csv 3 "$i"_2comp.profile >"$i"_2comp.reliability.txt
done

############## Run BayesR #####################################
#Build all.annote.txt with one column of value 1 for BayesR
awk '{print "1"}' ../Reference19575.bta9.QC.bim > all.annote.txt
for i in {1..5}
do
#The phenotype of milk production traits were put in the 6~10th column of Reference19575.bta9.QC.fam
bayesRCO -bfile Reference19575.bta9.QC -n "$i" -ncat 1 -catfile all.annote.txt -out BayesR."$i"  -burnin 5000 -numit 20000
bayesRCO -bfile Validation3991.bta9.QC -predict -out BayesR."$i" -model BayesR."$i".model -freq BayesR."$i".frq -param BayesR."$i".param -ncat 1 -catfile all.annote.txt
##The phenotype were located in the third column of Val.phe.csv
Rscript BayesR_correlation.R Val.phe.csv 3 BayesR."$i".gv >BayesR."$i".reliability.txt
done
############## Run BayesRC #####################################
#2comp.annote.txt with two columns of value 0 or 1 based on function annotation for BayesRC
for i in {1..5}
do
#The phenotype of milk production traits were put in the 6~10th column of Reference19575.bta9.QC.fam
bayesRCO -bfile Reference19575.bta9.QC -n "$i" -ncat 2 -catfile 2comp.annote.txt  -out BayesRC."$i" -burnin 5000 -numit 20000
bayesRCO -bfile Validation3991.bta9.QC -predict -out BayesRC."$i" -model BayesRC."$i".model -freq BayesRC."$i".frq -param BayesRC."$i".param -ncat 2 -catfile 2comp.annote.txt
##The phenotype were located in the third column of Val.phe.csv
Rscript BayesR_correlation.R Val.phe.csv 3 BayesRC."$i".gv >BayesRC."$i".reliability.txt
done
