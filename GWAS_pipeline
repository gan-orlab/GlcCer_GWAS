# GWAS pipeline for GlcCer manuscript
### includes outliers and uses softcalls with maf > 0.05

cat ../../phenos/phenos.txt | while read line #Loops through all isoforms
do
## Linear regression
### Have to manually include N370S due to low maf
plink --bfile ../../ppmi_hg38_noindels --pheno ../../phenos/$line.txt --update-sex sex.txt --make-bed --out ppmi_sc_tm_$line
plink --bfile ppmi_sc_tm_$line --linear --ci .95 --maf 0.05 --covar pdhc_covar_2.txt --covar-name Sex,Age,Status,wbc,PC1,PC2,PC3,PC4,PC5 --out ppmi_sc_tm_maf05_$line
awk '{if ($12!="NA") print $0}' ppmi_sc_tm_maf05_$line.assoc.linear | grep 'ADD' | sort -gk12 > p.ppmi_sc_tm_maf05_$line.assoc.linear
plink --bfile ../../ppmi_hg38_noindels --pheno ../../phenos/$line.txt --update-sex sex.txt --extract GBA_N370S.txt --make-bed --out ppmi_sc_tm_GBA_$line
plink --bfile ppmi_sc_tm_GBA_$line --linear --ci .95 --covar pdhc_covar_2.txt --covar-name Sex,Age,Status,wbc,PC1,PC2,PC3,PC4,PC5 --out ppmi_sc_tm_maf05_GBA_$line
awk '{if ($12!="NA") print $0}' ppmi_sc_tm_maf05_GBA_$line.assoc.linear | grep 'ADD' | sort -gk12 > p.ppmi_sc_tm_maf05_GBA_$line.assoc.linear
cat p.ppmi_sc_tm_maf05_GBA_$line.assoc.linear p.ppmi_sc_tm_maf05_$line.assoc.linear > linear/p.ppmi_sc_tm_maf05_$line.assoc.linear
cat ppmi_sc_tm_maf05_GBA_$line.assoc.linear ppmi_sc_tm_maf05_$line.assoc.linear > linear/ppmi_sc_tm_maf05_$line.assoc.linear
mv *.log log/
rm ppmi_sc_tm*
rm p.*

## Generate plots and summary statistics
echo $line > marker.txt
cp linear/p.ppmi_sc_tm_maf05_$line.assoc.linear assoc
R < ~/runs/emsom/scripts/sumstats_from_assoc.R --no-save
paste lambda.txt >> plots/lambdas_all.txt
mv QQ.tiff plots/QQ_ppmi_sc_tm_maf05_$line.tiff
mv Metal.tab sumstats/Metal_ppmi_sc_tm_maf05_$line.txt
mv COJO.tab sumstats/COJO_ppmi_sc_tm_maf05_$line.txt
mv fullSTATS.tab sumstats/fullSTATS_ppmi_sc_tm_maf05_$line.txt
mv ManH.tiff plots/ManH_ppmi_sc_tm_maf05_$line.tiff

## Annotate summary statistics
head -n1 sumstats/fullSTATS_ppmi_sc_tm_maf05_$line.txt > header.txt
awk '{print $1,$2,$2,$5,$4}' sumstats/fullSTATS_ppmi_sc_tm_maf05_$line.txt | sed '1d' > input.annovar.txt
perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar.txt ~/runs/emsom/softwares/annovar/humandb/ \
    -buildver hg38 -out annotatedtable -remove -protocol refGene,avsnp150,clinvar_20210501,gnomad_genome \
    -operation g,f,f,f -nastring . -polish
paste sumstats/fullSTATS_ppmi_sc_tm_maf05_$line.txt annotatedtable.hg38_multianno.txt > temp.txt
tr ' ' '\t' < temp.txt > temp2.txt
sort -gk10 temp2.txt > annos/anno_ppmi_sc_tm_maf05_$line.txt

rm ppmi_sc_tm_$line
rm shapiro.txt
rm marker.txt
rm lambda.txt
rm header.txt
rm temp.txt
rm input.annovar.txt
rm annotatedtable.hg38_multianno.txt
done











