####################7/2/20###################
###########DADES GCAT#########################
#Preparaci? de files pgen pvar psam
#plink
#we first install the needed packages
install.packages("caTools")
library(tidyr)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(reshape)
library(numpy)

#go to the directory in which we have the data
setwd("I:/RDecid/Ana_Sanchez/questionari_taules")

#Read the file
tot <- fread("conditions_nivell1.csv" ,na.strings = "NULL")
#We keep the "participants" and discard the "familiars"
tot2 <- filter(tot, FAMILIAR=="PARTICIPANT")
#Take off the repeated "entity_id"
tot3 <- select(tot2,entity_id, TIPO) %>% unique()
#We create the file  tipos_malalties that consists in two columns: One with the diseases (CIE codes) and the other column has 
#the number of persons from our set of "participants that suffer it
tipos_malalties <- data.frame(tot3$TIPO) %>% group_by(tot3$TIPO) %>% summarise(N = n()) %>% `colnames<-`(c('TIPO', 'N'))

#Filternby those diseases that have more than 40 participants
#Create a file with those (55 different diseases)
disease_2pmil <- (filter(tipos_malalties, N > 40 & !is.na(TIPO)))

#Create the file with the filtered diseases and its entity_id
tot2_filtre <- filter(tot3, TIPO %in% disease_2pmil$TIPO)

#Create a dataframe of a matrix of 0s and 1s. 1s for the SNPs that have the disease
matriu<- data.frame(tot2_filtre, model.matrix(~TIPO+0,tot2_filtre))

#Read "consulta" in order to also take into account the healthy participants, that will be the control variants.
consulta <- fread("consulta_nivell2.csv", na.strings = "NULL")

#Take the healthy individuals
sans = filter(consulta, !entity_id %in% matriu$entity_id) #no malalties

#Create a matrix of 0s with the healthy individuals
matriu_sans = data.frame(matrix(data = c(sans$entity_id, rep(0,56*6777)), ncol = 57, nrow = 6777)) %>% 
  `colnames<-`(colnames(matriu))

#We join the two matrices
matriu_all <- rbind(matriu, matriu_sans) %>% select(-TIPO)

matriu_all <- cbind(matriu_all[1], data.frame(sapply(matriu_all[-1], as.numeric))) %>% group_by(entity_id) %>% summarise_all(sum)

matriu_all_final <- merge(matriu_all, consulta %>% select(entity_id, GENOTYPED_SAMPLE), by = 'entity_id')%>% filter(!is.na(GENOTYPED_SAMPLE)) %>% select(-entity_id)

#Verify that our  matrix is OK
sum(matriu_all_final$TIPO558)
colSums(matriu_all_final)

#Create a new column on the matrix with the diseases IDs
matriu1 = as.matrix(matriu_all_final[1:55]+1)
IID = matriu_all_final[56]
#Put the column as first column
matriu_ultima <- cbind(IID, matriu1)

#List of the CIE 9 codes of the 55 diseases we are goiong to analyse
#TIPO041,TIPO242,TIPO244,TIPO246,TIPO250,TIPO268,TIPO272,TIPO274,TIPO275,TIPO280,TIPO285,TIPO296,TIPO300,TIPO311,TIPO346,TIPO365,TIPO401,TIPO427,TIPO432,TIPO459,TIPO472,TIPO493,TIPO496,TIPO530,TIPO536,TIPO537,TIPO553,TIPO558,TIPO564,TIPO569,TIPO600,TIPO602,TIPO626,TIPO627,TIPO692,TIPO696,TIPO704,TIPO710,TIPO715,TIPO716,TIPO719,TIPO722,TIPO724,TIPO728,TIPO729,TIPO733,TIPO780,TIPO784,TIPO788,TIPO799,TIPO995,TIPOC43,TIPOC44,TIPOC50,TIPOV25
table(matriu_ultima$TIPO)
lista_tipo <- c("TIPO041","TIPO242","TIPO244","TIPO246","TIPO250","TIPO268","TIPO272","TIPO274","TIPO275","TIPO280","TIPO285","TIPO296","TIPO300","TIPO311","TIPO346","TIPO365","TIPO401","TIPO427","TIPO432","TIPO459","TIPO472","TIPO493","TIPO496","TIPO530","TIPO536","TIPO537","TIPO553","TIPO558","TIPO564","TIPO569","TIPO600","TIPO602","TIPO626","TIPO627","TIPO692","TIPO696","TIPO704","TIPO710","TIPO715","TIPO716","TIPO719","TIPO722","TIPO724","TIPO728","TIPO729","TIPO733","TIPO780","TIPO784","TIPO788","TIPO799","TIPO995","TIPOC43","TIPOC44","TIPOC50","TIPOV25")

#We add to the IID column the disease name
for(TIPO in lista_tipo){
  l=""
  i <- 1
  while(i< 55) {
    l[[i]] <- sum(matriu_all_final$TIPO)
    i <- i + 1
  }
}

#Safe the matrix
write.table(matriu_ultima, "I:/RDecid/Ana_Sanchez/GWAS/matriu_ultima.txt", sep = " ", row.names= F, quote = F)



############################
#Terminal: 
#Join the psam file with "matriu final"

#Change the name "entity_id" by "IID" in order to have the same name for the columns
require(reshape)
matriu_all <- fread("matriu_ultima.txt", na.strings = "NULL")
matriu_all_IID <- rename(matriu_all, c("GENOTYPED_SAMPLE"="IID"))
write.table(matriu_all_IID, "matriu_all_IDD.txt", sep = " ", row.names= F, quote = F)

#Takes off the rows of the psam  we don't want.
cut -d ' ' -f 1,2,3,4,5,6,7,8,9 GCAT_imputed.psam > GCAT_imputed_nochol.psam

#Do the merge of "matriu_final" and the .psam file by the "IID" and keep the needed columns.
GCAT_imputed <- fread("GCAT_imputed.psam", na.strings = "NULL")
GCAT_imputed_nochol <- fread("GCAT_imputed_nochol.psam", na.strings = "NULL")
GCAT_imputed_IID <- left_join(GCAT__imputed_nochol, matriu_all_IID)
write.table(GCAT_imputed_IID, "GCAT_imputed_IID.psam", sep = " ", row.names= F, quote = F)

####################

apply(GCAT_imputed_IID[,10:64],2,function(x)table(x))


#Create a file with the data that doen't pass the filter= INFO < 0.7 and minor allele frequency < 0.01
exclude <- filter(pvar, INFO < 0.7 | maf < 0.01)
#Safe the exclude file
write.table(exclude, "/imppc/labs/dnalab/share/ana_test/output/gwas/exclude.txt", sep = "\t", row.names= F, quote = F)

######################
#PLINK
#Taula
#We first safe the directories we will need
plink2=/imppc/labs/dnalab/nblaym/gwas/plink2
output=/imppc/labs/dnalab/share/ana_test/output
output2=/imppc/labs/dnalab/share/ana_test/output/gwas
input=/imppc/labs/dnalab/share/ana_test/output/GCAT_imputed_QC
prefix='GCAT_imputed'
exclude=/imppc/labs/dnalab/share/ana_test/input/exclude.txt

#Use the exclude file to take off the variants that doesn't pass the  filter
cd $output
$plink2 --pfile $input --maf 0.01 --exclude $exclude --make-pgen --out GCAT_imputed_QC --threads 7

#safe the "IIDs"
names_pheno = paste(colnames(matriu_ultima)[-1],collapse = ",")

#We perforom the gwas analysis for the first 22 chromosomes by using plink
$plink2 --pfile $input --pheno-name TIPO041,TIPO242,TIPO244,TIPO246,TIPO250,TIPO268,TIPO272,TIPO274,TIPO275,TIPO280,TIPO285,TIPO296,TIPO300,TIPO311,TIPO346,TIPO365,TIPO401,TIPO427,TIPO432,TIPO459,TIPO472,TIPO493,TIPO496,TIPO530,TIPO536,TIPO537,TIPO553,TIPO558,TIPO564,TIPO569,TIPO600,TIPO602,TIPO626,TIPO627,TIPO692,TIPO696,TIPO704,TIPO710,TIPO715,TIPO716,TIPO719,TIPO722,TIPO724,TIPO728,TIPO729,TIPO733,TIPO780,TIPO784,TIPO788,TIPO799,TIPO995,TIPOC43,TIPOC44,TIPOC50,TIPOV25 --covar-name GENDER, PC1, PC2, PC3, PC4, AGE --ci 0.95 --glm firth cols=chrom,pos,ref,alt,a1count,a1countcc,nobs,a1freq,a1freqcc,test,orbeta,se,ci,tz,p hide-covar --out $output2/gwas --threads 7

#We make this analysis for the X chromosome as make it with all chromosomes together may lead to errors
$plink2 --pfile $input --chr 23 --pheno-name TIPO041,TIPO242,TIPO244,TIPO246,TIPO250,TIPO268,TIPO272,TIPO274,TIPO275,TIPO280,TIPO285,TIPO296,TIPO300,TIPO311,TIPO346,TIPO365,TIPO401,TIPO427,TIPO432,TIPO459,TIPO472,TIPO493,TIPO496,TIPO530,TIPO536,TIPO537,TIPO553,TIPO558,TIPO564,TIPO569,TIPO600,TIPO602,TIPO626,TIPO627,TIPO692,TIPO696,TIPO704,TIPO710,TIPO715,TIPO716,TIPO719,TIPO722,TIPO724,TIPO728,TIPO729,TIPO733,TIPO780,TIPO784,TIPO788,TIPO799,TIPO995,TIPOC43,TIPOC44,TIPOC50,TIPOV25 --covar-name PC1, PC2, PC3, PC4, AGE --xchr-model 1  --ci 0.95 --glm firth  cols=chrom,pos,ref,alt,a1count,a1countcc,nobs,a1freq,a1freqcc,test,orbeta,se,ci,tz,p hide-covar --out $output2/gwas_x --threads 7

####   R
library(data.table)
library(dplyr)
library(ggplot2)
library(qqman)
library(xlsx)
library(qqman)
library(ggplot2)

#uncompress the outputs
TIPO041,TIPO242,TIPO244,TIPO246,TIPO250,TIPO268,TIPO272,TIPO274,TIPO275,TIPO280,TIPO285,TIPO296,TIPO300,TIPO311,TIPO346,TIPO365,TIPO401,TIPO427,TIPO432,TIPO459,TIPO472,TIPO493,TIPO496,TIPO530,TIPO536,TIPO537,TIPO553,TIPO558,TIPO564,TIPO569,TIPO600,TIPO602,TIPO626,TIPO627,TIPO692,TIPO696,TIPO704,TIPO710,TIPO715,TIPO716,TIPO719,TIPO722,TIPO724,TIPO728,TIPO729,TIPO733,TIPO780,TIPO784,TIPO788,TIPO799,TIPO995,TIPOC43,TIPOC44,TIPOC50,TIPOV25
|gzip -c > gwas_all_$TIPO.glm.firth 

gwas <- fread("gwas.TIPO242.glm.firth")
sum(is.na(gwas$P))
grep "NA" gwas.TIPO242.glm.firth

#8       115421653       8-115421653     GTT     G       G       0       0       0       0       0       0       ADD     4988    NA      NA      NA      NA      NA      NA
#9       104764062       9-104764062     A       G       G       0       0       0       0       0       0       ADD     4988    NA      NA      NA      NA      NA      NA
#15      75782973        15-75782973     C       G       G       0       0       0       0       0       0       ADD     4988    NA      NA      NA      NA      NA      NA

#keep the columns of interest and safe them
grep -v ^#CHROM $input | cut -f 1,2,4,5,20 > /imppc/labs/dnalab/share/ana_test/output/cols_yes

#Add the columns as the header
grep ^#CHROM gwas.TIPO041.glm.firth | cut -f 1,2,3,4,5,20 > header

cat $gwas.TIPO $input1x > $output1

696 728 729 733
TIPO041 TIPO242 TIPO244 TIPO246 TIPO250 TIPO268 TIPO272 TIPO274 TIPO275 TIPO280 TIPO285 TIPO296 TIPO300 TIPO311 TIPO346 TIPO365 TIPO401 

cd /imppc/labs/dnalab/share/ana_test/output/gwas
#we are doing a for for all the diseases
for TIPO in TIPO041,TIPO242,TIPO244,TIPO246,TIPO250,TIPO268,TIPO272,TIPO274,TIPO275,TIPO280,TIPO285,TIPO296,TIPO300,TIPO311,
TIPO346,TIPO365,TIPO401,TIPO427,TIPO432,TIPO459,TIPO472,TIPO493,TIPO496,TIPO530,TIPO536,TIPO537,TIPO553,TIPO558,TIPO564,TIPO569,
TIPO600,TIPO602,TIPO626,TIPO627,TIPO692,TIPO696,TIPO704,TIPO710,TIPO715,TIPO716,TIPO719,TIPO722,TIPO724,TIPO728,TIPO729,TIPO733,TIPO780,
TIPO784,TIPO788,TIPO799,TIPO995,TIPOC43,TIPOC44,TIPOC50,TIPOV25
do
gunzip gwas.$TIPO.glm.firth.gz #uncompress the files
cat gwas.$TIPO.glm.firth gwas_x.$TIPO.glm.firth > gwas_all_$TIPO.glm.firth #join the chromosomes files
grep -v ^#CHROM gwas_all_$TIPO.glm.firth | sort -k1,1V -k2,2n | cut -f 1,2,3,4,5,20 > gwas_all_sorted_$TIPO.glm.firth #sort the files
cat header gwas_all_sorted_$TIPO.glm.firth > gwas_output_$TIPO.glm.firth  
awk 'NR==FNR{A[$1]++;next}!($3 in A)' exclude_ok.txt gwas_output_$TIPO.glm.firth | gzip -c > gwas_QC_$TIPO.glm.firth.gz #take off the variants of the exclude files
rm gwas_all_$TIPO.glm.firth #remove the output files that we don't need
rm gwas_all_sorted_$TIPO.glm.firth
rm gwas_output_$TIPO.glm.firth 
echo $TIPO
done

#Change the names of the files
for TIPO in TIPO041,TIPO242,TIPO244,TIPO246,TIPO250,TIPO268,TIPO272,TIPO274,TIPO275,TIPO280,TIPO285,TIPO296,TIPO300,TIPO311,
TIPO346,TIPO365,TIPO401,TIPO427,TIPO432,TIPO459,TIPO472,TIPO493,TIPO496,TIPO530,TIPO536,TIPO537,TIPO553,TIPO558,TIPO564,TIPO569,
TIPO600,TIPO602,TIPO626,TIPO627,TIPO692,TIPO696,TIPO704,TIPO710,TIPO715,TIPO716,TIPO719,TIPO722,TIPO724,TIPO728,TIPO729,TIPO733,TIPO780,
TIPO784,TIPO788,TIPO799,TIPO995,TIPOC43,TIPOC44,TIPOC50,TIPOV25
do
mv gwas_QC_TIPO$TIPO.glm.firth.gz $TIPO.gz
echo $TIPO
done

#keep the chromosome
grep -v ^#CHROM gwas_x.TIPO472.glm.firth > gwas_472
  #join the two files
cat gwas.TIPO472.glm.firth gwas_x.TIPO472.glm.firth > gwas_all_472.glm.firth

#Results look good?
awk 'NR==FNR{A[$1]++;next}$3' include.txt gwas.TIPO041.glm.firth | head #YES
