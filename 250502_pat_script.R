# PAT script
# this is currently assuming that the discovery GWAS is based on samples of European ancestry. needs to be adapated otherwise.

# INPUT NEEDED
# full summary statistics for EUR discovery GWAS
path_eur_sumstats=/SAN/ugi/mdd/trans_ethnic_ma/results/FE_EUR_MDD2_ASDG_qced_rsid.txt.gz 
# genome-wide summary statistics for your non-EUR target data set
nonEURgwas=/SAN/ugi/mdd/trans_ethnic_ma/results/FE_eas_qced_rsid.txt.gz
# working directory for this: 
WD=/SAN/ugi/mdd/trans_ethnic_ma/kk_testdir
# number of cases and controls for the non-European target GWAS


# optional
# instead of creating our own list of independent snps, one can also put in an existing list of gwas hits, e.g a table from a paper 
gwashits=/SAN/ugi/mdd/trans_ethnic_ma/data/EUR_agds_mvp_208_SNPs_data.txt 



######################################################################
# create list if independent gwas hits 
# this is only necessary if you dont already have a list of independent hits and their gwas results (e.g. a table from the paper)

# this needs an ld ref
ldref=/SAN/ugi/mdd/trans_ethnic_ma/data/1000G/EUR/eur_autosomes_hwe_1e-6_rmdup_rmmissing #the EUR LD reference panel from the 1000 Genomes project. dont change this unless you are using a different ancestry 


# print header. assuming zipped
gzip -dc $path_eur_sumstats | head
# # # # # MODIFY HERE 
# define the column names of the rsid and p-value fields in this file
field_p=P
field_rsid=rsid
/share/apps/plink-1.90-beta-6.10/plink  --clump $path_eur_sumstats --bfile $ldref --clump-p1 5e-8 --clump-r2 0.001 --clump-snp-field $field_rsid --clump-field $field_p --out $WD/clumped_hits --clump-best

# save only the rsids of lead variants
awk 'NR > 1 {print $3}' $WD/clumped_hits.clumped > $WD/clumped_hits.txt 




######################################################################
# for these gwashits, create summary stats files following our naming conventions

/share/apps/R-4.2.2/bin/R
if (!require("tidyverse")) install.packages('tidyverse')
if (!require("data.table")) install.packages('data.table')

# # # # # MODIFY HERE 
wd<-"/SAN/ugi/mdd/trans_ethnic_ma/kk_testdir"


#####################################
#### OPTION1: there is an existing list of independent gwas hits and their results 
path_gwashits<-"/SAN/ugi/mdd/trans_ethnic_ma/data/EUR_agds_mvp_208_SNPs_data.txt"
#path_noneur<-"/SAN/ugi/mdd/trans_ethnic_ma/results/FE_eas_qced_rsid.txt.gz"
data <- read.table(path_gwashits, header = TRUE, sep = "", stringsAsFactors = FALSE)
# Show current column names
cat("\nCurrent column names:\n")
print(colnames(data))
# # # # # MODIFY HERE 
data <- data %>% 
  rename(
  	rsid = rsid,      # Rename rsid column name to rsid
    chr = chr,      # Rename chromosome column name to chr
    pos = pos,      # Rename position column name to pos
    ea = EA,           # Rename effect allele column name to ea
    nea = NEA,         # Rename non-effect (reference) allele to nea
    beta = BETA,       # Rename BETA to beta
    se = SE,           # Rename SE to se
    p = P              # Rename P to p
  ) %>% 
  select(rsid, chr, pos, ea, nea, beta, se, p)
head(data) # View the renamed data
#####################################

#### OPTION2: there is only a genome-wide sum stats file for EUR and you have just created your list of independent hits iwht plink in the step above. 
path_eur_sumstats<-"/SAN/ugi/mdd/trans_ethnic_ma/results/FE_EUR_MDD2_ASDG_qced_rsid.txt.gz" # full summary statistics for EUR discovery GWAS
path_clumped_hits<-paste(wd,"/clumped_hits.txt",sep="") # full summary statistics for EUR discovery GWAS
data <- fread(path_eur_sumstats)
hits<-read.table(path_clumped_hits,header=F,stringsAsFactors=F)
# Show current column names
cat("\nCurrent column names:\n")
print(colnames(data))
# # # # # MODIFY HERE 
data <- data %>% 
  rename(
  	rsid = rsid,      # Rename rsid column name to rsid
    chr = Chromosome,      # Rename chromosome column name to chr
    pos = Position,      # Rename position column name to pos
    ea = EA,           # Rename effect allele column name to ea
      eaf = EAF,           # Rename effect allele frequency column name to eaf
    nea = NEA,         # Rename non-effect (reference) allele to nea
    beta = BETA,       # Rename BETA to beta
    se = SE,           # Rename SE to se
    p = P              # Rename P to p
  ) %>% 
  select(rsid, chr, pos, ea, nea, eaf, beta, se, p)
# View the renamed data
head(data)
dim(data)
table(hits$V1 %in% data$rsid)
data<-data[which(data$rsid %in% hits$V1),]
#####################################



# checks and filters (add as required)
#1. are there any pvalues larger than genomewide significance? (could be ok)
table(data$p>5*10^-8)
if(any(data$p>5*10^-8)){
	print(data[which(data$p>5*10^-8),])
}

#missing any key values?
if(any(is.na(data$p))){
data[which(is.na(data$p)),]
cat("N hits before removing missing p: ")
print(dim(data)[1])
data<-data[which(!is.na(data$p)),]
cat("N hits after removing missing p: ")
print(dim(data)[1])
}

if(any(is.na(data$se))){
data[which(is.na(data$se)),]
cat("N hits before removing missing se: ")
print(dim(data)[1])
data<-data[which(!is.na(data$se)),]
cat("N hits after removing missing se: ")
print(dim(data)[1])
}




# Save to file 
outfile<-paste(wd,"/gwashits.txt",sep="")
write_delim(
  data,
  file = outfile, delim = "\t",na = "NA",quote = "none",
  escape = "none",  col_names = TRUE
)
# save a list with the rsids
outfile<-paste(wd,"/gwashits_rsids.txt",sep="")
write.table(data$rsid, file=outfile,quote=F,row.names=F,col.names=F)




##########################################################
# create credible sets
# this is using plink 

# define your working directory
# # # # # MODIFY HERE 
WD=/SAN/ugi/mdd/trans_ethnic_ma/kk_testdir # project working directory


# path to the file with the rs-ids of the lead variants of previously publsihed hits from the (EUR) GWAS
# this was created above, in r, and saved in the wd 
gwashits=$WD/gwashits_rsids.txt #File name of published GWAS in Europeans


#####1. Create credible set based on Established SNP list with LD calculated based on EUR ref panel

ldref=/SAN/ugi/mdd/trans_ethnic_ma/data/1000G/EUR/eur_autosomes_hwe_1e-6_rmdup_rmmissing #the EUR LD reference panel from the 1000 Genomes project. dont change this unless you are using a different ancestry 

#plink command to generate the credible set for each loci based on LD from European population
/share/apps/plink-1.90-beta-6.10/plink --bfile $ldref --r2 --ld-snp-list $gwashits --ld-window-kb 50 --ld-window 99999 --ld-window-r2 0.8 --out $WD/crediblesets

##reformat the credible set result from plink
sed q $WD/crediblesets.ld > $WD/head.txt
grep -v -f $WD/head.txt $WD/crediblesets.ld > $WD/body.txt
sed 's/SNP_A/SNP.Lead/g' $WD/head.txt > $WD/temp.txt
sed 's/SNP_B/SNP.proxy/g' $WD/temp.txt > $WD/header.txt
sed 's/CHR_B/CHR/g' $WD/header.txt > $WD/temp.txt
sed 's/BP_B/BP/g' $WD/temp.txt > $WD/header.txt
cat $WD/header.txt $WD/body.txt > $WD/crediblesets.txt  
rm $WD/body.txt $WD/header.txt $WD/head.txt $WD/temp.txt 
# The final file with the credible sets is $WD/crediblesets.txt
## Subsequent steps implemented in R








######################################################################
# pat ratios

/share/apps/R-4.2.2/bin/R
if (!require("tidyverse")) install.packages('tidyverse')
if (!require("data.table")) install.packages('data.table')
if (!require("dplyr")) install.packages('dplyr')


# # # # # MODIFY HERE 
path_noneur<-"/SAN/ugi/mdd/trans_ethnic_ma/results/FE_eas_qced_rsid.txt.gz" # full summary statistics for nonEUR target GWAS
wd<-"/SAN/ugi/mdd/trans_ethnic_ma/kk_testdir"
path_eur_sumstats<-"/SAN/ugi/mdd/trans_ethnic_ma/results/FE_EUR_MDD2_ASDG_qced_rsid.txt.gz" # full summary statistics for EUR discovery GWAS
N_ctrl <- 360956 #Number of controls in the non-European (target) GWAS
N_cas <- 21980 #Number of cases in the non-European (target) GWAS

# read in the EUR gwas hits file we created in step 1 
loci<-read.table(paste(wd,"/gwashits.txt",sep=""),header=T,stringsAsFactors=F)
head(loci)




loci$Pthresh=loci$p*100 ## NOTE: Established SNP list summary stats to include additional column with P-value threshold (100*pSentinel) for filtering credible set proxy SNPs
loci$N.COJO=NULL
loci$Location=NULL
loci$IMPACT=NULL
loci$Allele=NULL
dim(loci) #number of loci which we are able to assess

##Read in credible set that we created in the plink step 
path_cred<-paste(wd,"/crediblesets.txt",sep="")
cred=read.table(path_cred, header=T) #EUR credible set file from PLINK, see Transferability.sh
cred$Lead_Proxy="Proxy"
cred$Lead_Proxy=replace(cred$Lead_Proxy, cred$SNP.Lead==cred$SNP.proxy, "Lead") #Assign lead snps 
cred$Lead_Proxy=as.factor(cred$Lead_Proxy) 

loci$Locus <- loci$rsid
cred2=merge(loci[c("rsid", "Locus", "Pthresh")], cred, by.x="rsid", by.y="SNP.Lead", all=T)


#Read in EUR study sumstats and merge with credset to extract credset sumstats and filter based on P
eur <- fread(path_eur_sumstats)
# Show current column names
cat("\nCurrent column names:\n")
print(colnames(data))
# # # # # MODIFY HERE 
eur <- eur %>% 
  rename(
  	rsid = rsid,      # Rename rsid column name to rsid
    chr = Chromosome,      # Rename chromosome column name to chr
    pos = Position,      # Rename position column name to pos
    ea = EA,           # Rename effect allele column name to ea
      eaf = EAF,           # Rename effect allele frequency column name to eaf
    nea = NEA,         # Rename non-effect (reference) allele to nea
    beta = BETA,       # Rename BETA to beta
    se = SE,           # Rename SE to se
    p = P              # Rename P to p
  ) %>% 
  select(rsid, chr, pos, ea, nea, eaf, beta, se, p)
head(eur)# View the renamed data

cred2=merge(cred2, eur, by.x=c("SNP.proxy"), by.y=c("rsid"), all.x=T)
head(cred2)

#cred2 <- cred2[cred2$SNP %in% loci$SNP,] #exclude redundant entries
cred2 <- cred2[!is.na(cred2$p),] #this removes the ones without results in the EUR GWAS. 
length(unique(cred2$Locus))
dim(cred2)


lead=subset(cred2, Lead_Proxy=="Lead")
proxy=subset(cred2, Lead_Proxy=="Proxy")
dim(lead)
dim(proxy)

cred_fin=subset(proxy, proxy$p<proxy$Pthresh | proxy$p <5e-900 | is.na(proxy$p))
cred_fin=rbind(lead,cred_fin)

#cred_fin <- cred_fin[cred_fin$SNP %in% loci$SNP,]
#length(unique(cred_fin$SNP))
#[1] 205

cred_fin$Pthresh=NULL
cred_fin$SNP.y=NULL
cred_fin$A1=NULL
cred_fin$TEST=NULL
cred_fin$A2=NULL
cred_fin$STAT=NULL
cred_fin$NMISS=NULL
cred_fin$CHR_A=NULL
cred_fin$BP_A=NULL
cred_fin$CHR=NULL
cred_fin$BP=NULL

names(cred_fin)[2]="SNP.Lead"
names(cred_fin)[1]="SNP"




#Calculate tranferability P threshold in each non-European ancestry with the Pf factor
#Pf factor is based on imperical estimation, please see manuscript method section for details
Pf_afr <- 0.008341
Pf_eas <- 0.007378
Pf_sas <- 0.006847
Pf_his <- 0.003147

df <- count(cred_fin,Locus)
df$P_thre_afr <- 10^(log(0.05,base=10)- Pf_afr*(df$n-1))
df$P_thre_eas <- 10^(log(0.05,base=10)- Pf_eas*(df$n-1))
df$P_thre_his <- 10^(log(0.05,base=10)- Pf_his*(df$n-1))
df$P_thre_sas <- 10^(log(0.05,base=10)- Pf_sas*(df$n-1))
names(df)[2] <- "N_credset"
cred_fin <- merge(cred_fin,df,all.x=T,order=F)

cred_fin <- cred_fin[order(cred_fin$chr,cred_fin$pos),]
#write.table(cred_fin, "MDD.EUR.credset_final.txt", col.names=T, row.names=F, quote=F, sep= "\t") #save the cred_fin file


# now read in the sum stats from the non-European ancestry study (target data set)
noneur<-fread(path_noneur) ## this could take a while 
noneur <- noneur %>% 
  rename(
  	rsid = rsid,      # Rename rsid column name to rsid
    chr_n = Chromosome,      # Rename chromosome column name to chr
    pos_n = Position,      # Rename position column name to pos
    ea_n = EA,           # Rename effect allele column name to ea
    nea_n = NEA,         # Rename non-effect (reference) allele to nea
    eaf_n = EAF,         # Rename effect allele frequency to eaf
    beta_n = BETA,       # Rename BETA to beta
    se_n = SE,           # Rename SE to se
    p_n = P              # Rename P to p
  ) %>% 
  select(rsid, chr_n, pos_n, ea_n, nea_n, eaf_n, beta_n, se_n, p_n)
  class(noneur)

names(noneur)<-c("rsid", "chr_n", "pos_n", "ea_n", "nea_n", "eaf_n", "beta_n", "se_n", "p_n")



cred_fin_noneur=merge(cred_fin, noneur, by.x="SNP", by.y="rsid", all.x=T) #merge AFR sumtats with the credible set file

dim(cred_fin_noneur)
cred_fin_noneur<- cred_fin_noneur[!is.na(cred_fin_noneur$p_n),] #exclude columns without P values
dim(cred_fin_noneur)

# access whether EA/NEA match between two populations
eur_alleles <- paste0(cred_fin_noneur$ea,cred_fin_noneur$nea) #paste EA and NEA into a new column, to check if EUR and AFR sumstats used the same alleles ad effect alleles
noneur_alleles <- paste0(cred_fin_noneur$ea_n,cred_fin_noneur$nea_n) 
check <- eur_alleles == noneur_alleles
sum(!check) #number of rows where noneur and EUR used different alleles as effect alleles
if(sum(!check)>0){
cred_fin_noneur[!check,]
#### be careful: make sure alleles are shown as a,t,c,g - otherwise recode 
# recode to matching alleles
if(any(cred_fin_noneur[which(cred_fin_noneur$ea==cred_fin_noneur$nea_n),])){
cred_fin_noneur[which(cred_fin_noneur$ea==cred_fin_noneur$ea_n),]$beta_n<-cred_fin_noneur[which(cred_fin_noneur$ea==cred_fin_noneur$ea_n),]$beta_n*(-1)
}
if(any(cred_fin_noneur$ea!=cred_fin_noneur$nea_n & cred_fin_noneur$ea!=cred_fin_noneur$ea_n)){
cred_fin_noneur[which(cred_fin_noneur$ea!=cred_fin_noneur$nea_n & cred_fin_noneur$ea!=cred_fin_noneur$ea_n),]<-NULL
}

length(unique(cred_fin_noneur$Locus)) #number of loci available for transferability analysis
#[1] 162

#Check direction of effect sizes
#Mismatched to be assigned as not transferred
cred_fin_noneur$b.direction=NA
cred_fin_noneur$b.direction=replace(cred_fin_noneur$b.direction, cred_fin_noneur$beta>=0 & cred_fin_noneur$beta_n>=0, "match")
cred_fin_noneur$b.direction=replace(cred_fin_noneur$b.direction, cred_fin_noneur$beta<=0 & cred_fin_noneur$beta_n<=0, "match")
cred_fin_noneur$b.direction=replace(cred_fin_noneur$b.direction, cred_fin_noneur$beta<=0 & cred_fin_noneur$beta_n>=0, "mismatch")
cred_fin_noneur$b.direction=replace(cred_fin_noneur$b.direction, cred_fin_noneur$beta>=0 & cred_fin_noneur$beta_n<=0, "mismatch")


#define transferable loci (using the specific threshold)
# # # # # MODIFY HERE: depdening on ancestry of target GWAS, decide whether to use P_thres_eas or P_thres_afr etc 
# one can also use a specific value, e.g. 5% 
cred_fin_noneur$transferred=NA
cred_fin_noneur$transferred=replace(cred_fin_noneur$transferred, cred_fin_noneur$p_n < cred_fin_noneur$P_thre_eas, "Yes")
cred_fin_noneur$transferred=replace(cred_fin_noneur$transferred, cred_fin_noneur$p_n >= cred_fin_noneur$P_thre_eas, "No")
cred_fin_noneur$transferred=replace(cred_fin_noneur$transferred, cred_fin_noneur$b.direction=="mismatch", "No")





#Calculate power (to derive N expected transferable)

alpha <- 0.05 #set the Type I error rate level for power calculation

N_tot <- N_cas+N_ctrl
n <- N_tot

# ratio cases controls
r <- N_cas/(N_cas+N_ctrl)

# proportion of exposed controls african
p_ctrl <- cred_fin_noneur$eaf_n
f <- p_ctrl

# # # # beta europeans
b <- cred_fin_noneur$beta

phi <- r 

POWER<-pchisq(qchisq(alpha,df=1,lower = F), df=1, ncp = 2*f*(1-f)*n*phi*(1-phi)*b^2, lower = F)

#sum(POWER,na.rm=T) #sum up power and divide by number of loci/variants for total expected power
#[1] 2376.05
#2376.05/204
#[1] 11.6473

cred_fin_noneur$power=POWER

cred_fin_noneur$powered=NA
cred_fin_noneur$powered=replace(cred_fin_noneur$powered, cred_fin_noneur$power>=0.8, "Yes") #power greater than 80% defined as sufficient statistical power
cred_fin_noneur$powered=replace(cred_fin_noneur$powered, cred_fin_noneur$power<0.8, "No")

outfile<-paste(wd,"/trasferability_results_crediblesets.txt",sep="")
write.table(cred_fin_noneur, outfile, col.names=T, row.names=F, quote=F, sep= "\t") #save the intermediate data frame



## Define transferable loci and ancestry specific loci
## Define transferable loci
my_loci <- unique(cred_fin_noneur$Locus) #extract the vector of unique loci available for assessment in AFR ancestry

trans <- c() # create the vector for transferable loci
# loop through each loci to look for transferable ones
for (i in 1:length(my_loci)){
t <- cred_fin_noneur[cred_fin_noneur$Locus == my_loci[i],]
if(sum(t$transferred=="Yes")>0) {trans <- c(trans,my_loci[i])}
}

## Define specific
spec <- c() # create the vector for ancestry specific loci
#Loop through each loci to look for non-transferable ones

for (i in 1:length(my_loci)){
t <- cred_fin_noneur[cred_fin_noneur$Locus == my_loci[i],]
# # # # # MODIFY HERE: depdening on ancestry of target GWAS, decide whether to use P_thres_eas or P_thres_afr etc 
if(sum(t$powered=="Yes")>0 & sum(t$p_n<t$P_thre_eas)==0) {spec <- c(spec,my_loci[i])} #loci with at least one powered variants, while have zero transferable variant will be defined as ancestry specific (i.e., non-transferable)
}

intersect(trans,spec) # check the intersection between trans and spec vectors, expected to be 0
#[1] 0

## Finalize transferability results 
loci$transferable <- "underpowered" #the final column for transferability result in target population
loci$transferable[loci$rsid %in% trans] <- "transferable"
loci$transferable[loci$rsid %in% spec] <- "EURspecific"
table(loci$transferable)


outfile<-paste(wd,"/trasferability_results_leadsnps.txt",sep="")
write.table(loci, outfile, col.names=T, row.names=F, quote=F, sep= "\t") #save the intermediate data frame


trans[!trans %in% unique(cred_fin_noneur[which(cred_fin_noneur$transferred=="Yes"),]$Locus)]

##Power Adjusted Transferability (PAT) ratio calculation


cred_fin_noneur[which(cred_fin_noneur$Lead_Proxy),]

loci <- fread("GWAS_EUR_published.txt") #list of GWAS significant loci in EUR
afr <- fread("EUR_loci_transferability_AFR_final.txt") #transferablity result with all SNPs available in African sumstats

afr <- afr[afr$Locus %in% loci$Locus,] #only keep loci available in both Africans and Europeans


max_power <- aggregate(power ~ Locus, data = cred_fin_noneur, max) #retrive SNP with maximum power in each locus
expected<-sum(max_power$power) #sum of maximum power in each locus: this is the expected number of transferable loci
table(loci$transferable)
observed<-length(unique(cred_fin_noneur[which(cred_fin_noneur$transferred=="Yes"),]$Locus))

cat("The PAT ratio is: ")
round(observed/expected,3)









