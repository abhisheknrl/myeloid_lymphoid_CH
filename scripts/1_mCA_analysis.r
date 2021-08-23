
## Preparing data for mCA analysis

library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)

ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")

## UKB subset data
ukbFile<-"data/ukb_subset.tab"
## colnames mapping
nameCodeFile<-"meta/name_code_map.txt"
## ICD10 codes
hemeICDFile<-"meta/heme_malignancies_icd10.txt"
hemeICDFile9<-"meta/heme_malignancies_icd9.txt"
## all WES samples
wesSamp<-"meta/WES_samples.txt"
## panel of normals used for Mutect2
ponSamp<-"meta/UKB_PON.txt"
## UKB withdrawn participants
withdrawlList<-"meta/withdrawl_list.txt"

## Sampels that failed mCA calling
mCAExclusion<-"mCA_calls/sample_inclusion_Ebert.txt"

## mCA calls
cnvFile<-"mCA_calls/CNV_data_mapped_annot_500K.txt"


## M-CHIP and L-CHIP genes
mygeneFile<-"CHIP_calls/CHIP_genes_18June2020.txt"
lygeneFile<-"CHIP_calls/Lymphoid_geneList.txt"

## Updated death data
deathFile<-"meta/death.txt"
codFile<-"meta/death_cause.txt"

## Related individuals
relatedFile<-"meta/related_individuals.dat"

## BMI
bmiFile<-"data/body_measurement.tab"

figDir<-""


##===================================================================
##===================================================================
##ICD10codes
hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")
hicd9<-read.delim(hemeICDFile9,head=TRUE,sep="\t")


##===================================================================
## Name code

nc<-read.table(nameCodeFile,head=FALSE,sep="\t")
cName<-as.character(nc[,1])
names(cName)<-paste("f.",as.character(nc[,2]),sep="")

##===================================================================
## PON

#slist<-readLines(wesSamp)
slist<-paste(figDir,"1_WES_eidlist.txt",sep="")
w<-unique(readLines(withdrawlList))

## LOAD UKB data
ukb<-data.frame(fread(ukbFile))

CN<-colnames(ukb)
CN1<-CN[!CN %in% cName]
names(CN1)<-CN1
cName<-c(cName,CN1)

colnames(ukb)<-cName[colnames(ukb)]


## Withdrawn samples
ukb<-ukb[!ukb$eid %in% w,]

## mCA detection inclusion
mCAINC<-unique(readLines(mCAExclusion))
ukb<-ukb[ukb$eid %in% mCAINC,]


## Sex
ukb<-ukb[!is.na(ukb$genetic_sex),]
ukb<-ukb[ukb$Sex == ukb$genetic_sex,]


#=====================================================================
## BMI

bmi<-read.delim(bmiFile,head=TRUE,sep="\t")
BMI<-bmi$f.21001.0.0
names(BMI)<-as.character(bmi$f.eid)

ukb$BMI<-BMI[as.character(ukb$eid)]


##===================================================================
## Cancer registry data

hicdlist<-as.character(hicd[hicd$Group %in% c("Myeloid","Lymphoid","Unspecified"),]$ICD)
hicdlist9<-as.character(hicd9[hicd9$Disease %in% c("Myeloid","Lymphoid","Unspecified"),]$Code)

## Diagnosis with dates
getDiag<-function(ukb,icdlist,columnNames,dateCols){
	K<-1
	tdf<-c()
	for (i in 1:length(columnNames)){
		CN<-columnNames[i]
		DCN<-dateCols[i]
		#print (CN)
		temp<-ukb[!is.na(ukb[,c(CN)]),]
		if (nrow(temp) == 0){
			next
		}
		temp<-temp[temp[,c(CN)] %in% icdlist,]
		if (nrow(temp)==0){
			next
		}
		## Get the date and eids
		tempdf<-data.frame(eid = temp$eid,
			ICD = temp[,c(CN)],
			Date = temp[,c(DCN)],
			DOAAC = temp$date_attending_assessment_center)
		if (K == 1){
			K<-2
			tdf<-tempdf
		}else{
			tdf<-rbind(tdf,tempdf)
		}
	}
	return (tdf)
}

## Cancer codes
getCancer<-function(ukb,icdlist,columnNames,dateCols){
	K<-1
	tdf<-c()
	for (i in 1:length(columnNames)){
		CN<-columnNames[i]
		DCN<-dateCols[i]
		#print (CN)
		temp<-ukb[!is.na(ukb[,c(CN)]),]
		if (nrow(temp) == 0){
			next
		}
		if (length(icdlist) > 0){
			temp<-temp[temp[,c(CN)] %in% icdlist,]
		}
		if (nrow(temp)==0){
			next
		}
		## Get the date and eids
		tempdf<-data.frame(eid = temp$eid,
			ICD = temp[,c(CN)],
			Date = temp[,c(DCN)],
			DOAAC = temp$date_attending_assessment_center)
		if (K == 1){
			K<-2
			tdf<-tempdf
		}else{
			tdf<-rbind(tdf,tempdf)
		}
	}
	return (tdf)
}


## Get unique samples
getUnique<-function(tdf){
	tdf$Diff<-as.numeric(difftime(tdf$Date, tdf$DOAAC,units=c("days")))
	## Finding the first diagnosis
	T<-table(tdf$eid)
	t1<-tdf[tdf$eid %in% names(T[T==1]),]
	t2<-tdf[tdf$eid %in% names(T[T>1]),]
	for (ID in names(T[T>1])){
		temp<-tdf[tdf$eid == ID,]
		temp<-temp[order(temp$Diff,decreasing=FALSE),]
		tnk<-data.frame(eid = ID,
			ICD = temp$ICD[1],
			Date = temp$Date[1],
			DOAAC = temp$DOAAC[1],
			Diff = temp$Diff[1])
		t1<-rbind(t1,tnk)
	}
	return (t1)
}



cancer10<-c("type_of_cancer_icd10",paste("type_of_cancer_icd10",1:16,sep="_"))
cancerDate10<-c("date_cancer_diagnosed",paste("date_cancer_diagnosed",1:16,sep="_"))

any10<-paste("f.41270.0.",0:212,sep="")
any9<-paste("f.41271.0.",0:46,sep="")
date10<-paste("f.41280.0.",0:212,sep="")
date9<-paste("f.41281.0.",0:46,sep="")


tdf1<-getDiag(ukb,hicdlist,any10,date10)
tdf2<-getCancer(ukb,hicdlist,cancer10,cancerDate10)
tdf3<-getCancer(ukb,hicdlist9,any9,date9)
tdf<-getUnique(rbind(tdf1,tdf2,tdf3))
ukb$Heme<-rep(0,nrow(ukb))
ukb[ukb$eid %in% tdf[tdf$Diff <= 365.25*0.5,]$eid,]$Heme<-1



#==========================================================
## CNV data
cnv<-read.table(cnvFile,head=TRUE,sep="\t")
cnv$CELL_FRAC<-as.numeric(as.character(cnv$CELL_FRAC))

ukb$CNV<-rep(0,nrow(ukb))
ukb[ukb$eid %in% cnv[cnv$CHR %in% 1:22,]$eid,]$CNV<-1
ukb$CNV_size<-ukb$CNV
ukb[ukb$eid %in% cnv[cnv$CHR %in% 1:22 & (!is.na(cnv$CELL_FRAC)) & cnv$CELL_FRAC >= 0.1,]$eid,]$CNV_size<-2
ukb$CNVX<-rep(0,nrow(ukb))
ukb[ukb$eid %in% cnv[cnv$CHR %in% c("X"),]$eid,]$CNVX<-1
ukb$CNVY<-rep(0,nrow(ukb))
ukb[ukb$eid %in% cnv[cnv$CHR %in% c("Y"),]$eid,]$CNVY<-1


mygeneList<-readLines(mygeneFile)
mygeneList<-mygeneList[!mygeneList %in% c("MYD88","STAT5B","STAT3")]

## Annotating CNV data with the genes
### https://ashpublications.org/blood/article/115/14/2731/27247/Copy-neutral-loss-of-heterozygosity-a-novel

myg<-c()

for (i in 1:nrow(cnv)){
	copy<-as.character(cnv$COPY_CHANGE)[i]
	if (copy != "neutral"){
		myg[i]<-""
		next
	}
	if (cnv$CHR[i] %in% c("X","Y")){
		myg[i]<-""
		next
	}
	gl<-unlist(strsplit(as.character(cnv$Genes)[i],"[,]"))
	I<-intersect(gl,mygeneList)
	if ("DNMT3A" %in% I){
		myg[i]<-"DNMT3A"
	}else if ("TET2" %in% I){
		myg[i]<-"TET2"
	}else if ("JAK2" %in% I){
		myg[i]<-"JAK2"
	}else if ("TP53" %in% I){
		myg[i]<-"TP53"
	}else if ("KIT" %in% I){
		myg[i]<-"KIT"
	}else if ("GNB1" %in% I){
		myg[i]<-"MPL/GNB1"
	}else if ("EP300" %in% I){
		myg[i]<-"EP300"
	}else if ("CBL" %in% I){
		myg[i]<-"CBL"
	}else if ("CTCF" %in% I){
		myg[i]<-"CTCF"
	}else if ("EZH2" %in% I){
		myg[i]<-"EZH2"
	}else{
		myg[i]<-paste(I,collapse=",")
	}
}

cnv$myGenes<-myg


## Lymphoid genes
lygeneList<-readLines(lygeneFile)
lygeneList<-lygeneList[!lygeneList %in% mygeneList]

lyg<-c()
for (i in 1:nrow(cnv)){
	copy<-as.character(cnv$COPY_CHANGE)[i]
	if (copy != "neutral"){
		lyg[i]<-""
		next
	}
	if (cnv$CHR[i] %in% c("X","Y")){
		lyg[i]<-""
		next
	}
	gl<-unlist(strsplit(as.character(cnv$Genes)[i],"[,]"))
	I<-intersect(gl,lygeneList)
	if ("ATM" %in% I){
		lyg[i]<-"ATM"
	}else if ("CD79B" %in% I){
		lyg[i]<-"CD79B"
	}else if ("CD79A" %in% I){
		lyg[i]<-"CD79A"
	}else if ("NOL9" %in% I){
		lyg[i]<-"NOL9"
	}else if ("ITPKB" %in% I | "NTRK1" %in% I){
		lyg[i]<-"ITPKB/NTRK1"
	}else if ("NCOR2" %in% I){
		lyg[i]<-"NCOR2"
	}else if ("MIR16-1" %in% I){
		lyg[i]<-"MIR16-1"
	}else if ("TCL1A" %in% I){
		lyg[i]<-"TCL1A"
	}else if ("CIITA" %in% I){
		lyg[i]<-"CIITA"
	}else if ("KMT2C" %in% I){
		lyg[i]<-"KMT2C"
	}else if ("NOTCH1" %in% I){
		lyg[i]<-"NOTCH1"
	}else if ("CDKN2A" %in% I){
		lyg[i]<-"CDKN2A"
	}else if ("CHD2" %in% I){
		lyg[i]<-"CHD2"
	}else if ("ZNF217" %in% I){
		lyg[i]<-"ZNF217"
	}else{
		lyg[i]<-paste(I,collapse=",")
	}
}

cnv$lyGenes<-lyg

#=====================================================================
#=====================================================================
## Myeloid
## CDKN2A ==> JAK2
## TCL1A ==> Myeloid, 14q LOH (CML)

## Overlapping
## GNB1/MPL -- NOL9
## CBL -- ATM
## TP53 -- TP53

## Lymphoid
## .. ==> NCOR2
## FLT3 ==> MIR16-1


LYG<-c("ATM","CIITA","ITPKB/NTRK1","KMT2C","MIR16-1","NCOR2","NOTCH1")
MYG<-c("CBL","CTCF","EP300","JAK2","TP53","MPL/GNB1")

## canonical myeloid/lymphoid
canLy<-c("gain12q","gain15q","gain17q","gain21q","gain22q13.32","gain2p","gain3q","gain8q","gain9q","del10p","del10q","del11q","del13q","del14q","del15q","del17p","del1p","del1q","del22q","del6q","del7q","del8p","tri12","tri18","tri19")
canMy<-c("gain1q","gain21q","gain9p","del12q","del20q","del5q","tri8")


## del17p is common in lymphoid not in myeloid
## del21q is 1 each
## del7q is also lymphoid
## del22q11.22 (IGLV) (Lymphoid)
## del14q32 (IGHV) (Lymphoid)
## Del14q (before IGHV) (Myeloid)

## gain9p (Myeloid)
## 9q34 gain (Myeloid)
## gain9q (Lymphoid)
## tri8 (Myeloid/lymphoid)
## gain8q (Lymphoid)
## gain5q34 (Lymphoid)
## 3q gain (Lymphoid)
## gain22q13.32 (Lymphoid)
## gain22q (Myeloid)
## gain21q (lymphoid)
## gain2p (Lymphoid)
## tri19 (Lymphoid)
## tri 18 (Lymphoid)
## gain17q (Lymphoid)
## gain15q (IGHV) (Lymphoid)
## gain1q (myeloid) (21-32)



#=====================================================================
##

autocnv<-cnv[cnv$CHR %in% 1:22,]

manMG<-c("TP53","CTCF","MPL/GNB1","CBL")
lcnv<-autocnv[autocnv$lyGenes %in% LYG | autocnv$myGenes %in% manMG | autocnv$Canonical_CA %in% canLy,]
lcnv1<-autocnv[!autocnv$eid %in% lcnv$eid,]
lcnv1<-lcnv1[lcnv1$lyGenes != "" & lcnv1$myGenes == "" & lcnv1$Canonical_CA == "",]
ukb$LCNV<-rep(0,nrow(ukb))
ukb[ukb$eid %in% lcnv$eid,]$LCNV<-1
ukb$LCNV_size<-ukb$LCNV
ukb[ukb$eid %in% lcnv[lcnv$CELL_FRAC >= 0.1,]$eid,]$LCNV_size<-2
ukb$LCNVP<-rep(0,nrow(ukb))
ukb[ukb$eid %in% lcnv1$eid,]$LCNVP<-1

manLG<-c("TCL1A","ATM","HEATR3/PLCG2/IRF8")
mcnv<-autocnv[autocnv$myGenes %in% MYG | autocnv$lyGenes %in% manLG | autocnv$Canonical_CA %in% canMy,]
mcnv1<-autocnv[!autocnv$eid %in% mcnv$eid,]
mcnv1<-mcnv1[mcnv1$myGenes != "" & mcnv1$lyGenes == "" & mcnv1$Canonical_CA == "",]
ukb$MCNV<-rep(0,nrow(ukb))
ukb[ukb$eid %in% mcnv$eid,]$MCNV<-1
ukb$MCNV_size<-ukb$MCNV
ukb[ukb$eid %in% mcnv[mcnv$CELL_FRAC >= 0.1,]$eid,]$MCNV_size<-2
ukb$MCNVP<-rep(0,nrow(ukb))
ukb[ukb$eid %in% mcnv1$eid,]$MCNVP<-1


#=====================================================================
#=====================================================================

## Death data
death<-read.table(deathFile,head=TRUE,sep="\t")
death<-death[death$eid %in% ukb$eid,]
DOD<-as.character(death$date_of_death)
DOD<-as.Date(DOD,"%d/%m/%Y")
names(DOD)<-as.character(death$eid)

death2<-read.table(codFile,head=TRUE,sep="\t")
death2<-death2[death2$eid %in% ukb$eid,]

## Adding the date of death
deathCensorDate<-"2020-03-31"
ukb$dead<-rep(0,nrow(ukb))
ukb[ukb$eid %in% names(DOD),]$dead<-1
ukb$death_censor_date<-DOD[as.character(ukb$eid)]
ukb[is.na(ukb$death_censor_date),]$death_censor_date<-deathCensorDate


## Change the death status to 0 if died after the follow-up time
ukb$Diff<-as.numeric(difftime(deathCensorDate,ukb$death_censor_date,units=c("days")))
ukb[ukb$Diff < 0,]$dead<-0
ukb[ukb$Diff < 0,]$death_censor_date<-deathCensorDate

## Time of death since recruitment
ukb$death_year<-as.numeric(difftime(ukb$death_censor_date,ukb$date_attending_assessment_center,units=c("days")))/365.25


## Adding covariates
ukb$EverSmoke<-ukb$ever_smoked
ukb<-ukb[!is.na(ukb$EverSmoke),]

#=====================================================================
#=====================================================================
#=====================================================================

wesSlist<-readLines(paste(figDir,"1_WES_eidlist.txt",sep=""))

hemeukb<-ukb[ukb$Heme == 1,]
ukb<-ukb[ukb$Heme == 0,]
## Excluding one of the related individuals from the analysis
## relatedness

ukb<-ukb[is.na(ukb$genetic_relatedness_exclusion),]

relat<-read.table(relatedFile,head=TRUE,sep=" ")
relat1<-relat[relat$ID1 %in% ukb$eid & relat$ID2 %in% ukb$eid,]
ukb<-ukb[ukb$genetic_kinship_to_other_participants != "Participant excluded from kinship inference process",]

## Add age info to the related individual df
AGE<-ukb$age_attended_assessment_center
names(AGE)<-as.character(ukb$eid)
relat1$Age1<-AGE[as.character(relat1$ID1)]
relat1$Age2<-AGE[as.character(relat1$ID2)]


excl<-c()
for (i in 1:nrow(relat1)){
	eid1<-relat1$ID1[i]
	eid2<-relat1$ID2[i]
	## Check if already in filter list
	if (eid1 %in% excl | eid2 %in% excl){
		next
	}
	## Check if the sample is present in WES
	if (eid1 %in% wesSlist){
		excl<-c(excl,eid2)
		next
	}else if (eid2 %in% wesSlist){
		excl<-c(excl,eid1)
		next
	}
	age1<-relat1$Age1[i]
	age2<-relat1$Age2[i]
	##
	if (age1 >= age2){
		excl<-c(excl,eid2)
	}else if (age2 > age1){
		excl<-c(excl,eid1)
	}
}


ukb<-ukb[!ukb$eid %in% excl,]

#=====================================================================
## CNV heme malignancies
cnv$Heme<-rep(0,nrow(cnv))
cnv[cnv$eid %in% hemeukb$eid,]$Heme<-1

write.table(cnv,file=paste(figDir,"data/UKB_500K_CNV.txt",sep=""),
	col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


#=====================================================================
## heme malignancies

write.table(hemeukb,file=paste(figDir,"data/UKB_500K_hemeMalignancies.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


## Function to prepare censor dates
annotateCensorDates<-function(tukb,icdList,DOD,death2,DIS){
	forcedCensorDate<-as.Date("2020-03-31")
	deathCensorDate<-as.Date("2020-03-31")
	## diagnosis based on icd codes
	tdf1<-getDiag(tukb,icdList,any10,date10)
	tdf2<-getCancer(tukb,icdList,cancer10,cancerDate10)
	tdf<-getUnique(rbind(tdf1,tdf2))
	Date<-tdf$Date
	names(Date)<-as.character(tdf$eid)
	## censorDate
	tukb$has_disease<-rep(0,nrow(tukb))
	tukb[tukb$eid %in% tdf$eid,]$has_disease<-1
	## adding censor dates
	tukb1<-tukb[tukb$has_disease == 1,]
	tukb1$Censor_date<-as.Date(Date[as.character(tukb1$eid)])
	tukb2<-tukb[!tukb$eid %in% tukb1$eid,]
	tukb2$Censor_date<-as.Date(forcedCensorDate)
	##
	tukb<-rbind(tukb1,tukb2)
	#print (head(tukb$Censor_date))
	#print (head(tukb$death_censor_date))
	## Adding the death date as a censor date
	tukb1<-tukb[tukb$dead == 1,]
	tukb1$Diff<-as.numeric(difftime(tukb1$death_censor_date,tukb1$Censor_date,units=c("days")))
	tukb11<-tukb1[tukb1$Diff < 0,]
	tukb11$Censor_date<-tukb11$death_censor_date
	tukb12<-tukb1[!tukb1$eid %in% tukb11$eid,]
	tukb2<-tukb[!tukb$eid %in% tukb1$eid,]
	##
	tukb<-rbind(tukb2,tukb11[,colnames(tukb2)],tukb12[,colnames(tukb2)])
	## Adding disease
	tukb$Disease<-rep(DIS,nrow(tukb))
	## Return
	return (tukb)
}


icdList<-as.character(hicd[hicd$Group %in% c("Lymphoid","Plasma_cell"),]$ICD)
icdListMyeloid<-as.character(hicd[hicd$Group %in% c("Myeloid"),]$ICD)


## Myeloid malignancies
myheme<-annotateCensorDates(hemeukb,icdListMyeloid,DOD,death2,"Myeloid_malignancies")
myheme<-myheme[myheme$has_disease == 1,]

lyheme<-annotateCensorDates(hemeukb,icdList,DOD,death2,"Lymphoid_malignancies")
lyheme<-lyheme[lyheme$has_disease == 1,]


##
cnv$myHeme<-rep(0,nrow(cnv))
cnv$lyHeme<-rep(0,nrow(cnv))
cnv[cnv$eid %in% myheme$eid,]$myHeme<-1
cnv[cnv$eid %in% lyheme$eid,]$lyHeme<-1

cnv1<-cnv[cnv$myHeme == 1 | cnv$lyHeme == 1,]

write.table(cnv1,file=paste(figDir,"scripts/Final_figs/UKB_prevalent_heme_CNV.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")



#==============================================================
##

## PCA
CN<-colnames(ukb)
CN<-gsub("genetic_pca","PC",CN)
colnames(ukb)<-CN

## Age
ukb$Age<-ukb$age_attended_assessment_center
ukb$Age_sq<-ukb$Age*ukb$Age


write.table(ukb,file=paste(figDir,"data/UKB_500K.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

