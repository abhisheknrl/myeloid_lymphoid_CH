
## Script to test association between CHIP categories and malignancy sub-types
## WESL, WESM, L500, M500, FIG

args = commandArgs(trailingOnly=TRUE)

CMD<-args[1]



library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)

library(survival)
library(survminer)

ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"

##===================================================================
##ICD10codes
hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")


#=====================================================================
## Reload data
## Function to extract samples with diagnosis

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
date10<-paste("f.41280.0.",0:212,sep="")


#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================

if (CMD == "WES"){


## Load data
ukbc<-data.frame(fread(paste(figDir,"1_UKB_data.txt",sep="")))


## Adding covariates
ukbc$EverSmoke<-ukbc$ever_smoked
## PCA
CN<-colnames(ukbc)
CN<-gsub("genetic_pca","PC",CN)
colnames(ukbc)<-CN
## Age
ukbc$Age<-ukbc$age_attended_assessment_center
ukbc$Age_sq<-ukbc$Age*ukbc$Age

ukbc$Age_categ<-rep("40-49",nrow(ukbc))
ukbc[ukbc$age_attended_assessment_center >= 50 & ukbc$age_attended_assessment_center < 60,]$Age_categ<-"50-59"
ukbc[ukbc$age_attended_assessment_center >= 60,]$Age_categ<-">60"
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("40-49","50-59",">60"))


ukbc<-ukbc[!is.na(ukbc$EverSmoke),]

## Caucasian
ukbc$Caucasian<-rep(1,nrow(ukbc))
ukbc[is.na(ukbc$genetic_ethnic_grouping),]$Caucasian<-0


## Convert into factors

ukbc$CHIP<-as.factor(ukbc$CHIP)
ukbc$CHIP_size<-as.factor(ukbc$CHIP_size)
ukbc$MCHIP<-as.factor(ukbc$MCHIP)
ukbc$Msize<-as.factor(ukbc$Msize)


#=====================================================================
## Enrichment of specific lymphoid malignancies
DOD<-c()
death2<-c()
forcedCensorDate<-"2020-03-31"
deathCensorDate<-"2020-03-31"



icdList<-as.character(hicd[hicd$Group %in% c("Lymphoid","Plasma_cell"),]$ICD)
CLL<-c("C911","C830")
PC<-c("C900","D472")
HL<-c("C81",paste("C81",0:9,sep=""))
FL<-c("C82",paste("C82",0:9,sep=""))
DLBCL<-c("C833")
NHL<-c("C85",paste("C85",0:9,sep=""))
WS<-c("C880")
CIRC<-c("C901","C910","C913","C914","C831","C841","C915","C916")
nhl<-icdList[!icdList %in% c(PC,CLL,HL,FL,DLBCL,NHL,WS,CIRC)]

## Myeloid
icdList<-as.character(hicd[hicd$Group %in% c("Myeloid"),]$ICD)
mds_icdList<-c(paste("D46",0:9,sep=""),"C946")
aml_icdList<-paste("C92",0,sep="")
mpn_icdList<-c("D752","D45","C945","D471","D473")
other_myeloid<-icdList[!icdList %in% c(mds_icdList,aml_icdList,mpn_icdList)]


## Prepare status and follow-up date in the function
getStatusYear<-function(ukbc,allICD,icdlist,EXCL=""){
	## Lymphoid malignancies
	#tdf1<-getDiag(ukbc,c(icdlist,EXCL),any10,date10)
	#tdf2<-getCancer(ukbc,c(icdlist,EXCL),cancer10,cancerDate10)
	tdf1<-getDiag(ukbc,allICD,any10,date10)
	tdf2<-getCancer(ukbc,allICD,cancer10,cancerDate10)
	tempdf<-rbind(tdf1,tdf2)
	#exclEid<-tempdf[tempdf$ICD %in% EXCL,]$eid
	#tempdf<-tempdf[!tempdf$eid %in% exclEid,]
	tempdf$Diff2<-as.numeric(difftime(forcedCensorDate,tempdf$Date,units=c("days")))
	tempdf<-tempdf[tempdf$Diff2 >= 0,]
	tdf<-getUnique(tempdf[,colnames(tdf1)])
	tdf<-tdf[tdf$ICD %in% icdlist,]
	## Positive cases
	cenDate<-tdf$Date
	names(cenDate)<-as.character(tdf$eid)
	## dead individuals
	t1<-ukbc[(!ukbc$eid %in% tdf$eid) & ukbc$dead == 1,c("eid","death_censor_date")]
	cenDate1<-t1$death_censor_date
	names(cenDate1)<-as.character(t1$eid)
	## individuals who did not die
	t2<-ukbc[(!ukbc$eid %in% tdf$eid) & ukbc$dead == 0,]$eid
	cenDate2<-rep(forcedCensorDate,length(t2))
	names(cenDate2)<-as.character(t2)
	## Has disease
	Status<-c(rep(1,length(cenDate)),rep(0,length(cenDate1)+length(cenDate2)))
	names(Status)<-c(names(cenDate),names(cenDate1),names(cenDate2))
	DATE<-c(as.character(cenDate),as.character(cenDate1),cenDate2)
	return (list(Status,DATE))
}


#=====================================================================
## Cox model
MYL<-list(CLL,PC,HL,FL,DLBCL,NHL,WS,CIRC,nhl,mds_icdList,aml_icdList,mpn_icdList,other_myeloid)
names(MYL)<-c("CLL","MM","HL","FL","DLBCL","NHL","WS","CIRC","Other_lymphoid","MDS","AML","MPN","Other_myeloid")

allICD<-c(CLL,PC,HL,FL,DLBCL,NHL,WS,CIRC,nhl,mds_icdList,aml_icdList,mpn_icdList,other_myeloid)

K<-1
for (MALIG in names(MYL)){
	print (MALIG)
	dlist<-MYL[[MALIG]]
	ukbc$Status<-rep(0,nrow(ukbc))
	ukbc$Date<-rep(forcedCensorDate,nrow(ukbc))
	ukbc$Year<-NULL
	if (MALIG %in% c("MM","HL","FL","DLBCL","NHL","WS","CIRC","Other_lymphoid")){
		SD<-getStatusYear(ukbc,allICD,dlist,CLL)
	}else{
		SD<-getStatusYear(ukbc,allICD,dlist,"")
	}
	ukbc$Status<-SD[[1]][as.character(ukbc$eid)]
	ukbc$Date<-SD[[2]][as.character(ukbc$eid)]
	ukbc$Year<-as.numeric(difftime(ukbc$Date,ukbc$date_attending_assessment_center,units=c("days")))/365.25
	q95<-as.numeric(quantile(ukbc$Year,0.95))
	ukbc[ukbc$Year > q95,]$Status<-0
	ukbc[ukbc$Year > q95,]$Year<-q95
	## Now do a cox model for L-CHIP
	NE<-nrow(ukbc[ukbc$CHIP == 1 & ukbc$Status == 1 & ukbc$MCHIP == 0,])
	if (NE >= 4){
		cox1<-coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CHIP, data = ukbc[ukbc$MCHIP == 0,])
		x1<-summary(cox1)
		mvdf1<-as.data.frame(x1$conf.int)
		colnames(mvdf1)<-c("Haz_ratio","negHR","CI_low","CI_high")
		mvdf1$P<-x1$coefficients[,5]
		mvdf1$Variable<-rownames(mvdf1)
		##
		Variable<-"CHIP1"
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "L-CHIP",
			Total = nrow(ukbc[ukbc$CHIP == 1 & ukbc$MCHIP == 0,]),
			Event = nrow(ukbc[ukbc$CHIP == 1 & ukbc$Status == 1 & ukbc$MCHIP == 0,]),
			HR = mvdf1[mvdf1$Variable %in% Variable,]$Haz_ratio,
			CI_low = mvdf1[mvdf1$Variable %in% Variable,]$CI_low,
			CI_high = mvdf1[mvdf1$Variable %in% Variable,]$CI_high,
			P = mvdf1[mvdf1$Variable %in% Variable,]$P)
	}else{
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "L-CHIP",
			Total = nrow(ukbc[ukbc$CHIP == 1 & ukbc$MCHIP == 0,]),
			Event = nrow(ukbc[ukbc$CHIP == 1 & ukbc$Status == 1 & ukbc$MCHIP == 0,]),
			HR = NA,
			CI_low = NA,
			CI_high = NA,
			P = NA)
	}
	if (K == 1){
		mydf<-tempdf
		K<-2
	}else{
		mydf<-rbind(mydf,tempdf)
	}
	## Now do a cox model for M-CHIP
	NE<-nrow(ukbc[ukbc$MCHIP == 1 & ukbc$Status == 1 & ukbc$CHIP == 0,])
	if (NE >= 4){
		cox1<-coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + MCHIP, data = ukbc[ukbc$CHIP == 0,])
		x1<-summary(cox1)
		mvdf1<-as.data.frame(x1$conf.int)
		colnames(mvdf1)<-c("Haz_ratio","negHR","CI_low","CI_high")
		mvdf1$P<-x1$coefficients[,5]
		mvdf1$Variable<-rownames(mvdf1)
		##
		Variable<-"MCHIP1"
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "M-CHIP",
			Total = nrow(ukbc[ukbc$MCHIP == 1 & ukbc$CHIP == 0,]),
			Event = nrow(ukbc[ukbc$MCHIP == 1 & ukbc$Status == 1 & ukbc$CHIP == 0,]),
			HR = mvdf1[mvdf1$Variable %in% Variable,]$Haz_ratio,
			CI_low = mvdf1[mvdf1$Variable %in% Variable,]$CI_low,
			CI_high = mvdf1[mvdf1$Variable %in% Variable,]$CI_high,
			P = mvdf1[mvdf1$Variable %in% Variable,]$P)
	}else{
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "M-CHIP",
			Total = nrow(ukbc[ukbc$MCHIP == 1 & ukbc$CHIP == 0,]),
			Event = nrow(ukbc[ukbc$MCHIP == 1 & ukbc$Status == 1 & ukbc$CHIP == 0,]),
			HR = NA,
			CI_low = NA,
			CI_high = NA,
			P = NA)
	}
	if (K == 1){
		mydf<-tempdf
		K<-2
	}else{
		mydf<-rbind(mydf,tempdf)
	}
}


write.table(mydf,file=paste(figDir,"scripts/data/WES_malignancy_subtypes_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

}

#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
## 500K cases

if (CMD == 500){

## Load data
ukbc<-data.frame(fread(paste(figDir,"data/UKB_500K.txt",sep="")))

## Adding covariates
ukbc$EverSmoke<-ukbc$ever_smoked
## PCA
CN<-colnames(ukbc)
CN<-gsub("genetic_pca","PC",CN)
colnames(ukbc)<-CN
## Age
ukbc$Age<-ukbc$age_attended_assessment_center
ukbc$Age_sq<-ukbc$Age*ukbc$Age

## Caucasian
ukbc$Caucasian<-rep(1,nrow(ukbc))
ukbc[is.na(ukbc$genetic_ethnic_grouping),]$Caucasian<-0

## Factorize the variables
ukbc$Caucasian<-as.factor(ukbc$Caucasian)


ukbc$Age_categ<-rep("40-49",nrow(ukbc))
ukbc[ukbc$age_attended_assessment_center >= 50 & ukbc$age_attended_assessment_center < 60,]$Age_categ<-"50-59"
ukbc[ukbc$age_attended_assessment_center >= 60,]$Age_categ<-">60"
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("40-49","50-59",">60"))


ukbc<-ukbc[!is.na(ukbc$EverSmoke),]

## CH group
ukbc$Mmca<-rep(0,nrow(ukbc))
ukbc[ukbc$MCNV == 1 & ukbc$LCNV == 0,]$Mmca<-1
ukbc$Lmca<-rep(0,nrow(ukbc))
ukbc[ukbc$MCNV == 0 & ukbc$LCNV == 1,]$Lmca<-1
ukbc$Amca<-rep(0,nrow(ukbc))
ukbc[ukbc$MCNV == 1 & ukbc$LCNV == 1,]$Amca<-1
ukbc$Umca<-rep(0,nrow(ukbc))
ukbc[ukbc$MCNV == 0 & ukbc$LCNV == 0 & ukbc$CNV == 1,]$Umca<-1

ukbc$Mmca<-as.factor(ukbc$Mmca)
ukbc$Lmca<-as.factor(ukbc$Lmca)
ukbc$Amca<-as.factor(ukbc$Amca)
ukbc$Umca<-as.factor(ukbc$Umca)


#=====================================================================
## Enrichment of specific lymphoid malignancies

DOD<-c()
death2<-c()
forcedCensorDate<-"2020-03-31"
deathCensorDate<-"2020-03-31"


icdList<-as.character(hicd[hicd$Group %in% c("Lymphoid","Plasma_cell"),]$ICD)
CLL<-c("C911","C830")
PC<-c("C900","D472")
HL<-c("C81",paste("C81",0:9,sep=""))
FL<-c("C82",paste("C82",0:9,sep=""))
DLBCL<-c("C833")
NHL<-c("C85",paste("C85",0:9,sep=""))
WS<-c("C880")
CIRC<-c("C901","C910","C913","C914","C831","C841","C915","C916")
nhl<-icdList[!icdList %in% c(PC,CLL,HL,FL,DLBCL,NHL,WS,CIRC)]

## Myeloid
icdList<-as.character(hicd[hicd$Group %in% c("Myeloid"),]$ICD)
mds_icdList<-c(paste("D46",0:9,sep=""),"C946")
aml_icdList<-paste("C92",0,sep="")
mpn_icdList<-c("D752","D45","C945","D471","D473")
other_myeloid<-icdList[!icdList %in% c(mds_icdList,aml_icdList,mpn_icdList)]


## Prepare status and follow-up date in the function
getStatusYear<-function(ukbc,allICD,icdlist,EXCL=""){
	## Lymphoid malignancies
	#tdf1<-getDiag(ukbc,c(icdlist,EXCL),any10,date10)
	#tdf2<-getCancer(ukbc,c(icdlist,EXCL),cancer10,cancerDate10)
	tdf1<-getDiag(ukbc,allICD,any10,date10)
	tdf2<-getCancer(ukbc,allICD,cancer10,cancerDate10)
	tempdf<-rbind(tdf1,tdf2)
	#exclEid<-tempdf[tempdf$ICD %in% EXCL,]$eid
	#tempdf<-tempdf[!tempdf$eid %in% exclEid,]
	tempdf$Diff2<-as.numeric(difftime(forcedCensorDate,tempdf$Date,units=c("days")))
	tempdf<-tempdf[tempdf$Diff2 >= 0,]
	tdf<-getUnique(tempdf[,colnames(tdf1)])
	tdf<-tdf[tdf$ICD %in% icdlist,]
	## Positive cases
	cenDate<-tdf$Date
	names(cenDate)<-as.character(tdf$eid)
	## dead individuals
	t1<-ukbc[(!ukbc$eid %in% tdf$eid) & ukbc$dead == 1,c("eid","death_censor_date")]
	cenDate1<-t1$death_censor_date
	names(cenDate1)<-as.character(t1$eid)
	## individuals who did not die
	t2<-ukbc[(!ukbc$eid %in% tdf$eid) & ukbc$dead == 0,]$eid
	cenDate2<-rep(forcedCensorDate,length(t2))
	names(cenDate2)<-as.character(t2)
	## Has disease
	Status<-c(rep(1,length(cenDate)),rep(0,length(cenDate1)+length(cenDate2)))
	names(Status)<-c(names(cenDate),names(cenDate1),names(cenDate2))
	DATE<-c(as.character(cenDate),as.character(cenDate1),cenDate2)
	return (list(Status,DATE))
}


#=====================================================================
## Cox model
MYL<-list(CLL,PC,HL,FL,DLBCL,NHL,WS,CIRC,nhl,mds_icdList,aml_icdList,mpn_icdList,other_myeloid)
names(MYL)<-c("CLL","MM","HL","FL","DLBCL","NHL","WS","CIRC","Other_lymphoid","MDS","AML","MPN","Other_myeloid")


allICD<-c(CLL,PC,HL,FL,DLBCL,NHL,WS,CIRC,nhl,mds_icdList,aml_icdList,mpn_icdList,other_myeloid)


K<-1
for (MALIG in names(MYL)){
	print (MALIG)
	dlist<-MYL[[MALIG]]
	ukbc$Status<-rep(0,nrow(ukbc))
	ukbc$Date<-rep(forcedCensorDate,nrow(ukbc))
	ukbc$Year<-NULL
	if (MALIG %in% c("MM","HL","FL","DLBCL","NHL","WS","CIRC","Other_lymphoid")){
		SD<-getStatusYear(ukbc,allICD,dlist,CLL)
	}else{
		SD<-getStatusYear(ukbc,allICD,dlist,"")
	}
	ukbc$Status<-SD[[1]][as.character(ukbc$eid)]
	ukbc$Date<-SD[[2]][as.character(ukbc$eid)]
	ukbc$Year<-as.numeric(difftime(ukbc$Date,ukbc$date_attending_assessment_center,units=c("days")))/365.25
	q95<-as.numeric(quantile(ukbc$Year,0.95))
	ukbc[ukbc$Year > q95,]$Status<-0
	ukbc[ukbc$Year > q95,]$Year<-q95
	## Now do a cox model for L-mCA
	NE<-nrow(ukbc[ukbc$Lmca == 1 & ukbc$Status == 1,])
	if (NE >= 4){
		cox1<-coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + Lmca, data = ukbc[ukbc$Lmca == 1 | ukbc$CNV == 0,])
		x1<-summary(cox1)
		mvdf1<-as.data.frame(x1$conf.int)
		colnames(mvdf1)<-c("Haz_ratio","negHR","CI_low","CI_high")
		mvdf1$P<-x1$coefficients[,5]
		mvdf1$Variable<-rownames(mvdf1)
		##
		Variable<-"Lmca1"
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "L-mCA",
			Total = nrow(ukbc[ukbc$Lmca == 1,]),
			Event = nrow(ukbc[ukbc$Lmca == 1 & ukbc$Status == 1,]),
			HR = mvdf1[mvdf1$Variable %in% Variable,]$Haz_ratio,
			CI_low = mvdf1[mvdf1$Variable %in% Variable,]$CI_low,
			CI_high = mvdf1[mvdf1$Variable %in% Variable,]$CI_high,
			P = mvdf1[mvdf1$Variable %in% Variable,]$P)
	}else{
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "L-mCA",
			Total = nrow(ukbc[ukbc$Lmca == 1,]),
			Event = nrow(ukbc[ukbc$Lmca == 1 & ukbc$Status == 1,]),
			HR = NA,
			CI_low = NA,
			CI_high = NA,
			P = NA)
	}
	if (K == 1){
		mydf<-tempdf
		K<-2
	}else{
		mydf<-rbind(mydf,tempdf)
	}
	## Now do a cox model for M-mCA
	## Now do a cox model for M-mCA
	NE<-nrow(ukbc[ukbc$Mmca == 1 & ukbc$Status == 1,])
	if (NE >= 4){
		cox1<-coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + Mmca, data = ukbc[ukbc$Mmca == 1 | ukbc$CNV == 0,])
		x1<-summary(cox1)
		mvdf1<-as.data.frame(x1$conf.int)
		colnames(mvdf1)<-c("Haz_ratio","negHR","CI_low","CI_high")
		mvdf1$P<-x1$coefficients[,5]
		mvdf1$Variable<-rownames(mvdf1)
		##
		Variable<-"Mmca1"
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "M-mCA",
			Total = nrow(ukbc[ukbc$Mmca == 1,]),
			Event = nrow(ukbc[ukbc$Mmca == 1 & ukbc$Status == 1,]),
			HR = mvdf1[mvdf1$Variable %in% Variable,]$Haz_ratio,
			CI_low = mvdf1[mvdf1$Variable %in% Variable,]$CI_low,
			CI_high = mvdf1[mvdf1$Variable %in% Variable,]$CI_high,
			P = mvdf1[mvdf1$Variable %in% Variable,]$P)
	}else{
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "M-mCA",
			Total = nrow(ukbc[ukbc$Mmca == 1,]),
			Event = nrow(ukbc[ukbc$Mmca == 1 & ukbc$Status == 1,]),
			HR = NA,
			CI_low = NA,
			CI_high = NA,
			P = NA)
	}
	if (K == 1){
		mydf<-tempdf
		K<-2
	}else{
		mydf<-rbind(mydf,tempdf)
	}
	## Now do a cox model for A-mCA
	NE<-nrow(ukbc[ukbc$Amca == 1 & ukbc$Status == 1,])
	if (NE >= 4){
		cox1<-coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + Amca, data = ukbc[ukbc$Amca == 1 | ukbc$CNV == 0,])
		x1<-summary(cox1)
		mvdf1<-as.data.frame(x1$conf.int)
		colnames(mvdf1)<-c("Haz_ratio","negHR","CI_low","CI_high")
		mvdf1$P<-x1$coefficients[,5]
		mvdf1$Variable<-rownames(mvdf1)
		##
		Variable<-"Amca1"
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "A-mCA",
			Total = nrow(ukbc[ukbc$Amca == 1,]),
			Event = nrow(ukbc[ukbc$Amca == 1 & ukbc$Status == 1,]),
			HR = mvdf1[mvdf1$Variable %in% Variable,]$Haz_ratio,
			CI_low = mvdf1[mvdf1$Variable %in% Variable,]$CI_low,
			CI_high = mvdf1[mvdf1$Variable %in% Variable,]$CI_high,
			P = mvdf1[mvdf1$Variable %in% Variable,]$P)
	}else{
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "A-mCA",
			Total = nrow(ukbc[ukbc$Amca == 1,]),
			Event = nrow(ukbc[ukbc$Amca == 1 & ukbc$Status == 1,]),
			HR = NA,
			CI_low = NA,
			CI_high = NA,
			P = NA)
	}
	if (K == 1){
		mydf<-tempdf
		K<-2
	}else{
		mydf<-rbind(mydf,tempdf)
	}
	## Now do a cox model for M-mCA
	NE<-nrow(ukbc[ukbc$Umca == 1 & ukbc$Status == 1,])
	if (NE >= 4){
		cox1<-coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + Umca, data = ukbc[ukbc$Umca == 1 | ukbc$CNV == 0,])
		x1<-summary(cox1)
		mvdf1<-as.data.frame(x1$conf.int)
		colnames(mvdf1)<-c("Haz_ratio","negHR","CI_low","CI_high")
		mvdf1$P<-x1$coefficients[,5]
		mvdf1$Variable<-rownames(mvdf1)
		##
		Variable<-"Umca1"
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "U-mCA",
			Total = nrow(ukbc[ukbc$Umca == 1,]),
			Event = nrow(ukbc[ukbc$Umca == 1 & ukbc$Status == 1,]),
			HR = mvdf1[mvdf1$Variable %in% Variable,]$Haz_ratio,
			CI_low = mvdf1[mvdf1$Variable %in% Variable,]$CI_low,
			CI_high = mvdf1[mvdf1$Variable %in% Variable,]$CI_high,
			P = mvdf1[mvdf1$Variable %in% Variable,]$P)
	}else{
		tempdf<-data.frame(Malignancy = MALIG,
			CHIP = "U-mCA",
			Total = nrow(ukbc[ukbc$Umca == 1,]),
			Event = nrow(ukbc[ukbc$Umca == 1 & ukbc$Status == 1,]),
			HR = NA,
			CI_low = NA,
			CI_high = NA,
			P = NA)
	}
	if (K == 1){
		mydf<-tempdf
		K<-2
	}else{
		mydf<-rbind(mydf,tempdf)
	}
}



write.table(mydf,file=paste(figDir,"scripts/data/UKB500_malignancy_subtypes_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

}

#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
## Generating the plots

if (CMD == "FIG"){


library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)

ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"


#=====================================================================
## Histogram Lymphoid

tdf1<-data.frame(fread(paste(figDir,"scripts/data/WES_malignancy_subtypes_df.txt",sep="")))
tdf2<-data.frame(fread(paste(figDir,"scripts/data/UKB500_malignancy_subtypes_df.txt",sep="")))

tdf<-rbind(tdf1,tdf2)


tdf$CHIP <-factor(tdf$CHIP,levels=c("U-mCA","A-mCA","L-mCA","M-mCA","L-CHIP","M-CHIP"))


tdf$Malignancy<-as.character(tdf$Malignancy)
tdf$Malignancy<-gsub("CIRC","Circulating\nlymphoma",tdf$Malignancy)
tdf$Malignancy<-gsub("MM","MM/MGUS",tdf$Malignancy)
tdf$Malignancy<-gsub("_","\n",tdf$Malignancy)

tdf$Malignancy <-factor(tdf$Malignancy,levels=rev(c("Other\nlymphoid","HL","WS","FL","DLBCL","NHL","Circulating\nlymphoma","CLL","MM/MGUS","Other\nmyeloid","AML","MDS","MPN")))

p1 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8,shape=15,fatten=3, position=position_dodge(width=0.9)) +
	facet_grid(rows = vars(Malignancy)) +
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid.major.x=element_line(),
		axis.line.y=element_blank(),
		axis.text.x=element_text(color="black"),
		axis.text.y=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none",
		strip.text.y = element_text(angle = 0),
		panel.border = element_rect(colour = "black", fill=NA))+
	scale_color_manual(values=c("grey","orange","blue","red","blue","red"))+
	scale_y_log10()+
	coord_flip()


tdf<-tdf[order(tdf$Malignancy,tdf$CHIP),]

tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("CHIP","Malignancy","Total","Event","HR_ci","P")]
colnames(tdf1)<-c("CHIP","Malignancy","Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))

F1<-ggarrange(p1,ptdf,ncol=2,nrow=1,widths=c(1.5,2))


pdf(paste(figDir,"scripts/Final_figs/Fig_ED3.pdf",sep=""),width=10,height=18)
print (F1)
dev.off()


}
