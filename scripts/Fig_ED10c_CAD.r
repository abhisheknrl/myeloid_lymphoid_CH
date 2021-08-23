
## CVD analysis

## WES, 500, FIG


args = commandArgs(trailingOnly=TRUE)

CMD<-args[1]


if (CMD == "WES"){

library(ggplot2)
library(data.table)
library(ggpubr)


setwd("")
hemeICDFile<-"meta/heme_malignancies_icd10.txt"

## CAD data from Mike
cadFile<-"CAD/202006.Coronary_Artery_Disease_INTERMEDIATE.tab.tsv"
idMapFile<-"ukb_eid_map_Ebert_Natarajan.txt"

## BMI
bmiFile<-"data/body_measurement.tab"

## Updated death data
deathFile<-"meta/death.txt"
codFile<-"meta/death_cause.txt"


figDir<-""



##===================================================================

##ICD10codes
hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")

##===================================================================
## Mapping the ids
mf<-read.table(idMapFile,head=TRUE,sep="\t")
EID<-mf$eid
names(EID)<-as.character(mf$Natarajan_eid)

#=====================================================================
## load UKB File

pn<-data.frame(fread(cadFile))
pn<-pn[pn$sample_id %in% names(EID),]
pn$eid<-EID[as.character(pn$sample_id)]
pn$CAD_year<-round(as.numeric(difftime(pn$censor_date,pn$enroll_date,units=c("days")))/365.25,3)

pn1<-pn[,c("eid","has_disease","CAD_year")]
colnames(pn1)<-c("eid","CAD","CAD_year")


#=====================================================================
## Reload data

ukbc<-data.frame(fread(paste(figDir,"1_UKB_data.txt",sep="")))
ukbc$Caucasian<-rep(1,nrow(ukbc))
ukbc[is.na(ukbc$genetic_ethnic_grouping),]$Caucasian<-0

## Adding covariates
ukbc$EverSmoke<-ukbc$ever_smoked
## PCA
CN<-colnames(ukbc)
CN<-gsub("genetic_pca","PC",CN)
colnames(ukbc)<-CN

## Age
ukbc$Age<-ukbc$age_attended_assessment_center
ukbc$Age_sq<-ukbc$Age*ukbc$Age


#=====================================================================

## Death data
death<-read.table(deathFile,head=TRUE,sep="\t")
death<-death[death$eid %in% ukbc$eid,]
DOD<-as.character(death$date_of_death)
DOD<-as.Date(DOD,"%d/%m/%Y")
names(DOD)<-as.character(death$eid)

death2<-read.table(codFile,head=TRUE,sep="\t")
death2<-death2[death2$eid %in% ukbc$eid,]

## Adding the date of death
deathCensorDate<-"2020-03-31"
ukbc$dead<-rep(0,nrow(ukbc))
ukbc[ukbc$eid %in% names(DOD),]$dead<-1
ukbc$death_censor_date<-DOD[as.character(ukbc$eid)]
ukbc[is.na(ukbc$death_censor_date),]$death_censor_date<-deathCensorDate

## Censoring individuals who died after the censorDate
ukbc$Diff<-as.numeric(difftime(deathCensorDate,ukbc$death_censor_date,units=c("days")))
ukbc[ukbc$Diff < 0,]$dead<-0
ukbc[ukbc$Diff < 0,]$death_censor_date<-deathCensorDate

ukbc$death_year<-as.numeric(difftime(ukbc$death_censor_date,ukbc$date_attending_assessment_center,units=c("days")))/365.25



#===========================================================================
#===========================================================================

forcedCensorDate<-"2020-03-31"

tukb<-ukbc

## CAD from Pradeep
tukb<-merge(tukb,pn1)

## Renaming CAD as CVD to map the rest of the scriptt
tukb$CVD<-tukb$CAD
tukb$CVD_year<-tukb$CAD_year

#========================================================================
## Hypertension status
ht_icd9<-c(401, 4010, 4011, 4019,402, 4020, 4021, 4029, 403,4030, 4031, 4039, 404, 4040,4041, 4049, 405, 4050, 4051,4059)
ht_icd10<-c("I10","I11","I110","I119","I12","I120","I129","I13","I130","I131","I132","I139","I15","I150","I151","I152","I158","I159")
t2d_icd9<-c(2500, 25000, 25001, 25009, 25011, 25019, 2503, 2504, 2505, 25099)
t2d_icd10<-c("E10","E101","E102","E103","E104","E105","E106","E107","E108","E109","E11","E110","E111","E112","E113","E114","E115","E116","E117","E118","E119","E12","E121","E128","E129","E13","E131","E132","E133","E135","E136","E137","E138","E139","E14","E141","E142","E143","E144","E145","E146","E147","E148","E149")


any10<-paste("f.41270.0.",0:212,sep="")
any9<-paste("f.41271.0.",0:46,sep="")
date10<-paste("f.41280.0.",0:212,sep="")
date9<-paste("f.41281.0.",0:46,sep="")


## Function to extract samples with diagnosis
getDiagStatus<-function(ukb,icdlist,columnNames){
	eidList<-c()
	for (CN in columnNames){
		#print (CN)
		temp<-ukb[!is.na(ukb[,c(CN)]),]
		if (nrow(temp) == 0){
			next
		}
		temp<-temp[temp[,c(CN)] %in% icdlist,]
		if (nrow(temp)==0){
			next
		}
		eidList<-c(eidList,temp$eid)
	}
	return (unique(eidList))
}

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


## H2T
diag1<-getDiag(tukb,ht_icd9,any9,date9)
diag2<-getDiag(tukb,ht_icd10,any10,date10)
ht<-rbind(diag1,diag2)
ht$Diff<-as.numeric(difftime(ht$Date,ht$DOAAC,units=c("days")))
ht<-ht[ht$Diff < 0,]

## Add column
tukb$HT<-rep(0,nrow(tukb))
tukb[tukb$eid %in% ht$eid,]$HT<-1

## T2D
t2d1<-getDiag(tukb,t2d_icd9,any9,date9)
t2d2<-getDiag(tukb,t2d_icd10,any10,date10)
t2d<-rbind(t2d1,t2d2)
t2d$Diff<-as.numeric(difftime(t2d$Date,t2d$DOAAC,units=c("days")))
t2d<-t2d[t2d$Diff < 0,]

## Add column
tukb$T2D<-rep(0,nrow(tukb))
tukb[tukb$eid %in% t2d$eid,]$T2D<-1

## Statin users
## simvastatin, fluvastatin, pravastatin,	atorvastatin,	rosuvastatin

codeList<-c(1140861958,1140888594,1140888648,1141146234,1141192410)
colList<-paste("f.20003.0.",0:47,sep="")
Statin<-getDiagStatus(tukb,codeList,colList)

## Add column
tukb$Statin<-rep(0,nrow(tukb))
tukb[tukb$eid %in% Statin,]$Statin<-1

## HDL & LDL
tukb$HDL<- tukb$hdl_cholesterol
tukb$LDL <- tukb$ldl_direct

tmp1<-tukb[tukb$Statin == 0,]
tmp2<-tukb[tukb$Statin == 1,]
tmp2$LDL_c<-tmp2$LDL/0.68
tmp1$LDL_c<-tmp1$LDL
tukb<-rbind(tmp1,tmp2)

## BMI

bmi<-read.table(bmiFile,head=TRUE,sep="\t")
BMI<-bmi$f.21001.0.0
names(BMI)<-as.character(bmi$f.eid)

tukb$BMI<-BMI[as.character(tukb$eid)]

CN1<-c("eid","Age","Age_sq","Sex","EverSmoke","Caucasian",paste("PC",1:5,sep=""),"HT","T2D","LDL_c","HDL","BMI")
CN2<-c("dead","death_year","MI","MI_year","ST","ST_year","CAD","CAD_year","CVD","CVD_year")
CN3<-c("MCHIP","Msize","CHIP","CHIP_size","CNV","CNV_size")

write.table(tukb[,c(CN1,CN2,CN3)],file=paste(figDir,"scripts/data/CVD_WES_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


}

#================================================================================
#================================================================================
#================================================================================
#================================================================================
#================================================================================
#================================================================================
#================================================================================
#================================================================================
## Copy number alteration data

if (CMD == 500){

library(ggplot2)
library(data.table)
library(ggpubr)


setwd("")
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"

## Updated death data
deathFile<-"meta/death.txt"
codFile<-"meta/death_cause.txt"

## Related individuals
relatedFile<-"meta/related_individuals.dat"

## CAD file from Mike
cadFile<-"CAD/Coronary_Artery_Disease.tab.tsv"
idMapFile<-"ukb_eid_map_7089_50834_s487K.txt"

## BMI
bmiFile<-"data/body_measurement.tab"


##===================================================================
## UKB

ukb<-data.frame(fread(paste(figDir,"data/UKB_500K.txt",sep="")))

ukb$Caucasian<-rep(1,nrow(ukb))
ukb[is.na(ukb$genetic_ethnic_grouping),]$Caucasian<-0


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


## Censoring individuals who died after the censorDate
ukb$Diff<-as.numeric(difftime(deathCensorDate,ukb$death_censor_date,units=c("days")))
ukb[ukb$Diff < 0,]$dead<-0
ukb[ukb$Diff < 0,]$death_censor_date<-deathCensorDate

ukb$death_year<-as.numeric(difftime(ukb$death_censor_date,ukb$date_attending_assessment_center,units=c("days")))/365.25


##===================================================================
## Mapping the ids
mf<-read.table(idMapFile,head=TRUE,sep="\t")
mf<-mf[mf$ID_1_50834 %in% ukb$eid,]
EID<-mf$ID_1_50834
names(EID)<-as.character(mf$ID_1_7089)

#=====================================================================
## load UKB File

pn<-data.frame(fread(cadFile))
pn<-pn[pn$sample_id %in% names(EID),]
pn$eid<-EID[as.character(pn$sample_id)]
pn$CAD_year<-round(as.numeric(difftime(pn$censor_date,pn$enroll_date,units=c("days")))/365.25,3)

pn1<-pn[,c("eid","has_disease","CAD_year")]
colnames(pn1)<-c("eid","CAD","CAD_year")

ukb<-ukb[ukb$eid %in% pn1$eid,]
## CAD from Pradeep
tukb<-merge(ukb,pn1)


#=====================================================================
##

forcedCensorDate<-"2020-03-31"

## Renaming CAD as CVD
tukb$CVD<-tukb$CAD
tukb$CVD_year<-tukb$CAD_year


#========================================================================
## Hypertension status
forcedCensorDate<-"2020-03-31"


ht_icd9<-c(401, 4010, 4011, 4019,402, 4020, 4021, 4029, 403,4030, 4031, 4039, 404, 4040,4041, 4049, 405, 4050, 4051,4059)
ht_icd10<-c("I10","I11","I110","I119","I12","I120","I129","I13","I130","I131","I132","I139","I15","I150","I151","I152","I158","I159")
t2d_icd9<-c(2500, 25000, 25001, 25009, 25011, 25019, 2503, 2504, 2505, 25099)
t2d_icd10<-c("E10","E101","E102","E103","E104","E105","E106","E107","E108","E109","E11","E110","E111","E112","E113","E114","E115","E116","E117","E118","E119","E12","E121","E128","E129","E13","E131","E132","E133","E135","E136","E137","E138","E139","E14","E141","E142","E143","E144","E145","E146","E147","E148","E149")


any10<-paste("f.41270.0.",0:212,sep="")
any9<-paste("f.41271.0.",0:46,sep="")
date10<-paste("f.41280.0.",0:212,sep="")
date9<-paste("f.41281.0.",0:46,sep="")


## Function to extract samples with diagnosis
getDiagStatus<-function(ukb,icdlist,columnNames){
	eidList<-c()
	for (CN in columnNames){
		#print (CN)
		temp<-ukb[!is.na(ukb[,c(CN)]),]
		if (nrow(temp) == 0){
			next
		}
		temp<-temp[temp[,c(CN)] %in% icdlist,]
		if (nrow(temp)==0){
			next
		}
		eidList<-c(eidList,temp$eid)
	}
	return (unique(eidList))
}

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


## H2T
diag1<-getDiag(tukb,ht_icd9,any9,date9)
diag2<-getDiag(tukb,ht_icd10,any10,date10)
ht<-rbind(diag1,diag2)
ht$Diff<-as.numeric(difftime(ht$Date,ht$DOAAC,units=c("days")))
ht<-ht[ht$Diff < 0,]

## Add column
tukb$HT<-rep(0,nrow(tukb))
tukb[tukb$eid %in% ht$eid,]$HT<-1

## T2D
t2d1<-getDiag(tukb,t2d_icd9,any9,date9)
t2d2<-getDiag(tukb,t2d_icd10,any10,date10)
t2d<-rbind(t2d1,t2d2)
t2d$Diff<-as.numeric(difftime(t2d$Date,t2d$DOAAC,units=c("days")))
t2d<-t2d[t2d$Diff < 0,]

## Add column
tukb$T2D<-rep(0,nrow(tukb))
tukb[tukb$eid %in% t2d$eid,]$T2D<-1

## Statin users
## simvastatin, fluvastatin, pravastatin,	atorvastatin,	rosuvastatin

codeList<-c(1140861958,1140888594,1140888648,1141146234,1141192410)
colList<-paste("f.20003.0.",0:47,sep="")
Statin<-getDiagStatus(tukb,codeList,colList)

## Add column
tukb$Statin<-rep(0,nrow(tukb))
tukb[tukb$eid %in% Statin,]$Statin<-1

## HDL & LDL
tukb$HDL<- tukb$hdl_cholesterol
tukb$LDL <- tukb$ldl_direct

tmp1<-tukb[tukb$Statin == 0,]
tmp2<-tukb[tukb$Statin == 1,]
tmp2$LDL_c<-tmp2$LDL/0.68
tmp1$LDL_c<-tmp1$LDL
tukb<-rbind(tmp1,tmp2)

## BMI

bmi<-read.table(bmiFile,head=TRUE,sep="\t")
BMI<-bmi$f.21001.0.0
names(BMI)<-as.character(bmi$f.eid)

tukb$BMI<-BMI[as.character(tukb$eid)]



CN1<-c("eid","Age","Age_sq","Sex","EverSmoke","Caucasian",paste("PC",1:5,sep=""),"HT","T2D","LDL_c","HDL","BMI")
CN2<-c("dead","death_year","MI","MI_year","ST","ST_year","CAD","CAD_year","CVD","CVD_year")
CN3<-c("LCNV","LCNV_size","MCNV","MCNV_size","CNV","CNV_size")

write.table(tukb[,c(CN1,CN2,CN3)],file=paste(figDir,"scripts/data/CVD_500K_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
}

#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================
## GEnerate figures 
## WES
## Status and death competition

if (CMD == "FIG"){


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


## Function to get forest data
extractCox<-function(cox1,Disease,Variable,CHIP,Group){
	x1<-summary(cox1)
	mvdf1<-as.data.frame(x1$conf.int)
	colnames(mvdf1)<-c("Haz_ratio","negHR","CI_low","CI_high")
	mvdf1$P<-x1$coefficients[,5]
	mvdf1$Variable<-rownames(mvdf1)
	##
	df1<-data.frame(Malignancy = Disease,
		CHIP = CHIP,
		HR = mvdf1[mvdf1$Variable %in% Variable,]$Haz_ratio,
		CI_low = mvdf1[mvdf1$Variable %in% Variable,]$CI_low,
		CI_high = mvdf1[mvdf1$Variable %in% Variable,]$CI_high,
		P = mvdf1[mvdf1$Variable %in% Variable,]$P,
		VAF = Group)
	return (df1)
}


#===========================================================================
## 

tukb<-data.frame(fread(paste(figDir,"scripts/data/CVD_WES_df.txt",sep="")))


#tukb1<-tukb[tukb$Caucasian == 1,]
tukb1<-tukb

tukb1$MCHIP<-as.factor(tukb1$MCHIP)
tukb1$Msize<-as.factor(tukb1$Msize)
tukb1$CHIP<-as.factor(tukb1$CHIP)
tukb1$CHIP_size<-as.factor(tukb1$CHIP_size)
tukb1$Caucasian<-as.factor(tukb1$Caucasian)
tukb1$HT<-as.factor(tukb1$HT)
tukb1$T2D<-as.factor(tukb1$T2D)

tukb1$Status<-tukb1$CVD
tukb1$Year<-tukb1$CVD_year

## Death
temp1<-tukb1[tukb1$death_year < tukb1$Year & tukb1$Status == 0,]
temp2<-tukb1[!tukb1$eid %in% temp1$eid,]
temp1$Year<-temp1$death_year

tukb1<-rbind(temp1,temp2)


## Exclude prevalent cases
tukb1<-tukb1[tukb1$Year > 0,]

## q95
q95<-as.numeric(quantile(tukb1$Year,0.95))
tukb1[tukb1$Year >= q95,]$Status<-0
tukb1[tukb1$Year >= q95,]$Year<-q95


#===============================================================
tukb1$CH<-rep("Control",nrow(tukb1))
tukb1[tukb1$CHIP == 1 & tukb1$MCHIP == 1,]$CH<-"ML"
tukb1[tukb1$CHIP == 0 & tukb1$MCHIP == 1,]$CH<-"CHMD"
tukb1[tukb1$CHIP == 1 & tukb1$MCHIP == 0,]$CH<-"CHLD"

tukb1<-tukb1[tukb1$CH != "ML",]

tukb1$CH<-factor(tukb1$CH,levels=c("Control","CHMD","CHLD"))

## Covariats missing
tukb1<-tukb1[!is.na(tukb1$EverSmoke) & !is.na(tukb1$BMI),]
tukb1$Age_categ<-rep("<50",nrow(tukb1))
tukb1[tukb1$Age >= 50 & tukb1$Age < 60,]$Age_categ<-"50-59"
tukb1[tukb1$Age >= 60,]$Age_categ<-">60"

## Analyzing separately
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + HT + T2D + BMI + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1)
tukb2<-tukb1[tukb1$Msize != 1 & tukb1$CHIP_size != 1,]
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + HT + T2D + BMI + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb2)



mvdf1<-extractCox(cox1,"CVD","CHCHMD","CHMD","CHMD")
mvdf2<-extractCox(cox2,"CVD","CHCHMD","CHMD-high","CHMD")
mvdf3<-extractCox(cox1,"CVD","CHCHLD","CHLD","CHLD")
mvdf4<-extractCox(cox2,"CVD","CHCHLD","CHLD-high","CHLD")
BL<-data.frame(Malignancy = "", CHIP = "Baseline", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Control")
temp3<-rbind(BL,mvdf1,mvdf3,mvdf2,mvdf4)
temp3$at_risk<-c(as.numeric(table(tukb1$CH)[c("Control","CHMD","CHLD")]),as.numeric(table(tukb2$CH)[c("CHMD","CHLD")]))
temp3$event<-c(as.numeric(table(tukb1[tukb1$Status == 1,]$CH)[c("Control","CHMD","CHLD")]),as.numeric(table(tukb2[tukb2$Status == 1,]$CH)[c("CHMD","CHLD")]))
rownames(temp3)<-as.character(temp3$CHIP)
temp3<-temp3[c("Baseline","CHMD","CHMD-high","CHLD","CHLD-high"),]

## Table
tdf<-temp3
tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("VAF","at_risk","event","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf1<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 11,
		padding = unit(c(2, 2), "mm")))

temp3$CHIP<-factor(temp3$CHIP,levels=c("CHLD-high","CHLD","CHMD-high","CHMD","Baseline"))
p1 <- ggplot(data=temp3, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	#facet_grid(CHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[1],2),rep(ctList[2],2),"black"))

Fa<-ggarrange(ptdf1,p1,ncol=2,nrow=1,widths=c(2,1))

mydf<-temp3

print ("CH done")


#===========================================================================
#===========================================================================
#===========================================================================
## For mCA

tukb<-data.frame(fread(paste(figDir,"scripts/data/CVD_500K_df.txt",sep="")))

cnv<-data.frame(fread(paste(figDir,"data/UKB_500K_CNV.txt",sep="")))
cnv[is.na(cnv$CELL_FRAC),]$CELL_FRAC<-0.001
cnv<-cnv[cnv$CHR %in% 1:22,]

## mCA groups
tukb$CH<-rep("Control",nrow(tukb))
tukb[tukb$MCNV == 0 & tukb$LCNV == 1,]$CH<-"mCALD"
tukb[tukb$MCNV == 1 & tukb$LCNV == 0,]$CH<-"mCAMD"
tukb[tukb$MCNV == 1 & tukb$LCNV == 1,]$CH<-"mCAAD"
tukb[tukb$MCNV == 0 & tukb$LCNV == 0 & tukb$CNV == 1,]$CH<-"mCAUD"
tukb$CH<-factor(tukb$CH,levels=c("Control","mCAMD","mCALD","mCAAD","mCAUD"))

tukb$MCNV<-as.factor(tukb$MCNV)
tukb$LCNV<-as.factor(tukb$LCNV)
tukb$MCNV_size<-as.factor(tukb$MCNV_size)
tukb$LCNV_size<-as.factor(tukb$LCNV_size)


## For the analysis
tukb1<-tukb

tukb1$Status<-tukb1$CVD
tukb1$Year<-tukb1$CVD_year

## Death
temp1<-tukb1[tukb1$death_year < tukb1$Year & tukb1$Status == 0,]
temp2<-tukb1[!tukb1$eid %in% temp1$eid,]
temp1$Year<-temp1$death_year

tukb1<-rbind(temp1,temp2)


## Exclude prevalent cases
tukb1<-tukb1[tukb1$Year > 0,]

## q95
q95<-as.numeric(quantile(tukb1$Year,0.95))
tukb1[tukb1$Year >= q95,]$Status<-0
tukb1[tukb1$Year >= q95,]$Year<-q95


## Missing covariates
tukb1<-tukb1[!is.na(tukb1$EverSmoke) & !is.na(tukb1$BMI),]

## Categorical age
tukb1$Age_categ<-rep("<50",nrow(tukb1))
tukb1[tukb1$Age >= 50 & tukb1$Age < 60,]$Age_categ<-"50-59"
tukb1[tukb1$Age >= 60,]$Age_categ<-">60"
tukb1$Age_categ<-factor(tukb1$Age_categ,levels=c("<50","50-59",">60"))

## COX
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + HT + EverSmoke + T2D + BMI + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1)
cox4 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + HT + EverSmoke + T2D + BMI + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1[tukb1$CNV_size != 1,])

mvdf1<-extractCox(cox2,"Lymphoid",c("CHmCAMD","CHmCALD"),c("MmCA","LmCA"),"Any")
mvdf5<-extractCox(cox2,"Lymphoid",c("CHmCAAD","CHmCAUD"),c("AmCA","UmCA"),"Any")
mvdf2<-extractCox(cox4,"Lymphoid",c("CHmCAMD","CHmCALD"),c("MmCA-high","LmCA-high"),"Any")
mvdf4<-extractCox(cox4,"Lymphoid",c("CHmCAAD","CHmCAUD"),c("AmCA-high","UmCA-high"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "No mCA", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5,mvdf2,mvdf4)
tdf$at_risk<-c(as.numeric(table(tukb1$CH)[c("Control","mCAMD","mCALD","mCAAD","mCAUD")]),as.numeric(table(tukb1[tukb1$CNV_size != 1,]$CH)[c("mCAMD","mCALD","mCAAD","mCAUD")]))
tdf$event<-c(as.numeric(table(tukb1[tukb1$Status ==1,]$CH)[c("Control","mCAMD","mCALD","mCAAD","mCAUD")]),as.numeric(table(tukb1[tukb1$CNV_size != 1 & tukb1$Status ==1,]$CH)[c("mCAMD","mCALD","mCAAD","mCAUD")]))
rownames(tdf)<-as.character(tdf$CHIP)
tdf<-tdf[c("No mCA","MmCA","MmCA-high","LmCA","LmCA-high","AmCA","AmCA-high","UmCA","UmCA-high"),]


tdf$CHIP<-factor(tdf$CHIP,levels=c("UmCA-high","UmCA","AmCA-high","AmCA","LmCA-high","LmCA","MmCA-high","MmCA","No mCA"))
p2 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	#facet_grid(CHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep("grey",2),rep("orange",2),rep(ctList[1],2),rep(ctList[2],2),"black"))



temp3<-tdf
## Tables
tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("CHIP","at_risk","event","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf2<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 11,
		padding = unit(c(2, 2), "mm")))


Fb<-ggarrange(ptdf2,p2,ncol=2,nrow=1,widths=c(2,1))


#===========================================================================
## Combined forest plots for somatic variants and mCAs

tdf<-rbind(mydf,temp3)


tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)

print (head(tdf))
print (tdf$CHIP)

tdf1<-tdf[,c("CHIP","at_risk","event","HR_ci","P","HR","CI_low","CI_high")]
tdf1$CHIP<-as.character(tdf1$CHIP)
tdf1$CHIP<-gsub("Baseline","No CHIP",tdf1$CHIP)
tdf1$CHIP<-gsub("CHMD","M-CHIP",tdf1$CHIP)
tdf1$CHIP<-gsub("M-CHIP-high","M-CHIP, VAF>=0.1",tdf1$CHIP)
tdf1$CHIP<-gsub("CHLD","L-CHIP",tdf1$CHIP)
tdf1$CHIP<-gsub("L-CHIP-high","L-CHIP, VAF>=0.1",tdf1$CHIP)
tdf1$CHIP<-gsub("MmCA","M-mCA",tdf1$CHIP)
tdf1$CHIP<-gsub("M-mCA-high","M-mCA, CF>=0.1",tdf1$CHIP)
tdf1$CHIP<-gsub("LmCA","L-mCA",tdf1$CHIP)
tdf1$CHIP<-gsub("L-mCA-high","L-mCA, CF>=0.1",tdf1$CHIP)
tdf1$CHIP<-gsub("AmCA","A-mCA",tdf1$CHIP)
tdf1$CHIP<-gsub("A-mCA-high","A-mCA, CF>=0.1",tdf1$CHIP)
tdf1$CHIP<-gsub("UmCA","U-mCA",tdf1$CHIP)
tdf1$CHIP<-gsub("U-mCA-high","U-mCA, CF>=0.1",tdf1$CHIP)

tdf1[tdf1$HR_ci == "1 (1, 1)",]$P<-""
tdf1[tdf1$HR_ci == "1 (1, 1)",]$HR_ci<-""


tdf2<-tdf1[,c("CHIP","at_risk","event","HR_ci","P")]
colnames(tdf2)<-c("","Total","Event","HR (95% CI)","P-value")

ptdf1<-ggtexttable(tdf2,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))

#print (tdf1)

tdf1$CHIP<-factor(as.character(tdf1$CHIP),levels=c("U-mCA, CF>=0.1","U-mCA","A-mCA, CF>=0.1","A-mCA","L-mCA, CF>=0.1","L-mCA","M-mCA, CF>=0.1","M-mCA","No mCA","L-CHIP, VAF>=0.1","L-CHIP","M-CHIP, VAF>=0.1","M-CHIP","No CHIP"))
F2b<- ggplot(data=tdf1, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey",size=0.5) +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep("grey",2),rep("orange",2),rep(ctList[1],2),rep(ctList[2],2),"black",rep(ctList[1],2),rep(ctList[2],2),"black"))


P0<-ggarrange(F2b,ptdf1,ncol=2,widths=c(1.8,2))


pdf(paste(figDir,"scripts/Final_figs/Fig_ED10c_CAD.pdf",sep=""),width=8,height=3)
print (P0)
dev.off()



}
