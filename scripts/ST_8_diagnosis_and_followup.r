
##===================================================================
##===================================================================


library(ggplot2)
library(data.table)
library(ggpubr)


setwd("")
hemeICDFile<-"meta/heme_malignancies_icd10.txt"
figDir<-""



## Incidence analysis
prepCox<-function(tukb){
	DT<-c()
	## Individuals who got diagnosis
	t0<-tukb[tukb$has_disease == 1,]
	DT0<-t0$Censor_date
	names(DT0)<-as.character(t0$eid)
	## Individuals who died
	t1<-tukb[tukb$has_disease == 0 & tukb$dead == 1, ]
	DT1<-t1$death_censor_date
	names(DT1)<-as.character(t1$eid)
	## Individuals who are alive but did not develop malignancy
	t2<-tukb[tukb$has_disease == 0 & tukb$dead == 0,]$eid
	DT2<-rep(forcedCensorDate,length(t2))
	names(DT2)<-as.character(t2)
	## Combine
	DT<-c(DT0,DT1,DT2)
	names(DT)<-c(names(DT0),names(DT1),names(DT2))
	## Censor date
	tukb$Censor_date<-DT[as.character(tukb$eid)]
	## Follow up start date
	tukb$followStart<-as.Date(tukb$date_attending_assessment_center)
	## Status
	tukb$Status<-tukb$has_disease
	tukb$Year<-as.numeric(difftime(tukb$Censor_date,tukb$followStart,units=c("days")))/365.25
	q95<-round(as.numeric(quantile(tukb$Year,0.95))-0.5)
	tukb[tukb$Year > q95,]$Status<-0
	tukb[tukb$Year > q95,]$Year<-q95
	## return
	return (tukb)
}



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


forcedCensorDate<-as.Date("2020-03-31")


##===================================================================
##ICD10codes
hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")
hicdlist<-as.character(hicd[hicd$Group %in% c("Myeloid","Lymphoid","Unspecified"),]$ICD)


#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
## FIG

ukbc<-data.frame(fread(paste(figDir,"Fig_1_WES_malignancy_df.txt",sep="")))


ukbc$Caucasian<-as.factor(ukbc$Caucasian)
ukbc$ever_smoked<-factor(ukbc$ever_smoked,levels=c("No","Yes"))
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))


## CH category
ukbc$CH<-rep("Control",nrow(ukbc))
ukbc[ukbc$CHIP == 1 & ukbc$MCHIP == 0,]$CH<-"Lymphoid"
ukbc[ukbc$MCHIP == 1 & ukbc$CHIP == 0,]$CH<-"Myeloid"
ukbc[ukbc$CHIP == 1 & ukbc$MCHIP == 1,]$CH <- "ML"


#=====================================================================
## Analyzing the effect of three groups of CH

# FUNC
getdfWES<-function(tukb,M,AGE){
	tdf1<-data.frame(Malignancy = M,
	Cohort = "WES",
	Age = AGE,
	Followup = c(median(tukb$Year),NA,NA,NA,NA),
	Group = c("All","Control","M-CHIP","L-CHIP","M-CHIP + L-CHIP"),
	N = c(NA,
		nrow(tukb[tukb$Status == 1 & tukb$CH == "Control",]),
		nrow(tukb[tukb$Status == 1 & tukb$CH == "Myeloid",]),
		nrow(tukb[tukb$Status == 1 & tukb$CH == "Lymphoid",]),
		nrow(tukb[tukb$Status == 1 & tukb$CH == "ML",])),
	Total = c(NA,
		nrow(tukb[tukb$CH == "Control",]),
		nrow(tukb[tukb$CH == "Myeloid",]),
		nrow(tukb[tukb$CH == "Lymphoid",]),
		nrow(tukb[tukb$CH == "ML",])),
	timeToDiagnosis = c(NA,
		median(tukb[tukb$CH == "Control" & tukb$Status == 1,]$Year),
		median(tukb[tukb$CH == "Myeloid" & tukb$Status == 1,]$Year),
		median(tukb[tukb$CH == "Lymphoid" & tukb$Status == 1,]$Year),
		median(tukb[tukb$CH == "ML" & tukb$Status == 1,]$Year)))
	return (tdf1)
}


tukb<-ukbc
tukb$has_disease<-tukb$Myeloid
tukb$Censor_date<-tukb$Myeloid_date

tukb<-prepCox(tukb)

tdf1<-getdfWES(tukb,"Myeloid","All")
tdf11<-getdfWES(tukb[tukb$Age < 50,],"Myeloid","<50")
tdf21<-getdfWES(tukb[tukb$Age >= 50 & tukb$Age <60,],"Myeloid","50-59")
tdf31<-getdfWES(tukb[tukb$Age >= 60,],"Myeloid",">=60")

#=====================================================================
#=====================================================================
#=========================================
## Lymphoid malignancies


tukb<-ukbc
tukb$has_disease<-tukb$Lymphoid
tukb$Censor_date<-tukb$Lymphoid_date

tukb<-prepCox(tukb)


tdf2<-getdfWES(tukb,"Lymphoid","All")
tdf12<-getdfWES(tukb[tukb$Age < 50,],"Lymphoid","<50")
tdf22<-getdfWES(tukb[tukb$Age >= 50 & tukb$Age < 60,],"Lymphoid","50-59")
tdf32<-getdfWES(tukb[tukb$Age >= 60,],"Lymphoid",">=60")


#=========================================
## CLL

tukb<-ukbc
tukb$has_disease<-tukb$CLL
tukb$Censor_date<-tukb$CLL_date

tukb<-prepCox(tukb)

tdf7<-getdfWES(tukb,"CLL/SLL","All")
tdf17<-getdfWES(tukb[tukb$Age < 50,],"CLL/SLL","<50")
tdf27<-getdfWES(tukb[tukb$Age >= 50 & tukb$Age < 60,],"CLL/SLL","50-59")
tdf37<-getdfWES(tukb[tukb$Age >= 60,],"CLL/SLL",">=60")



#=====================================================================
## Analyzing the effect of three groups of CH
## mCA

# FUNC
getdfSNP<-function(tukb,M,AGE){
	tdf1<-data.frame(Malignancy = M,
	Cohort = "SNP-array",
	Age = AGE,
	Followup = c(median(tukb$Year),NA,NA,NA,NA,NA),
	Group = c("All","Control","M-mCA","L-mCA","A-mCA","U-mCA"),
	N = c(NA,
		nrow(tukb[tukb$Status == 1 & tukb$CH == "Control",]),
		nrow(tukb[tukb$Status == 1 & tukb$CH == "Myeloid",]),
		nrow(tukb[tukb$Status == 1 & tukb$CH == "Lymphoid",]),
		nrow(tukb[tukb$Status == 1 & tukb$CH == "ML",]),
		nrow(tukb[tukb$Status == 1 & tukb$CH == "Unk",])),
	Total = c(NA,
		nrow(tukb[tukb$CH == "Control",]),
		nrow(tukb[tukb$CH == "Myeloid",]),
		nrow(tukb[tukb$CH == "Lymphoid",]),
		nrow(tukb[tukb$CH == "ML",]),
		nrow(tukb[tukb$CH == "Unk",])),
	timeToDiagnosis = c(NA,
		median(tukb[tukb$CH == "Control" & tukb$Status == 1,]$Year),
		median(tukb[tukb$CH == "Myeloid" & tukb$Status == 1,]$Year),
		median(tukb[tukb$CH == "Lymphoid" & tukb$Status == 1,]$Year),
		median(tukb[tukb$CH == "ML" & tukb$Status == 1,]$Year),
		median(tukb[tukb$CH == "Unk" & tukb$Status == 1,]$Year)))
	return (tdf1)
}

ukbc<-data.frame(fread(paste(figDir,"Fig_1_UKB500_malignancy_df.txt",sep="")))


ukbc$Caucasian<-as.factor(ukbc$Caucasian)
ukbc$ever_smoked<-factor(ukbc$ever_smoked,levels=c("No","Yes"))
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))


## CH category
ukbc$CH<-rep("Control",nrow(ukbc))
ukbc[ukbc$LCNV == 1 & ukbc$MCNV == 0,]$CH<-"Lymphoid"
ukbc[ukbc$MCNV == 1 & ukbc$LCNV == 0,]$CH<-"Myeloid"
ukbc[ukbc$LCNV == 1 & ukbc$MCNV == 1,]$CH <- "ML"
ukbc[ukbc$LCNV == 0 & ukbc$MCNV == 0 & ukbc$CNV == 1,]$CH <- "Unk"

#================================================================================
## Myeloid

tukb1<-ukbc
tukb1$has_disease<-tukb1$Myeloid
tukb1$Censor_date<-tukb1$Myeloid_date

tukb1<-prepCox(tukb1)


tdf3<-getdfSNP(tukb1,"Myeloid","All")
tdf13<-getdfSNP(tukb1[tukb1$Age < 50,],"Myeloid","<50")
tdf23<-getdfSNP(tukb1[tukb1$Age >= 50 & tukb1$Age < 60,],"Myeloid","50-59")
tdf33<-getdfSNP(tukb1[tukb1$Age >= 60,],"Myeloid",">=60")

#=====================================================================
#=====================================================================
## Lymphoid

tukb1<-ukbc
tukb1$has_disease<-tukb1$Lymphoid
tukb1$Censor_date<-tukb1$Lymphoid_date

tukb1<-prepCox(tukb1)


tdf4<-getdfSNP(tukb1,"Lymphoid","All")
tdf14<-getdfSNP(tukb1[tukb1$Age < 50,],"Lymphoid","<50")
tdf24<-getdfSNP(tukb1[tukb1$Age >= 50 & tukb1$Age < 60,],"Lymphoid","50-59")
tdf34<-getdfSNP(tukb1[tukb1$Age >= 60,],"Lymphoid",">=60")

#=====================================================================
## CLL/SLL

tukb1<-ukbc
tukb1$has_disease<-tukb1$CLL
tukb1$Censor_date<-tukb1$CLL_date

tukb1<-prepCox(tukb1)


tdf8<-getdfSNP(tukb1,"CLL/SLL","All")
tdf18<-getdfSNP(tukb1[tukb1$Age < 50,],"CLL/SLL","<50")
tdf28<-getdfSNP(tukb1[tukb1$Age >= 50 & tukb1$Age < 60,],"CLL/SLL","50-59")
tdf38<-getdfSNP(tukb1[tukb1$Age >= 60,],"CLL/SLL",">=60")


t<-rbind(tdf1,tdf11,tdf21,tdf31,tdf2,tdf12,tdf22,tdf32,tdf7,tdf17,tdf27,tdf37,tdf3,tdf13,tdf23,tdf33,tdf4,tdf14,tdf24,tdf34,tdf8,tdf18,tdf28,tdf38)

write.table(t,file = paste(figDir,"scripts/Final_figs/ST_X_followup_diagnosis.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

