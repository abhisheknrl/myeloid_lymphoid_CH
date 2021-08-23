

## Mortality analysis
## WES, 500, FIG, FIG_v2

args = commandArgs(trailingOnly=TRUE)

CMD<-args[1]


if (CMD == "WES"){

library(ggplot2)
library(data.table)


setwd("")
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"

## Updated death data
deathFile<-"meta/death.txt"
codFile<-"meta/death_cause.txt"


#=====================================================================
## Reload data

ukbc<-data.frame(fread(paste(figDir,"1_UKB_data.txt",sep="")))

## Adding covariates
ukbc$EverSmoke<-ukbc$ever_smoked
#Pack years
ukbc$Pack_year<-ukbc$pack_years_of_smoking
## isEuropean
ukbc$Caucasian<-rep(1,nrow(ukbc))
ukbc[is.na(ukbc$genetic_ethnic_grouping),]$Caucasian<-0

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
## Mortality analysis
## Exome publication data. Any individual who underwent MRI after this date would not be selected for exome sequencing based on available MRI
exomeDate<-"2019-03-31"

## Subtract the age difference between imaging visit and first visit

temp1<-ukbc[is.na(ukbc$date_attending_assessment_center_2),]
temp2<-ukbc[!is.na(ukbc$date_attending_assessment_center_2),]
temp2$Diff<-as.numeric(difftime(exomeDate, temp2$date_attending_assessment_center_2,units=c("days")))/365.25
temp11<-temp2[temp2$Diff <= 0,]
temp12<-temp2[temp2$Diff > 0,]
temp12$death_year<-temp12$death_year - as.numeric(difftime(temp12$date_attending_assessment_center_2,temp12$date_attending_assessment_center,units=c("days")))/365.25

ukbc<-rbind(temp1,temp11[,colnames(temp1)],temp12[,colnames(temp1)])

#===========================================================================
## 

tukb<-ukbc
tukb$Status<-tukb$dead
tukb$Year<-tukb$death_year

## Add the CHIP Data
tukb$CH<-rep("Control",nrow(tukb))
tukb[tukb$CHIP == 0 & tukb$MCHIP == 1,]$CH<-"CHMD"
tukb[tukb$CHIP == 1 & tukb$MCHIP == 0,]$CH<-"CHLD"
tukb[tukb$CHIP == 1 & tukb$MCHIP == 1,]$CH<-"ML"

tukb<-tukb[tukb$CH != "ML",]

## Quantile
q95<-as.numeric(quantile(tukb$Year,0.95))
tukb[tukb$Year >= q95,]$Status<-0
tukb[tukb$Year >= q95,]$Year<-q95


#=================================================================
## Annotate with malignancy data

hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")
myeloid_icdList<-as.character(hicd[hicd$Group == "Myeloid",]$ICD)
lymphoid_icdList<-as.character(hicd[hicd$Group == "Lymphoid",]$ICD)
unspe_icdList<-as.character(hicd[hicd$Group == "Unspecified",]$ICD)


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

## MY/LY
tdf1<-getDiag(tukb,c(myeloid_icdList,lymphoid_icdList,unspe_icdList),any10,date10)
tdf2<-getCancer(tukb,c(myeloid_icdList,lymphoid_icdList,unspe_icdList),cancer10,cancerDate10)
tdf<-getUnique(rbind(tdf1,tdf2))
## Myeloid/lmphoid
tukb$Heme<-rep(0,nrow(tukb))
tukb[tukb$eid %in% tdf$eid,]$Heme<-1
DT<-tdf$Diff/365.25
names(DT)<-as.character(tdf$eid)
tukb$Heme_date<-DT[as.character(tukb$eid)]

#=================================================================


### Write these data and quit R session and load the data again
CN<-c("eid","Age","Age_sq","Sex","EverSmoke","Caucasian",paste("PC",1:5,sep=""),"Status","Year","CH","CHIP","CHIP_size","MCHIP","Msize","lymphocyte_count","CNV","CNV_size","Heme","Heme_date")
write.table(tukb[,CN],file=paste(figDir,"scripts/data/Fig2_mortality_WES_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

}

#=================================================================
#=================================================================
#=================================================================
#=================================================================
## Mortality data for the 500K population

if (CMD == 500){

library(ggplot2)
library(data.table)


setwd("")
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"

## Updated death data
deathFile<-"meta/death.txt"
codFile<-"meta/death_cause.txt"

#=====================================================================
## Reload data

ukbc<-data.frame(fread(paste(figDir,"data/UKB_500K.txt",sep="")))

## Adding covariates
ukbc$EverSmoke<-ukbc$ever_smoked
#Pack years
ukbc$Pack_year<-ukbc$pack_years_of_smoking
## isEuropean
ukbc$Caucasian<-rep(1,nrow(ukbc))
ukbc[is.na(ukbc$genetic_ethnic_grouping),]$Caucasian<-0

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
## Mortality analysis

exomeDate<-"2019-03-31"

## Subtract the age difference between imaging visit and first visit

temp1<-ukbc[is.na(ukbc$date_attending_assessment_center_2),]
temp2<-ukbc[!is.na(ukbc$date_attending_assessment_center_2),]
temp2$Diff<-as.numeric(difftime(exomeDate, temp2$date_attending_assessment_center_2,units=c("days")))/365.25
temp11<-temp2[temp2$Diff <= 0,]
temp12<-temp2[temp2$Diff > 0,]
temp12$death_year<-temp12$death_year - as.numeric(difftime(temp12$date_attending_assessment_center_2,temp12$date_attending_assessment_center,units=c("days")))/365.25

ukbc<-rbind(temp1,temp11[,colnames(temp1)],temp12[,colnames(temp1)])

#===========================================================================
## 

tukb<-ukbc
tukb$Status<-tukb$dead
tukb$Year<-tukb$death_year

## Add the CHIP Data
tukb$CH<-rep("Control",nrow(tukb))
tukb[tukb$LCNV == 0 & tukb$MCNV == 1,]$CH<-"MmCA"
tukb[tukb$LCNV == 1 & tukb$MCNV == 0,]$CH<-"LmCA"
tukb[tukb$LCNV == 1 & tukb$MCNV == 1,]$CH<-"AmCA"
tukb[tukb$LCNV == 0 & tukb$MCNV == 0 & tukb$CNV == 1,]$CH<-"UmCA"


## Quantile
q95<-as.numeric(quantile(tukb$Year,0.95))
tukb[tukb$Year >= q95,]$Status<-0
tukb[tukb$Year >= q95,]$Year<-q95


#=================================================================
## Annotate with malignancy data

hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")
myeloid_icdList<-as.character(hicd[hicd$Group == "Myeloid",]$ICD)
lymphoid_icdList<-as.character(hicd[hicd$Group == "Lymphoid",]$ICD)
unspe_icdList<-as.character(hicd[hicd$Group == "Unspecified",]$ICD)


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

## MY/LY
tdf1<-getDiag(tukb,c(myeloid_icdList,lymphoid_icdList,unspe_icdList),any10,date10)
tdf2<-getCancer(tukb,c(myeloid_icdList,lymphoid_icdList,unspe_icdList),cancer10,cancerDate10)
tdf<-getUnique(rbind(tdf1,tdf2))
## Myeloid/lmphoid
tukb$Heme<-rep(0,nrow(tukb))
tukb[tukb$eid %in% tdf$eid,]$Heme<-1
DT<-tdf$Diff/365.25
names(DT)<-as.character(tdf$eid)
tukb$Heme_date<-DT[as.character(tukb$eid)]

#=================================================================


### Write these data and quit R session and load the data again
CN<-c("eid","Age","Age_sq","Sex","EverSmoke","Caucasian",paste("PC",1:5,sep=""),"Status","Year","CH","MCNV","MCNV_size","LCNV","LCNV_size","lymphocyte_count","CNV","CNV_size","Heme","Heme_date")

write.table(tukb[,CN],file=paste(figDir,"scripts/data/Fig2_mortality_500K_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

}


#=================================================================
#=================================================================
#=================================================================
#=================================================================
## Making the plots


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



#=================================================================
## WES data

tukb<-data.frame(fread(paste(figDir,"scripts/data/Fig2_mortality_WES_df.txt",sep="")))

## Factorize
tukb$CH<-factor(tukb$CH,levels=c("Control","CHMD","CHLD"))
tukb$Caucasian<-as.factor(tukb$Caucasian)

tukb<-tukb[!is.na(tukb$EverSmoke),]

# #####
YLAB<-"Mortality"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb)
p86<-ggsurvplot(fit,data=tukb,fun="event",
		legend.title = "CH", legend=c(0.2,0.75), legend.labs=c("Control","CHMD","CHLD"),
		palette=c("black",ctList[2],ctList[1]),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")



tukb$Age_categ<-rep("<50",nrow(tukb))
tukb[tukb$Age >= 50 & tukb$Age < 60,]$Age_categ<-"50-59"
tukb[tukb$Age >= 60,]$Age_categ<-">60"

## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb)
tukb2<-tukb[tukb$Msize != 1 & tukb$CHIP_size != 1,]
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb2)


## Tables
mvdf1<-extractCox(cox1,"Lymphoid","CHCHMD","CHMD","Any")
mvdf5<-extractCox(cox1,"Lymphoid","CHCHLD","CHLD","Any")
mvdf2<-extractCox(cox2,"Lymphoid","CHCHMD","CHMD-high","Any")
mvdf4<-extractCox(cox2,"Lymphoid","CHCHLD","CHLD-high","Any")
BL<-data.frame(Malignancy = "", CHIP = "Baseline", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5,mvdf2,mvdf4)
tdf$Total<-c(as.numeric(table(tukb$CH)[c("Control","CHMD","CHLD")]),as.numeric(table(tukb2$CH)[c("CHMD","CHLD")]))
tdf$Positive<-c(as.numeric(table(tukb[tukb$Status ==1,]$CH)[c("Control","CHMD","CHLD")]),as.numeric(table(tukb2[tukb2$Status ==1,]$CH)[c("CHMD","CHLD")]))
rownames(tdf)<-as.character(tdf$CHIP)
tdf<-tdf[c("Baseline","CHMD","CHMD-high","CHLD","CHLD-high"),]

tdf$CHIP<-factor(tdf$CHIP,levels=c("CHLD-high","CHLD","CHMD-high","CHMD","Baseline"))

p87 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[1],2),rep(ctList[2],2),"black"))

##

tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1a<-p86
F1b<-ggarrange(p87,
	ggarrange(ptdf,p87,ncol=2,nrow=1,widths=c(2,1)),
	nrow=2,heights=c(1,2))



F2<-ggarrange(F1a$plot,F1b,ncol=2,nrow=1,widths=c(1,2))
mydf<-tdf

#=================================================================
## 500K CNV data

tukb1<-data.frame(fread(paste(figDir,"scripts/data/Fig2_mortality_500K_df.txt",sep="")))

## Factorize
tukb1$CH<-factor(tukb1$CH,levels=c("Control","MmCA","LmCA","AmCA","UmCA"))
tukb1$Caucasian<-as.factor(tukb1$Caucasian)


tukb1<-tukb1[!is.na(tukb1$EverSmoke),]

# #####
YLAB<-"Mortality"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb1)
p96<-ggsurvplot(fit,data=tukb1,fun="event",
		legend.title = "CH", legend=c(0.2,0.75), legend.labs=c("Control","MmCA","LmCA","AmCA","UmCA"),
		palette=c("black",ctList[2],ctList[1],"orange","grey"),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")




tukb1$Age_categ<-rep("<50",nrow(tukb1))
tukb1[tukb1$Age >= 50 & tukb1$Age < 60,]$Age_categ<-"50-59"
tukb1[tukb1$Age >= 60,]$Age_categ<-">60"
tukb1$Age_categ<-factor(tukb1$Age_categ,levels=c("<50","50-59",">60"))

## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1)
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex +EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1[tukb1$CNV_size != 1,])

## Tables
mvdf1<-extractCox(cox1,"Lymphoid",c("CHMmCA","CHLmCA"),c("MmCA","LmCA"),"Any")
print ("ok")
mvdf5<-extractCox(cox1,"Lymphoid",c("CHAmCA","CHUmCA"),c("AmCA","UmCA"),"Any")
mvdf2<-extractCox(cox2,"Lymphoid",c("CHMmCA","CHLmCA"),c("MmCA-high","LmCA-high"),"Any")
mvdf4<-extractCox(cox2,"Lymphoid",c("CHAmCA","CHUmCA"),c("AmCA-high","UmCA-high"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "No mCA", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5,mvdf2,mvdf4)
tdf$Total<-c(as.numeric(table(tukb1$CH)[c("Control","MmCA","LmCA","AmCA","UmCA")]),as.numeric(table(tukb1[tukb1$CNV_size != 1,]$CH)[c("MmCA","LmCA","AmCA","UmCA")]))
tdf$Positive<-c(as.numeric(table(tukb1[tukb1$Status ==1,]$CH)[c("Control","MmCA","LmCA","AmCA","UmCA")]),as.numeric(table(tukb1[tukb1$CNV_size != 1 & tukb1$Status ==1,]$CH)[c("MmCA","LmCA","AmCA","UmCA")]))
rownames(tdf)<-as.character(tdf$CHIP)

tdf<-tdf[c("No mCA","MmCA","MmCA-high","LmCA","LmCA-high","AmCA","AmCA-high","UmCA","UmCA-high"),]

tdf$CHIP<-factor(tdf$CHIP,levels=c("UmCA-high","UmCA","AmCA-high","AmCA","LmCA-high","LmCA","MmCA-high","MmCA","No mCA"))

p97 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep("grey",2),rep("orange",2),rep(ctList[1],2),rep(ctList[2],2),"black"))

##

tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf1<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1c<-p96
F1d<-ggarrange(ptdf1,p97,ncol=2,nrow=1,widths=c(2,1))

F3<-ggarrange(F1c$plot,F1d,ncol=2,nrow=1,widths=c(1,2))

#=================================================================
#=================================================================
## Combining the data

tdf<-rbind(mydf,tdf)

tdf1<-tdf[,c("CHIP","Total","Positive","HR_ci","P","HR","CI_low","CI_high")]
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


tdf2<-tdf1[,c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf2)<-c("","Total","Event","HR (95% CI)","P-value")


ptdf1<-ggtexttable(tdf2,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))

tdf1$CHIP<-factor(as.character(tdf1$CHIP),levels=c("U-mCA, CF>=0.1","U-mCA","A-mCA, CF>=0.1","A-mCA","L-mCA, CF>=0.1","L-mCA","M-mCA, CF>=0.1","M-mCA","No mCA","L-CHIP, VAF>=0.1","L-CHIP","M-CHIP, VAF>=0.1","M-CHIP","No CHIP"))
F2a<- ggplot(data=tdf1, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep("grey",2),rep("orange",2),rep(ctList[1],2),rep(ctList[2],2),"black",rep(ctList[1],2),rep(ctList[2],2),"black"))

P<-ggarrange(F2a,ptdf1,ncol=2,widths=c(1.8,2))



pdf(paste(figDir,"scripts/Final_figs/Fig_ED10a_mortality.pdf",sep=""),width=8,height=3)
print (P)
dev.off()




}


#=================================================================
## Making the plots


if (CMD == "FIG_v2"){

## Unrelated to heme malignancy

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

## Updated death data
deathFile<-"meta/death.txt"
codFile<-"meta/death_cause.txt"


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


cod<-read.delim(codFile,head=TRUE,sep="\t")

hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")
myeloid_icdList<-as.character(hicd[hicd$Group == "Myeloid",]$ICD)
lymphoid_icdList<-as.character(hicd[hicd$Group == "Lymphoid",]$ICD)
unspe_icdList<-as.character(hicd[hicd$Group == "Unspecified",]$ICD)

IL<-c(myeloid_icdList,lymphoid_icdList,unspe_icdList)
cod<-cod[cod$cause_icd10 %in% IL,]
hemeEID<-unique(cod$eid)

#=================================================================
## WES data

tukb<-data.frame(fread(paste(figDir,"scripts/data/Fig2_mortality_WES_df.txt",sep="")))

## Factorize
tukb$CH<-factor(tukb$CH,levels=c("Control","CHMD","CHLD"))
tukb$Caucasian<-as.factor(tukb$Caucasian)


tukb<-tukb[!is.na(tukb$EverSmoke) & tukb$Heme != 1 & (!tukb$eid %in% hemeEID),]

# #####
YLAB<-"Mortality"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb)
p86<-ggsurvplot(fit,data=tukb,fun="event",
		legend.title = "CH", legend=c(0.2,0.75), legend.labs=c("Control","CHMD","CHLD"),
		palette=c("black",ctList[2],ctList[1]),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")



tukb$Age_categ<-rep("<50",nrow(tukb))
tukb[tukb$Age >= 50 & tukb$Age < 60,]$Age_categ<-"50-59"
tukb[tukb$Age >= 60,]$Age_categ<-">60"

## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb)
tukb2<-tukb[tukb$Msize != 1 & tukb$CHIP_size != 1,]
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb2)


## Tables
mvdf1<-extractCox(cox1,"Lymphoid","CHCHMD","CHMD","Any")
mvdf5<-extractCox(cox1,"Lymphoid","CHCHLD","CHLD","Any")
mvdf2<-extractCox(cox2,"Lymphoid","CHCHMD","CHMD-high","Any")
mvdf4<-extractCox(cox2,"Lymphoid","CHCHLD","CHLD-high","Any")
BL<-data.frame(Malignancy = "", CHIP = "Baseline", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5,mvdf2,mvdf4)
tdf$Total<-c(as.numeric(table(tukb$CH)[c("Control","CHMD","CHLD")]),as.numeric(table(tukb2$CH)[c("CHMD","CHLD")]))
tdf$Positive<-c(as.numeric(table(tukb[tukb$Status ==1,]$CH)[c("Control","CHMD","CHLD")]),as.numeric(table(tukb2[tukb2$Status ==1,]$CH)[c("CHMD","CHLD")]))
rownames(tdf)<-as.character(tdf$CHIP)
tdf<-tdf[c("Baseline","CHMD","CHMD-high","CHLD","CHLD-high"),]

tdf$CHIP<-factor(tdf$CHIP,levels=c("CHLD-high","CHLD","CHMD-high","CHMD","Baseline"))

p87 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[1],2),rep(ctList[2],2),"black"))

##

tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1a<-p86
F1b<-ggarrange(p87,
	ggarrange(ptdf,p87,ncol=2,nrow=1,widths=c(2,1)),
	nrow=2,heights=c(1,2))



F2<-ggarrange(F1a$plot,F1b,ncol=2,nrow=1,widths=c(1,2))
mydf<-tdf

#=================================================================
## 500K CNV data

tukb1<-data.frame(fread(paste(figDir,"scripts/data/Fig2_mortality_500K_df.txt",sep="")))

## Factorize
tukb1$CH<-factor(tukb1$CH,levels=c("Control","MmCA","LmCA","AmCA","UmCA"))
tukb1$Caucasian<-as.factor(tukb1$Caucasian)


tukb1<-tukb1[!is.na(tukb1$EverSmoke) & tukb1$Heme != 1 & (!tukb1$eid %in% hemeEID),]

# #####
YLAB<-"Mortality"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb1)
p96<-ggsurvplot(fit,data=tukb1,fun="event",
		legend.title = "CH", legend=c(0.2,0.75), legend.labs=c("Control","MmCA","LmCA","AmCA","UmCA"),
		palette=c("black",ctList[2],ctList[1],"orange","grey"),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")




tukb1$Age_categ<-rep("<50",nrow(tukb1))
tukb1[tukb1$Age >= 50 & tukb1$Age < 60,]$Age_categ<-"50-59"
tukb1[tukb1$Age >= 60,]$Age_categ<-">60"
tukb1$Age_categ<-factor(tukb1$Age_categ,levels=c("<50","50-59",">60"))

## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1)
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1[tukb1$CNV_size != 1,])



## Tables
mvdf1<-extractCox(cox1,"Lymphoid",c("CHMmCA","CHLmCA"),c("MmCA","LmCA"),"Any")
print ("ok")
mvdf5<-extractCox(cox1,"Lymphoid",c("CHAmCA","CHUmCA"),c("AmCA","UmCA"),"Any")
mvdf2<-extractCox(cox2,"Lymphoid",c("CHMmCA","CHLmCA"),c("MmCA-high","LmCA-high"),"Any")
mvdf4<-extractCox(cox2,"Lymphoid",c("CHAmCA","CHUmCA"),c("AmCA-high","UmCA-high"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "No mCA", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5,mvdf2,mvdf4)
tdf$Total<-c(as.numeric(table(tukb1$CH)[c("Control","MmCA","LmCA","AmCA","UmCA")]),as.numeric(table(tukb1[tukb1$CNV_size != 1,]$CH)[c("MmCA","LmCA","AmCA","UmCA")]))
tdf$Positive<-c(as.numeric(table(tukb1[tukb1$Status ==1,]$CH)[c("Control","MmCA","LmCA","AmCA","UmCA")]),as.numeric(table(tukb1[tukb1$CNV_size != 1 & tukb1$Status ==1,]$CH)[c("MmCA","LmCA","AmCA","UmCA")]))
rownames(tdf)<-as.character(tdf$CHIP)

tdf<-tdf[c("No mCA","MmCA","MmCA-high","LmCA","LmCA-high","AmCA","AmCA-high","UmCA","UmCA-high"),]

tdf$CHIP<-factor(tdf$CHIP,levels=c("UmCA-high","UmCA","AmCA-high","AmCA","LmCA-high","LmCA","MmCA-high","MmCA","No mCA"))

p97 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep("grey",2),rep("orange",2),rep(ctList[1],2),rep(ctList[2],2),"black"))

##

tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf1<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1c<-p96
F1d<-ggarrange(ptdf1,p97,ncol=2,nrow=1,widths=c(2,1))

F3<-ggarrange(F1c$plot,F1d,ncol=2,nrow=1,widths=c(1,2))

#=================================================================
#=================================================================
## Combining the data

tdf<-rbind(mydf,tdf)

tdf1<-tdf[,c("CHIP","Total","Positive","HR_ci","P","HR","CI_low","CI_high")]
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



tdf2<-tdf1[,c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf2)<-c("","Total","Event","HR (95% CI)","P-value")

ptdf1<-ggtexttable(tdf2,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))

tdf1$CHIP<-factor(as.character(tdf1$CHIP),levels=c("U-mCA, CF>=0.1","U-mCA","A-mCA, CF>=0.1","A-mCA","L-mCA, CF>=0.1","L-mCA","M-mCA, CF>=0.1","M-mCA","No mCA","L-CHIP, VAF>=0.1","L-CHIP","M-CHIP, VAF>=0.1","M-CHIP","No CHIP"))
F2a<- ggplot(data=tdf1, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey",size=0.5) +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep("grey",2),rep("orange",2),rep(ctList[1],2),rep(ctList[2],2),"black",rep(ctList[1],2),rep(ctList[2],2),"black"))

P<-ggarrange(F2a,ptdf1,ncol=2,widths=c(1.8,2))


pdf(paste(figDir,"scripts/Final_figs/Fig_ED10b_noHeme.pdf",sep=""),width=8,height=3)
print (P)
dev.off()


}
