
## CNV, FIG

args = commandArgs(trailingOnly=TRUE)

CMD<-args[1]


if (CMD == "CNV"){

library(ggplot2)
library(data.table)
library(colortools)
library(survival)
library(survminer)
library(ggpubr)

ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
hemeICDFile<-"meta/heme_malignancies_icd10.txt"

figDir<-""


##===================================================================
##ICD10codes
hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")

cnv<-data.frame(fread(paste(figDir,"data/UKB_500K_CNV.txt",sep="")))
cnv<-cnv[cnv$CHR %in% 1:22,]


LYG<-c("ATM","CIITA","ITPKB/NTRK1","KMT2C","MIR16-1","NCOR2","NOTCH1")
MYG<-c("CBL","CTCF","EP300","JAK2","TP53","MPL/GNB1")

## canonical myeloid/lymphoid
canLy<-c("gain12q","gain15q","gain17q","gain21q","gain22q13.32","gain2p","gain3q","gain8q","gain9q","del10p","del10q","del11q","del13q","del14q","del15q","del17p","del1p","del1q","del22q","del6q","del7q","del8p","tri12","tri18","tri19")
canMy<-c("gain1q","gain21q","gain9p","del12q","del20q","del5q","tri8")

##
manMG<-c("TP53","CTCF","MPL/GNB1","CBL")
manLG<-c("TCL1A","ATM","HEATR3/PLCG2/IRF8")
cnv$LCNV<-rep(0,nrow(cnv))
cnv[cnv$lyGenes %in% c(LYG,manMG) | cnv$myGenes %in% c(LYG,manMG) | cnv$Canonical_CA %in% canLy,]$LCNV<-1
cnv$MCNV<-rep(0,nrow(cnv))
cnv[cnv$myGenes %in% c(MYG,manLG) | cnv$lyGenes %in% c(MYG,manLG) | cnv$Canonical_CA %in% canMy,]$MCNV<-1


## UKBC
ukbc<-data.frame(fread(paste(figDir,"data/UKB_500K.txt",sep="")))

print (dim(ukbc))

cnv<-cnv[cnv$eid %in% ukbc$eid,]

write.table(cnv,file=paste(figDir,"scripts/Final_figs/UKB_CNV_submission.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

}

##===================================================================
##===================================================================

if (CMD == "WES"){


library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)
library(mgcv)


library(survival)
library(survminer)


ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
hemeICDFile<-"meta/heme_malignancies_icd10.txt"
figDir<-""


##===================================================================
##ICD10codes
hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")
hicdlist<-as.character(hicd[hicd$Group %in% c("Myeloid","Lymphoid","Unspecified"),]$ICD)


#=====================================================================
## Reload data

ukbc<-data.frame(fread(paste(figDir,"1_UKB_data.txt",sep="")))


#=====================================================================
#=====================================================================
## Functions

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


forcedCensorDate<-as.Date("2020-03-31")


#=====================================================================
## Filtering out related samples and non-caucasian individuals

ukbc$Caucasian<-rep(1,nrow(ukbc))
ukbc[is.na(ukbc$genetic_ethnic_grouping),]$Caucasian<-0
ukbc$Caucasian<-as.factor(ukbc$Caucasian)


## Adding covariates
ukbc$EverSmoke<-ukbc$ever_smoked
ukbc<-ukbc[!is.na(ukbc$EverSmoke),]


## PCA
CN<-colnames(ukbc)
CN<-gsub("genetic_pca","PC",CN)
colnames(ukbc)<-CN

## Age
ukbc$Age<-ukbc$age_attended_assessment_center
ukbc$Age_sq<-ukbc$Age*ukbc$Age


## Factorize the variables
ukbc$CHIP<-as.factor(ukbc$CHIP)
ukbc$MCHIP<-as.factor(ukbc$MCHIP)
ukbc$CHIP_size<-as.factor(ukbc$CHIP_size)
ukbc$Msize<-as.factor(ukbc$Msize)


## Age categorical
ukbc$Age_categ<-rep("<50",nrow(ukbc))
ukbc[ukbc$Age >= 50 & ukbc$Age < 60,]$Age_categ<-"50-59"
ukbc[ukbc$Age >= 60,]$Age_categ<-">60"
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))

#=====================================================================
## Get diagnosis

myeloid_icdList<-as.character(hicd[hicd$Group == "Myeloid",]$ICD)
lymphoid_icdList<-as.character(hicd[hicd$Group == "Lymphoid",]$ICD)

## Get diagnosis
tdf1<-getDiag(ukbc,c(myeloid_icdList,lymphoid_icdList),any10,date10)
tdf2<-getCancer(ukbc,c(myeloid_icdList,lymphoid_icdList),cancer10,cancerDate10)
tempdf<-rbind(tdf1,tdf2)
tempdf$Diff2<-as.numeric(difftime(forcedCensorDate,tempdf$Date,units=c("days")))
tempdf<-tempdf[tempdf$Diff2 >= 0,]
tempdf<-tempdf[,colnames(tdf1)]

## Myeloid malignancies
tdf<-getUnique(tempdf)
tdf<-tdf[tdf$ICD %in% myeloid_icdList,]
ukbc$Myeloid<-rep(0,nrow(ukbc))
ukbc[ukbc$eid %in% tdf$eid,]$Myeloid<-1
DT<-tdf$Date
names(DT)<-as.character(tdf$eid)
ukbc$Myeloid_date<-DT[as.character(ukbc$eid)]


## Lymphoid diagnosis
tdf<-getUnique(tempdf)
tdf<-tdf[tdf$ICD %in% lymphoid_icdList,]
ukbc$Lymphoid<-rep(0,nrow(ukbc))
ukbc[ukbc$eid %in% tdf$eid,]$Lymphoid<-1
DT<-tdf$Date
names(DT)<-as.character(tdf$eid)
ukbc$Lymphoid_date<-DT[as.character(ukbc$eid)]


## Lymphoid diagnosis
CLL<-c("C911","C830")
tdf<-getUnique(tempdf)
tdf<-tdf[tdf$ICD %in% CLL,]
ukbc$CLL<-rep(0,nrow(ukbc))
ukbc[ukbc$eid %in% tdf$eid,]$CLL<-1
DT<-tdf$Date
names(DT)<-as.character(tdf$eid)
ukbc$CLL_date<-DT[as.character(ukbc$eid)]


write.table(ukbc,file=paste(figDir,"Fig_1_WES_malignancy_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

}


#=====================================================================
#=====================================================================
#=====================================================================

if (CMD == 500){



library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)
library(mgcv)


library(survival)
library(survminer)


ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
hemeICDFile<-"meta/heme_malignancies_icd10.txt"
figDir<-""


##===================================================================
##ICD10codes
hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")
hicdlist<-as.character(hicd[hicd$Group %in% c("Myeloid","Lymphoid","Unspecified"),]$ICD)


#=====================================================================
## Functions

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


forcedCensorDate<-as.Date("2020-03-31")



#=====================================================================
### Heme malignancies enrichment
## CNVs

ukbc<-data.frame(fread(paste(figDir,"data/UKB_500K.txt",sep="")))


ukbc$Caucasian<-rep(1,nrow(ukbc))
ukbc[is.na(ukbc$genetic_ethnic_grouping),]$Caucasian<-0
ukbc$Caucasian<-as.factor(ukbc$Caucasian)


## Adding covariates
ukbc$EverSmoke<-ukbc$ever_smoked
ukbc<-ukbc[!is.na(ukbc$EverSmoke),]


## PCA
CN<-colnames(ukbc)
CN<-gsub("genetic_pca","PC",CN)
colnames(ukbc)<-CN

## Age
ukbc$Age<-ukbc$age_attended_assessment_center
ukbc$Age_sq<-ukbc$Age*ukbc$Age


## Age categorical
ukbc$Age_categ<-rep("<50",nrow(ukbc))
ukbc[ukbc$Age >= 50 & ukbc$Age < 60,]$Age_categ<-"50-59"
ukbc[ukbc$Age >= 60,]$Age_categ<-">60"
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))

#=====================================================================
## Get diagnosis

myeloid_icdList<-as.character(hicd[hicd$Group == "Myeloid",]$ICD)
lymphoid_icdList<-as.character(hicd[hicd$Group == "Lymphoid",]$ICD)

## Get diagnosis
tdf1<-getDiag(ukbc,c(myeloid_icdList,lymphoid_icdList),any10,date10)
tdf2<-getCancer(ukbc,c(myeloid_icdList,lymphoid_icdList),cancer10,cancerDate10)
tempdf<-rbind(tdf1,tdf2)
tempdf$Diff2<-as.numeric(difftime(forcedCensorDate,tempdf$Date,units=c("days")))
tempdf<-tempdf[tempdf$Diff2 >= 0,]
tempdf<-tempdf[,colnames(tdf1)]

## Myeloid malignancies
tdf<-getUnique(tempdf)
tdf<-tdf[tdf$ICD %in% myeloid_icdList,]
ukbc$Myeloid<-rep(0,nrow(ukbc))
ukbc[ukbc$eid %in% tdf$eid,]$Myeloid<-1
DT<-tdf$Date
names(DT)<-as.character(tdf$eid)
ukbc$Myeloid_date<-DT[as.character(ukbc$eid)]


## Lymphoid diagnosis
tdf<-getUnique(tempdf)
tdf<-tdf[tdf$ICD %in% lymphoid_icdList,]
ukbc$Lymphoid<-rep(0,nrow(ukbc))
ukbc[ukbc$eid %in% tdf$eid,]$Lymphoid<-1
DT<-tdf$Date
names(DT)<-as.character(tdf$eid)
ukbc$Lymphoid_date<-DT[as.character(ukbc$eid)]

## Lymphoid diagnosis
CLL<-c("C911","C830")
tdf<-getUnique(tempdf)
tdf<-tdf[tdf$ICD %in% CLL,]
ukbc$CLL<-rep(0,nrow(ukbc))
ukbc[ukbc$eid %in% tdf$eid,]$CLL<-1
DT<-tdf$Date
names(DT)<-as.character(tdf$eid)
ukbc$CLL_date<-DT[as.character(ukbc$eid)]


write.table(ukbc,file=paste(figDir,"Fig_1_UKB500_malignancy_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

}

##===================================================================
##===================================================================


if (CMD == "FIG"){


library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)
library(mgcv)


library(survival)
library(survminer)


ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


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
## Reload data

ukbc<-data.frame(fread(paste(figDir,"1_UKB_data.txt",sep="")))

chip<-read.table(paste(figDir,"1_CHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
mychip<-read.table(paste(figDir,"1_MCHIP_UKB.txt",sep=""),head=TRUE,sep="\t")


##===================================================================
## Whole exome

#mod1<-gam(CHIP ~ s(age_attended_assessment_center,bs="cr"),data=ukb)
#mod2<-gam(MCHIP ~ s(age_attended_assessment_center,bs="cr"),data=ukb)
t<-ukbc[,c("age_attended_assessment_center","CHIP","MCHIP")]
colnames(t)<-c("Age","L_CHIP","M_CHIP")

p0<-ggplot(t,aes(Age,L_CHIP))+
	geom_smooth(method="gam",formula = y ~ s(x,bs="cr"),color=ctList[1],fill=ctList[1],alpha=0.2)+
	geom_smooth(data=t,aes(Age,M_CHIP),method="gam",formula = y ~ s(x,bs="cr"),color=ctList[2],fill=ctList[2],alpha=0.2)+
	annotate("text", x = 65, y=0.12, 
		label = paste("M-CHIP",sep=""),
		size=5,color=ctList[2])+
	annotate("text", x = 65, y=0.03, 
		label = paste("L-CHIP",sep=""),
		size=5,color=ctList[1])+
	theme_linedraw()+
	theme(panel.grid=element_blank())+
	xlab("Age")+
	ylab("Proportion with CHIP")


write.table(t,file = paste(figDir,"scripts/Final_figs/figure_df/F1_a.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================

getFig<-function(temp1,mycol,XLAB,isMY="F"){
	#==============================
	## genes
	T<-table(as.character(temp1$Hugo_Symbol))
	T<-T[order(T,decreasing=TRUE)]
	K<-1
	for (gene in names(T)[1:25]){
		temp<-temp1[temp1$Hugo_Symbol == gene,]
		##
		tempdf<-data.frame(Gene = gene,
			Count = as.numeric(T[gene]),
			N = nrow(temp))
		if (K==1){
			K<-2
			tdf<-tempdf
		}else{
			tdf<-rbind(tdf,tempdf)
		}
	}
	##
	tdf$Gene<-factor(tdf$Gene,levels=names(T[order(T,decreasing=TRUE)]))
	sf2d<-ggplot(tdf,aes(y=Count,x=Gene))+
		geom_bar(stat="identity",color=NA,fill=mycol)+
		theme_classic()+
		xlab(XLAB)+
		ylab("Number of mutations")+
		scale_y_log10()+
		theme(legend.position="none",
			panel.grid.major.y=element_line(color="gray80"))+
		theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=rel(0.8)),
			axis.text.y=element_text(color="black"),
			axis.ticks = element_line(color="black"))
	if (isMY == "T"){
		sf2d<-ggplot(tdf,aes(y=Count,x=Gene))+
		geom_hline(yintercept = 30,color="grey")+
		geom_bar(stat="identity",color=NA,fill=mycol)+
		theme_classic()+
		xlab(XLAB)+
		ylab("Number of mutations")+
		scale_y_log10()+
		theme(legend.position="none",
			panel.grid.major.y=element_line(color="gray80"))+
		theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=rel(0.8)),
			axis.text.y=element_text(color="black"),
			axis.ticks = element_line(color="black"))
	}
	return (list(sf2d,tdf))
}


F1b<-getFig(chip,ctList[1],"L-CHIP genes")[[1]]
F1a<-getFig(mychip,ctList[2],"M-CHIP genes","T")[[1]]

t1<-getFig(chip,ctList[1],"L-CHIP genes")[[2]]
t2<-getFig(mychip,ctList[2],"M-CHIP genes","T")[[2]]
t1$CHIP<-rep("L-CHIP",nrow(t1))
t2$CHIP<-rep("M-CHIP",nrow(t2))

t<-rbind(t1,t2)
t<-t[,c("CHIP","Gene","Count")]

write.table(t,file = paste(figDir,"scripts/Final_figs/figure_df/F1_b.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


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

ukbc<-ukbc[ukbc$CH != "ML",]
ukbc$CH<-factor(ukbc$CH,levels=c("Control","Myeloid","Lymphoid"))


#=====================================================================
## Analyzing the effect of three groups of CH

tukb<-ukbc
tukb$has_disease<-tukb$Myeloid
tukb$Censor_date<-tukb$Myeloid_date

tukb<-prepCox(tukb)

tmyeloid<-tukb[,c("eid","CH","Status","Year")]
colnames(tmyeloid)<-c("eid","CHIP","Status_MY","Year_MY")

# #####
YLAB<-"Cumulative incidence of\nmyeloid malignancies"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb)
p98<-ggsurvplot(fit,data=tukb,fun="event",
		legend.title = "CH", legend=c(0.2,0.75), legend.labs=c("Baseline","M-CHIP","L-CHIP"),
		palette=c("black",ctList[2],ctList[1]),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")


## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb[tukb$CH != "Lymphoid",])
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb[tukb$CH != "Myeloid",])



# ## Tables
mvdf1<-extractCox(cox1,"Myeloid","CHMyeloid","M-CHIP","Any")
mvdf5<-extractCox(cox2,"Myeloid","CHLymphoid","L-CHIP","Any")
#mvdf5<-data.frame(Malignancy = "", CHIP = "Lymphoid", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
BL<-data.frame(Malignancy = "", CHIP = "No CHIP", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5)
tdf$Total<-as.numeric(table(tukb$CH)[c("Control","Myeloid","Lymphoid")])
tdf$Positive<-as.numeric(table(tukb[tukb$Status ==1,]$CH)[c("Control","Myeloid","Lymphoid")])

tdf$CHIP<-factor(tdf$CHIP,levels=c("L-CHIP","M-CHIP","No CHIP"))

if (nrow(tdf[tdf$Positive == 1,]) >= 1){
	tdf[tdf$Positive == 1,]$HR<-NA
	tdf[tdf$Positive == 1,]$CI_high<-NA
	tdf[tdf$Positive == 1,]$CI_low<-NA
	tdf[tdf$Positive == 1,]$P<-NA
}

# #tdf[tdf == Inf]<-0

p99 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		#panel.grid.major.x=element_line(color="gray80",size=0.2),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(ctList[1],ctList[2],"black"))


forestdf1<-tdf[,c("CHIP","HR","CI_low","CI_high")]
forestdf1$Malignancy<-rep("Myeloid")


tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1c<-p98
F1d<-ggarrange(ggarrange(p99,ptdf,ncol=2,nrow=1,widths=c(1,1.5)),
	nrow=2,heights=c(2,1))


#=====================================================================
#=====================================================================
#=========================================
## Lymphoid malignancies


tukb<-ukbc
tukb$has_disease<-tukb$Lymphoid
tukb$Censor_date<-tukb$Lymphoid_date

tukb<-prepCox(tukb)


ST<-tukb$Status
YR<-tukb$Year
names(ST)<-names(YR)<-as.character(tukb$eid)

tmyeloid$Status_LY<-ST[as.character(tmyeloid$eid)]
tmyeloid$Year_LY<-YR[as.character(tmyeloid$eid)]
write.table(tmyeloid,file = paste(figDir,"scripts/Final_figs/figure_df/F1_cd_incidence.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

# #####
YLAB<-"Cumulative incidence of\nlymphoid malignancies"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb)
p88<-ggsurvplot(fit,data=tukb,fun="event",
		legend.title = "CH", legend=c(0.2,0.75), legend.labs=c("Control","M-CHIP","L-CHIP"),
		palette=c("black",ctList[2],ctList[1]),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")


## CHIP size
tukb$Caucasian<-as.factor(tukb$Caucasian)
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb[tukb$CH != "Lymphoid",])
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb[tukb$CH != "Myeloid",])


## Tables
mvdf1<-extractCox(cox1,"Lymphoid","CHMyeloid","M-CHIP","Any")
mvdf5<-extractCox(cox2,"Lymphoid","CHLymphoid","L-CHIP","Any")
BL<-data.frame(Malignancy = "", CHIP = "No CHIP", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5)
tdf$Total<-as.numeric(table(tukb$CH)[c("Control","Myeloid","Lymphoid")])
tdf$Positive<-as.numeric(table(tukb[tukb$Status ==1,]$CH)[c("Control","Myeloid","Lymphoid")])

tdf$CHIP<-factor(tdf$CHIP,levels=c("L-CHIP","M-CHIP","No CHIP"))


p89 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		#panel.grid.major.x=element_line(color="gray80",size=0.2),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(ctList[1],ctList[2],"black"))


forestdf2<-tdf[,c("CHIP","HR","CI_low","CI_high")]
forestdf2$Malignancy<-rep("Lymphoid")
##

tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1e<-p88
F1f<-ggarrange(ggarrange(p89,ptdf,ncol=2,nrow=1,widths=c(1,1.5)),
	nrow=2,heights=c(2,1))



t<-rbind(forestdf1,forestdf2)
write.table(t,file = paste(figDir,"scripts/Final_figs/figure_df/F1_cd_forestPlots.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


#=====================================================================
## Analyzing the effect of three groups of CH
## mCA

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

ukbc<-ukbc[!ukbc$CH %in% c("ML","Unk"),]
ukbc$CH<-factor(ukbc$CH,levels=c("Control","Myeloid","Lymphoid"))



#================================================================================
## Myeloid

tukb1<-ukbc
tukb1$has_disease<-tukb1$Myeloid
tukb1$Censor_date<-tukb1$Myeloid_date

tukb1<-prepCox(tukb1)


tmyeloid<-tukb1[,c("eid","CH","Status","Year")]
colnames(tmyeloid)<-c("eid","CHIP","Status_MY","Year_MY")

# #####
YLAB<-"Cumulative incidence of\nmyeloid malignancies"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb1)
p96<-ggsurvplot(fit,data=tukb1,fun="event",
		legend.title = "mCA", legend=c(0.2,0.75), legend.labs=c("Control","M-mCA","L-mCA"),
		palette=c("black",ctList[2],ctList[1]),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")


## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + EverSmoke + Sex + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1[tukb1$CH != "Lymphoid",])
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + EverSmoke + Sex + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1[tukb1$CH != "Myeloid",])


# ## Tables
mvdf1<-extractCox(cox1,"Myeloid","CHMyeloid","M-mCA","Any")
mvdf5<-extractCox(cox2,"Myeloid","CHLymphoid","L-mCA","Any")
#mvdf5<-data.frame(Malignancy = "", CHIP = "Lymphoid", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
BL<-data.frame(Malignancy = "", CHIP = "No mCA", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5)
tdf$Total<-as.numeric(table(tukb1$CH)[c("Control","Myeloid","Lymphoid")])
tdf$Positive<-as.numeric(table(tukb1[tukb1$Status ==1,]$CH)[c("Control","Myeloid","Lymphoid")])

tdf$CHIP<-factor(tdf$CHIP,levels=c("L-mCA","M-mCA","No mCA"))

# #tdf[tdf == Inf]<-0

p97 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,color="black") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		#panel.grid.major.x=element_line(color="gray80",size=0.2),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(ctList[1],ctList[2],"black"))

forestdf1<-tdf[,c("CHIP","HR","CI_low","CI_high")]
forestdf1$Malignancy<-rep("Myeloid")


tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1g<-p96
F1h<-ggarrange(ggarrange(p97,ptdf,ncol=2,nrow=1,widths=c(1,1.5)),
	nrow=2,heights=c(2,1))


#=====================================================================
#=====================================================================

## Myeloid

tukb1<-ukbc
tukb1$has_disease<-tukb1$Lymphoid
tukb1$Censor_date<-tukb1$Lymphoid_date

tukb1<-prepCox(tukb1)


ST<-tukb$Status
YR<-tukb$Year
names(ST)<-names(YR)<-as.character(tukb$eid)
tmyeloid$Status_LY<-ST[as.character(tmyeloid$eid)]
tmyeloid$Year_LY<-YR[as.character(tmyeloid$eid)]
write.table(tmyeloid,file = paste(figDir,"scripts/Final_figs/figure_df/F1_ef_incidence.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


# #####
YLAB<-"Cumulative incidence of\nlymphoid malignancies"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb1)
p86<-ggsurvplot(fit,data=tukb1,fun="event",
		legend.title = "mCA", legend=c(0.2,0.75), legend.labs=c("Control","M-mCA","L-mCA"),
		palette=c("black",ctList[2],ctList[1]),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")


## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1[tukb1$CH != "Lymphoid",])
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1[tukb1$CH != "Myeloid",])


## Tables
mvdf1<-extractCox(cox1,"Lymphoid","CHMyeloid","M-mCA","Any")
mvdf5<-extractCox(cox2,"Lymphoid","CHLymphoid","L-mCA","Any")
BL<-data.frame(Malignancy = "", CHIP = "No mCA", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5)
tdf$Total<-as.numeric(table(tukb1$CH)[c("Control","Myeloid","Lymphoid")])
tdf$Positive<-as.numeric(table(tukb1[tukb1$Status ==1,]$CH)[c("Control","Myeloid","Lymphoid")])

tdf$CHIP<-factor(tdf$CHIP,levels=c("L-mCA","M-mCA","No mCA"))


p87 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,color="black") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		#panel.grid.major.x=element_line(color="gray80",size=0.2),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(ctList[1],ctList[2],"black"))

forestdf2<-tdf[,c("CHIP","HR","CI_low","CI_high")]
forestdf2$Malignancy<-rep("Lymphoid")

##

tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1i<-p86
F1j<-ggarrange(ggarrange(p87,ptdf,ncol=2,nrow=1,widths=c(1,1.5)),
	nrow=2,heights=c(2,1))


t<-rbind(forestdf1,forestdf2)
write.table(t,file = paste(figDir,"scripts/Final_figs/figure_df/F1_ef_forestPlots.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


#==========================================================================
## Final figure
## Max A4 size is 8x11

F1<-ggarrange(ggarrange(p0,ggarrange(F1a,F1b,ncol=1,nrow=2),ncol=2,nrow=1,widths=c(0.8,1.2)),
	ggarrange(ggarrange(F1c$plot,F1d,ncol=1,nrow=2,heights=c(2,1.5)),
		ggarrange(F1e$plot,F1f,ncol=1,nrow=2,heights=c(2,1.5)),
		ncol=2,nrow=1),
	ggarrange(ggarrange(F1g$plot,F1h,ncol=1,nrow=2,heights=c(2,1.5)),
		ggarrange(F1i$plot,F1j,ncol=1,nrow=2,heights=c(2,1.5)),
		ncol=2,nrow=1),
	nrow=3,heights=c(2,2,2))



pdf(paste(figDir,"scripts/Final_figs/F1.pdf",sep=""),width=10,height=12)
print (F1)
dev.off()

}
