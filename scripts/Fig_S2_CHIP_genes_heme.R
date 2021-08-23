

library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)
library(ggrepel)

ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"

##===================================================================
## CHIP

chip<-data.frame(fread(paste(figDir,"1_CHIP_UKB.txt",sep="")))
mychip<-data.frame(fread(paste(figDir,"1_MCHIP_UKB.txt",sep="")))

## Myeloid malignancy df
tukb<-data.frame(fread(paste(figDir,"1_UKB_data.txt",sep="")))
tukb$CHIP<-as.factor(tukb$CHIP)
tukb$CHIP_size<-as.factor(tukb$CHIP_size)
tukb$MCHIP<-as.factor(tukb$MCHIP)
tukb$Msize<-as.factor(tukb$Msize)


## genelists
mychip<-mychip[mychip$eid %in% tukb$eid,]
chip<-chip[chip$eid %in% tukb$eid,]
T1<-table(mychip$Hugo_Symbol)
T1<-T1[T1>5]
T2<-table(chip$Hugo_Symbol)
T2<-T2[T2>5]


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

forcedCensorDate<-"2020/03/31"
## ICD list
lyICD<-as.character(hicd[hicd$Group %in% c("Lymphoid","Plasma_cell"),]$ICD)
myICD<-as.character(hicd[hicd$Group %in% c("Myeloid"),]$ICD)
icdList<-c(lyICD,myICD)

tdf1<-getDiag(tukb,icdList,any10,date10)
tdf2<-getCancer(tukb,icdList,cancer10,cancerDate10)
tempdf<-rbind(tdf1,tdf2)
tempdf$Diff2<-as.numeric(difftime(forcedCensorDate,tempdf$Date,units=c("days")))
tempdf<-tempdf[tempdf$Diff2 >= 0,]
tdf<-getUnique(tempdf[,colnames(tdf1)])


tukb$MY_disease<-rep(0,nrow(tukb))
tukb[tukb$eid %in% tdf[tdf$ICD %in% myICD,]$eid,]$MY_disease<-1

tukb$LY_disease<-rep(0,nrow(tukb))
tukb[tukb$eid %in% tdf[tdf$ICD %in% lyICD,]$eid,]$LY_disease<-1

## Age
tukb$Age<-tukb$age_attended_assessment_center
tukb$Age_sq<-tukb$Age*tukb$Age



K<-1
for (gene in c(names(T1),names(T2))){
	t1<-tukb
	t1$CH<-rep(0,nrow(t1))
	if (gene %in% names(T1)){
		eidlist<-mychip[mychip$Hugo_Symbol == gene,]$eid
		t1[t1$eid %in% eidlist,]$CH<-1
	}else{
		eidlist<-chip[chip$Hugo_Symbol == gene,]$eid
		t1[t1$eid %in% eidlist,]$CH<-1
	}
	if (nrow(t1[t1$CH == 1 & t1$MY_disease == 1,]) >= 3){
		t1$CH<-as.factor(t1$CH)
	}else{
		next
	}
	## Model
	m1<-glm(MY_disease ~ Age + Age_sq + ever_smoked + Sex + CH, data=t1,family="binomial")
	## 
	x<-summary(m1)$coefficients
	Est<-x[c("CH1"),1]
	Error<-x[c("CH1"),2]
	P<-x[c("CH1"),4]
	## 
	tempdf<-data.frame(Malignancy = "Myeloid",
		Gene = gene,
		N = length(eidlist),
		N_dis = nrow(t1[t1$eid %in% eidlist & t1$MY_disease == 1,]),
		Estimate = Est,
		Error = Error,
		P = P)
	if (K == 1){
		tdf<-tempdf
		K<-2
	}else{
		tdf<-rbind(tdf,tempdf)
	}
}

mydf<-tdf



### Lymphoid


K<-1
for (gene in c(names(T1),names(T2))){
	t1<-tukb
	t1$CH<-rep(0,nrow(t1))
	if (gene %in% names(T1)){
		eidlist<-mychip[mychip$Hugo_Symbol == gene,]$eid
		t1[t1$eid %in% eidlist,]$CH<-1
	}else{
		eidlist<-chip[chip$Hugo_Symbol == gene,]$eid
		t1[t1$eid %in% eidlist,]$CH<-1
	}
	if (nrow(t1[t1$CH == 1 & t1$LY_disease == 1,]) >= 3){
		t1$CH<-as.factor(t1$CH)
	}else{
		next
	}
	## Model
	m1<-glm(LY_disease ~ Age + Age_sq + ever_smoked + Sex + CH, data=t1,family="binomial")
	## 
	x<-summary(m1)$coefficients
	Est<-x[c("CH1"),1]
	Error<-x[c("CH1"),2]
	P<-x[c("CH1"),4]
	## 
	tempdf<-data.frame(Malignancy = "Lymphoid",
		Gene = gene,
		N = length(eidlist),
		N_dis = nrow(t1[t1$eid %in% eidlist & t1$LY_disease == 1,]),
		Estimate = Est,
		Error = Error,
		P = P)
	if (K == 1){
		tdf<-tempdf
		K<-2
	}else{
		tdf<-rbind(tdf,tempdf)
	}
}


mydf<-rbind(mydf,tdf)
mydf$q<-p.adjust(mydf$P,method="fdr")
mydf$Malignancy<-factor(mydf$Malignancy,levels=c("Myeloid","Lymphoid"))

## Odds ratio
mydf$Odds<-exp(mydf$Estimate)
mydf$Upper_ci<-round(exp(mydf$Estimate+1.96*mydf$Error),3)
mydf$Lower_ci<-round(exp(mydf$Estimate-1.96*mydf$Error),3)

p1 <- ggplot(data=mydf, aes(x=-log10(q), y=Odds, ymin=Lower_ci, ymax=Upper_ci, color=Malignancy)) +
	geom_hline(yintercept=1, col="black") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	coord_flip() + 
	xlab("-log10(adjusted p-value)") + ylab("Odds ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid.major.x=element_line(color="grey"),
		axis.text=element_text(color="black"),
		axis.ticks=element_line(color="black"),
		legend.position="none")+
	scale_color_manual(values=c(ctList[2],ctList[1]))+
	geom_text_repel(data = mydf,aes(label=Gene))+
	#geom_text_repel(data = subset(mydf, q < 0.01),aes(label=Gene))+
	scale_y_log10()


P<-ggarrange(p1)
pdf(paste(figDir,"scripts/Final_figs/Fig_S2_forest.pdf",sep=""),width=6,height=5)
print (P)
dev.off()

write.table(mydf,file=paste(figDir,"scripts/Final_figs/Fig_S2_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


