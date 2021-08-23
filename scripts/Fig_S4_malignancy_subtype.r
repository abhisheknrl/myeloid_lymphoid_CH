
## Histogram of frequencies of the malignancies
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

if (CMD %in% c("WESL","WESM")){


## function to extract a counts df
get_df<-function(tmp,icdlist,Disease){
	tmp$LMCHIP<-rep("Control",nrow(tmp))
	tmp[tmp$CHIP == 1 & tmp$MCHIP == 1,]$LMCHIP<-"ML"
	tmp[tmp$CHIP == 1 & tmp$MCHIP == 0,]$LMCHIP<-"Lymphoid"
	tmp[tmp$CHIP == 0 & tmp$MCHIP == 1,]$LMCHIP<-"Myeloid"
	##
	tmp1<-data.frame(Disease = Disease,
		Total_control = nrow(tmp[tmp$LMCHIP == "Control",]),
		Total_myeloid = nrow(tmp[tmp$LMCHIP == "Myeloid",]),
		Total_lymphoid = nrow(tmp[tmp$LMCHIP == "Lymphoid",]),
		Total_ML = nrow(tmp[tmp$LMCHIP == "ML",]),
		Disease_control = nrow(tmp[tmp$LMCHIP == "Control" & tmp$ICD %in% icdlist,]),
		Disease_myeloid = nrow(tmp[tmp$LMCHIP == "Myeloid"  & tmp$ICD %in% icdlist,]),
		Disease_lymphoid = nrow(tmp[tmp$LMCHIP == "Lymphoid"  & tmp$ICD %in% icdlist,]),
		Disease_ML = nrow(tmp[tmp$LMCHIP == "ML"  & tmp$ICD %in% icdlist,]))
	tmp1$Fraction_control<-round(tmp1$Disease_control*100/tmp1$Total_control,3)
	tmp1$Fraction_myeloid<-round(tmp1$Disease_myeloid*100/tmp1$Total_myeloid,3)
	tmp1$Fraction_lymphoid<-round(tmp1$Disease_lymphoid*100/tmp1$Total_lymphoid,3)
	tmp1$Fraction_ML<-round(tmp1$Disease_ML*100/tmp1$Total_ML,3)
	return (tmp1)
}

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

## Caucasian
ukbc$Caucasian<-rep(1,nrow(ukbc))
ukbc[is.na(ukbc$genetic_ethnic_grouping),]$Caucasian<-0

#=====================================================================
## Enrichment of specific lymphoid malignancies
DOD<-c()
death2<-c()
forcedCensorDate<-"2020-03-31"
deathCensorDate<-"2020-03-31"



if (CMD == "WESL"){

## Lymphoid ICD list
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


## Myeloid icd list
icdList<-as.character(hicd[hicd$Group %in% c("Myeloid"),]$ICD)
mds_icdList<-c(paste("D46",0:9,sep=""),"C946")
aml_icdList<-paste("C92",0,sep="")
mpn_icdList<-c("D752","D45","C945","D471","D473")
other_myeloid<-icdList[!icdList %in% c(mds_icdList,aml_icdList,mpn_icdList)]



icdList<-as.character(hicd[hicd$Group %in% c("Lymphoid","Plasma_cell","Myeloid"),]$ICD)
## Lymphoid malignancies
tdf1<-getDiag(ukbc,icdList,any10,date10)
tdf2<-getCancer(ukbc,icdList,cancer10,cancerDate10)
tempdf<-rbind(tdf1,tdf2)
tempdf$Diff2<-as.numeric(difftime(forcedCensorDate,tempdf$Date,units=c("days")))
tempdf<-tempdf[tempdf$Diff2 >= 0,]
tdf<-getUnique(tempdf[,colnames(tdf1)])


## For the unspecified types
# t<-getUnique(tempdf[!tempdf$ICD %in% c("C859","C851","C857"),colnames(tdf1)])
# t1<-tdf[!tdf$eid %in% t$eid,colnames(t)]
# tdf<-rbind(t,t1)


tempicd<-as.character(tdf$ICD)
names(tempicd)<-as.character(tdf$eid)
ukbc$ICD<-tempicd[as.character(ukbc$eid)]
ukbc[is.na(ukbc$ICD),]$ICD<-""

## Data
tdf<-rbind(get_df(ukbc,PC,"MM_MGUS"),
	get_df(ukbc,CLL,"CLL_SLL"),
	get_df(ukbc,CIRC,"Circulating_lymphoma"),
	get_df(ukbc,NHL,"Unspecified_NHL"),
	get_df(ukbc,DLBCL,"DLBCL"),
	get_df(ukbc,FL,"FL"),
	get_df(ukbc,WS,"WM"),
	get_df(ukbc,HL,"HL"),
	get_df(ukbc,nhl,"Other_lymphoma"),
	get_df(ukbc,mds_icdList,"MDS"),
	get_df(ukbc,aml_icdList,"AML"),
	get_df(ukbc,mpn_icdList,"MPN"),
	get_df(ukbc,other_myeloid,"Other_myeloid"))

tdf$Total_disease<-tdf$Disease_control + tdf$Disease_myeloid + tdf$Disease_lymphoid + tdf$Disease_ML

write.table(tdf,file=paste(figDir,"scripts/data/WES_malig_histogram_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

}

}

#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
## 500K cases

if (CMD %in% c("L500","M500")){

get_df<-function(tmp,icdlist,Disease){
	tmp$LMCHIP<-rep("Control",nrow(tmp))
	tmp[tmp$LCNV == 1 & tmp$MCNV == 1,]$LMCHIP<-"ML"
	tmp[tmp$LCNV == 1 & tmp$MCNV == 0,]$LMCHIP<-"Lymphoid"
	tmp[tmp$LCNV == 0 & tmp$MCNV == 1,]$LMCHIP<-"Myeloid"
	tmp[tmp$LCNV == 0 & tmp$MCNV == 0 & tmp$CNV == 1,]$LMCHIP<-"Unk"
	##
	tmp1<-data.frame(Disease = Disease,
		Total_control = nrow(tmp[tmp$LMCHIP == "Control",]),
		Total_myeloid = nrow(tmp[tmp$LMCHIP == "Myeloid",]),
		Total_lymphoid = nrow(tmp[tmp$LMCHIP == "Lymphoid",]),
		Total_ML = nrow(tmp[tmp$LMCHIP == "ML",]),
		Total_Unk = nrow(tmp[tmp$LMCHIP == "Unk",]),
		Disease_control = nrow(tmp[tmp$LMCHIP == "Control" & tmp$ICD %in% icdlist,]),
		Disease_myeloid = nrow(tmp[tmp$LMCHIP == "Myeloid"  & tmp$ICD %in% icdlist,]),
		Disease_lymphoid = nrow(tmp[tmp$LMCHIP == "Lymphoid"  & tmp$ICD %in% icdlist,]),
		Disease_ML = nrow(tmp[tmp$LMCHIP == "ML"  & tmp$ICD %in% icdlist,]),
		Disease_Unk = nrow(tmp[tmp$LMCHIP == "Unk"  & tmp$ICD %in% icdlist,]))
	tmp1$Fraction_control<-round(tmp1$Disease_control*100/tmp1$Total_control,3)
	tmp1$Fraction_myeloid<-round(tmp1$Disease_myeloid*100/tmp1$Total_myeloid,3)
	tmp1$Fraction_lymphoid<-round(tmp1$Disease_lymphoid*100/tmp1$Total_lymphoid,3)
	tmp1$Fraction_ML<-round(tmp1$Disease_ML*100/tmp1$Total_ML,3)
	tmp1$Fraction_Unk<-round(tmp1$Disease_Unk*100/tmp1$Total_Unk,3)
	return (tmp1)
}

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


#=====================================================================
## Enrichment of specific lymphoid malignancies

DOD<-c()
death2<-c()
forcedCensorDate<-"2020-03-31"
deathCensorDate<-"2020-03-31"


## Lymphoid malignancies
if (CMD == "L500"){

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


icdList<-as.character(hicd[hicd$Group %in% c("Lymphoid","Plasma_cell","Myeloid"),]$ICD)
## Lymphoid malignancies
tdf1<-getDiag(ukbc,icdList,any10,date10)
tdf2<-getCancer(ukbc,icdList,cancer10,cancerDate10)
tempdf<-rbind(tdf1,tdf2)
tempdf$Diff2<-as.numeric(difftime(forcedCensorDate,tempdf$Date,units=c("days")))
tempdf<-tempdf[tempdf$Diff2 >= 0,]
tdf<-getUnique(tempdf[,colnames(tdf1)])

##
tempicd<-as.character(tdf$ICD)
names(tempicd)<-as.character(tdf$eid)
ukbc$ICD<-tempicd[as.character(ukbc$eid)]
ukbc[is.na(ukbc$ICD),]$ICD<-""

## Data
tdf<-rbind(get_df(ukbc,PC,"MM_MGUS"),
	get_df(ukbc,CLL,"CLL_SLL"),
	get_df(ukbc,CIRC,"Circulating_lymphoma"),
	get_df(ukbc,NHL,"Unspecified_NHL"),
	get_df(ukbc,DLBCL,"DLBCL"),
	get_df(ukbc,FL,"FL"),
	get_df(ukbc,WS,"WM"),
	get_df(ukbc,HL,"HL"),
	get_df(ukbc,nhl,"Other_lymphoma"),
	get_df(ukbc,mds_icdList,"MDS"),
	get_df(ukbc,aml_icdList,"AML"),
	get_df(ukbc,mpn_icdList,"MPN"),
	get_df(ukbc,other_myeloid,"Other_myeloid"))

tdf$Total_disease<-tdf$Disease_control + tdf$Disease_myeloid + tdf$Disease_lymphoid + tdf$Disease_ML

write.table(tdf,file=paste(figDir,"scripts/data/UKB_500_malig_histogram_df.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

}

}


#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
## Figure


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

tdf<-data.frame(fread(paste(figDir,"scripts/data/WES_malig_histogram_df.txt",sep="")))
dlist<-as.character(tdf$Disease)
tdf1<-data.frame(Disease=c(dlist, dlist, dlist, dlist),
	N = c(tdf$Disease_control, tdf$Disease_myeloid, tdf$Disease_lymphoid, tdf$Disease_ML),
	Fraction = c(tdf$Fraction_control, tdf$Fraction_myeloid, tdf$Fraction_lymphoid, tdf$Fraction_ML),
	CHIP = c(rep("No CHIP",nrow(tdf)), rep("M-CHIP",nrow(tdf)), rep("L-CHIP",nrow(tdf)), rep("M+L",nrow(tdf))))

tdf1<-tdf1[!tdf1$Disease %in% c("MDS","MPN","AML","Other_myeloid"),]

tdf1$Disease<-gsub("_","\n",as.character(tdf1$Disease))
tdf1$Disease<-gsub("MM","MM/",as.character(tdf1$Disease))
tdf1$Disease<-gsub("CLL","CLL/",as.character(tdf1$Disease))
tdf1$CHIP<-factor(tdf1$CHIP,levels=c("No CHIP","M-CHIP","L-CHIP","M+L"))
tdf1$Disease<-factor(tdf1$Disease,levels=rev(unique(as.character(tdf1$Disease))))


temp<-tdf1[tdf1$CHIP != "M+L",]
p1<-ggplot(temp,aes(Fraction,Disease, Group=Disease,fill=CHIP))+
	geom_bar(stat="identity",position="dodge")+
	geom_text(aes(label=N), hjust=-0.2, size=3, position = position_dodge(0.9))+
	theme_bw()+
	theme(panel.grid.major.y=element_blank(),
		axis.ticks.y = element_blank(),
		panel.grid.minor=element_blank(),
		legend.position = c(0.8, 0.95),
		legend.key.size = unit(0.3, 'cm'),
		legend.title = element_blank(),
		axis.text=element_text(color="black"))+
	scale_fill_manual(values=c("black",ctList[2],ctList[1],"orange"))+
	ylab("")+
	xlim(0,3)+
	xlab("Percent with incident malignancies")


#=====================================================================
## Histogram myeloid

tdf1<-data.frame(Disease=c(dlist, dlist, dlist, dlist),
	N = c(tdf$Disease_control, tdf$Disease_myeloid, tdf$Disease_lymphoid, tdf$Disease_ML),
	Fraction = c(tdf$Fraction_control, tdf$Fraction_myeloid, tdf$Fraction_lymphoid, tdf$Fraction_ML),
	CHIP = c(rep("No CHIP",nrow(tdf)), rep("M-CHIP",nrow(tdf)), rep("L-CHIP",nrow(tdf)), rep("M+L",nrow(tdf))))

tdf1<-tdf1[tdf1$Disease %in% c("MDS","MPN","AML","Other_myeloid"),]

tdf1$Disease<-gsub("_","\n",as.character(tdf1$Disease))
tdf1$CHIP<-factor(tdf1$CHIP,levels=c("No CHIP","M-CHIP","L-CHIP","M+L"))
tdf1$Disease<-factor(tdf1$Disease,levels=c("Other\nmyeloid","MPN","MDS","AML"))

temp<-tdf1[tdf1$CHIP != "M+L",]
p2<-ggplot(temp,aes(Fraction,Disease,Group=Disease,fill=CHIP))+
	geom_bar(stat="identity",position="dodge")+
	geom_text(aes(label=N), hjust=-0.2, size=3, position = position_dodge(0.9))+
	theme_bw()+
	theme(panel.grid.major.y=element_blank(),
		axis.ticks.y=element_blank(),
		panel.grid.minor=element_blank(),
		legend.position = c(0.8, 0.88),
		legend.key.size = unit(0.3, 'cm'),
		legend.title = element_blank(),
		axis.text=element_text(color="black"))+
	scale_fill_manual(values=c("black",ctList[2],ctList[1],"orange"))+
	ylab("")+
	xlim(0,1)+
	xlab("Percent with incident malignancies")



#=====================================================================
## 500K WES
tdf<-data.frame(fread(paste(figDir,"scripts/data/UKB_500_malig_histogram_df.txt",sep="")))
dlist<-as.character(tdf$Disease)
tdf1<-data.frame(Disease=c(dlist,dlist,dlist,dlist,dlist),
	N = c(tdf$Disease_control,tdf$Disease_myeloid,tdf$Disease_lymphoid,tdf$Disease_ML,tdf$Disease_Unk),
	Fraction = c(tdf$Fraction_control,tdf$Fraction_myeloid,tdf$Fraction_lymphoid,tdf$Fraction_ML,tdf$Fraction_Unk),
	CHIP = c(rep("No mCA",nrow(tdf)),rep("M-mCA",nrow(tdf)),rep("L-mCA",nrow(tdf)),rep("A-mCA",nrow(tdf)),rep("U-mCA",nrow(tdf))))

tdf1<-tdf1[!tdf1$Disease %in% c("MDS","MPN","AML","Other_myeloid"),]

tdf1$Disease<-gsub("_","\n",as.character(tdf1$Disease))
tdf1$CHIP<-factor(tdf1$CHIP,levels=c("No mCA","M-mCA","L-mCA","A-mCA","U-mCA"))
tdf1$Disease<-factor(tdf1$Disease,levels=rev(unique(as.character(tdf1$Disease))))

p3<-ggplot(tdf1,aes(Fraction,Disease,Group=Disease,fill=CHIP))+
	geom_bar(stat="identity",position="dodge")+
	geom_text(aes(label=N), hjust=-0.2, size=3, position = position_dodge(0.9))+
	theme_bw()+
	theme(panel.grid.major.y=element_blank(),
		axis.ticks.y=element_blank(),
		panel.grid.minor=element_blank(),
		legend.position = c(0.8, 0.93),
		legend.key.size = unit(0.3, 'cm'),
		legend.title = element_blank(),
		axis.text=element_text(color="black"))+
	scale_fill_manual(values=c("black",ctList[2],ctList[1],"orange","grey"))+
	ylab("")+
	xlim(0,8)+
	xlab("Percent with incident malignancies")


#=====================================================================
## Histogram myeloid
tdf1<-data.frame(Disease=c(dlist,dlist,dlist,dlist,dlist),
	N = c(tdf$Disease_control,tdf$Disease_myeloid,tdf$Disease_lymphoid,tdf$Disease_ML,tdf$Disease_Unk),
	Fraction = c(tdf$Fraction_control,tdf$Fraction_myeloid,tdf$Fraction_lymphoid,tdf$Fraction_ML,tdf$Fraction_Unk),
	CHIP = c(rep("No mCA",nrow(tdf)),rep("M-mCA",nrow(tdf)),rep("L-mCA",nrow(tdf)),rep("A-mCA",nrow(tdf)),rep("U-mCA",nrow(tdf))))

tdf1<-tdf1[tdf1$Disease %in% c("MDS","MPN","AML","Other_myeloid"),]

tdf1$Disease<-gsub("_","\n",as.character(tdf1$Disease))
tdf1$CHIP<-factor(tdf1$CHIP,levels=c("No mCA","M-mCA","L-mCA","A-mCA","U-mCA"))
tdf1$Disease<-factor(tdf1$Disease,levels=c("Other\nmyeloid","MPN","MDS","AML"))

p4<-ggplot(tdf1,aes(Fraction,Disease,Group=Disease,fill=CHIP))+
	geom_bar(stat="identity",position="dodge")+
	geom_text(aes(label=N), hjust=-0.2, size=3, position = position_dodge(0.9))+
	theme_bw()+
	theme(panel.grid.major.y=element_blank(),
		axis.ticks.y=element_blank(),
		panel.grid.minor=element_blank(),
		legend.position = c(0.8, 0.83),
		legend.key.size = unit(0.3, 'cm'),
		legend.title = element_blank(),
		axis.text=element_text(color="black"))+
	scale_fill_manual(values=c("black",ctList[2],ctList[1],"orange","grey"))+
	ylab("")+
	xlim(0,8)+
	xlab("Percent with incident malignancies")


F1a<-ggarrange(ggarrange(p2,p1,ncol=1,nrow=2,heights=c(4,9)),
	ggarrange(p4,p3,ncol=1,nrow=2,heights=c(4,9)),
	nrow=1,ncol=2)


pdf(paste(figDir,"scripts/Final_figs/Fig_S4.pdf",sep=""),width=8,height=10)
print (F1a)
dev.off()


}

