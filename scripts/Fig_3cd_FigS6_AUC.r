
## AUC, PRED

args = commandArgs(trailingOnly=TRUE)
CMD<-args[1]


library(ggplot2)
library(data.table)
library(colortools)
library(survival)
library(survminer)
library(ggpubr)

ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"


forcedCensorDate<-as.Date("2020-03-31")

#=====================================================================
## Reload data

chip<-read.table(paste(figDir,"1_CHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
mychip<-read.table(paste(figDir,"1_MCHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
cnv<-data.frame(fread(paste(figDir,"data/UKB_500K_CNV.txt",sep="")))
cnv[is.na(cnv$CELL_FRAC),]$CELL_FRAC<-0.001
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


#=====================================================================
#=====================================================================
## Functions

## Blood count thresholds
ALC_upper<-4.25
ALC_lower<-0.65
WBC_upper<-9.57
WBC_lower<-3.53
RBC_upper<-5.5
RBC_lower<-3.96
PLT_upper<-397.1
PLT_lower<-169.06
ANC_upper<-7.06
ANC_lower<-1.47

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

#=====================================================================
## Load data

ukbc<-data.frame(fread(paste(figDir,"Fig_1_WES_malignancy_df.txt",sep="")))

ukbc$Caucasian<-as.factor(ukbc$Caucasian)
ukbc$ever_smoked<-factor(ukbc$ever_smoked,levels=c("No","Yes"))
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))


#=====================================================================
## CBC stratification

library(pROC)

if (CMD == "MY"){

ukbc1<-ukbc[(!is.na(ukbc$red_blood_cell_count)) & (!is.na(ukbc$platelet_count)) & (!is.na(ukbc$neutrophil_count)) & (!is.na(ukbc$lymphocyte_count)),]
ukbc11<-ukbc1[ukbc1$red_blood_cell_count < RBC_lower | ukbc1$platelet_count < PLT_lower | ukbc1$neutrophil_count < ANC_lower,]
ukbc12<-ukbc1[ukbc1$red_blood_cell_count > RBC_upper | ukbc1$platelet_count > PLT_upper | ukbc1$neutrophil_count > ANC_upper,]
ukbc12<-ukbc12[!ukbc12$eid %in% ukbc11$eid,]

ukbc1<-ukbc1[!ukbc1$eid %in% c(ukbc11$eid,ukbc12$eid),]


## UKBC
ukbc1$BCI<-rep("Normal",nrow(ukbc1))
ukbc11$BCI<-rep("Plt_low",nrow(ukbc11))
ukbc12$BCI<-rep("Plt_high",nrow(ukbc12))


## UKBC
ukbc1<-rbind(ukbc1,ukbc11,ukbc12)
ukbc1$BCI<-factor(ukbc1$BCI,levels=c("Normal","Plt_low","Plt_high"))


#=====================================================================
## Model

temp<-ukbc1
temp$has_disease<-temp$Myeloid
temp$Censor_date<-temp$Myeloid_date
temp<-prepCox(temp)

## Age, Sex, EverSmoke, Caucasian, lymphocyte count, RBC count, Platelet count, nCH, VAF CH
CN<-c("eid","Age","Age_categ","Sex","EverSmoke","Caucasian","lymphocyte_count","red_blood_cell_count","platelet_count","neutrophil_count","Myeloid","Lymphoid","CLL","Year","BCI")
tukb<-temp[,CN]
tukb$Myeloid<-as.factor(tukb$Myeloid)

#=====================================================================
## Add CHIP and mCA Data

Tchip<-table(mychip$eid)
Tcnv<-table(cnv[cnv$MCNV == 1 & cnv$eid %in% tukb$eid,]$eid)
Tch<-table(c(mychip$eid,cnv[cnv$MCNV == 1,]$eid))


## For individuals with multiple CH alterations
cvaf<-c()
mvaf<-c()
mcvaf<-c()
for (EID in names(Tch)){
	## CHIP
	if (EID %in% mychip$eid){
		v1<-mychip[mychip$eid == EID,]$VAF
	}else{
		v1<-c(0)
	}
	## mCA
	if (EID %in% cnv$eid){
		v2<-cnv[cnv$eid == EID & cnv$MCNV == 1,]$CELL_FRAC
	}else{
		v2<-c(0)
	}
	##
	cvaf[EID]<-max(v1)
	mvaf[EID]<-max(v2)
	mcvaf[EID]<-max(c(v1,v2))
}

## Adding VAf
tukb$CHIP_VAF<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(cvaf[cvaf > 0 & cvaf < 0.1]),]$CHIP_VAF<-1
tukb[tukb$eid %in% names(cvaf[cvaf >= 0.1]),]$CHIP_VAF<-2
tukb$CHIP_VAF<-as.factor(tukb$CHIP_VAF)

## Mca
tukb$mCA_VAF<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(mvaf[mvaf > 0 & mvaf < 0.1]),]$mCA_VAF<-1
tukb[tukb$eid %in% names(mvaf[mvaf >= 0.1]),]$mCA_VAF<-2
tukb$mCA_VAF<-as.factor(tukb$mCA_VAF)

## Combined
tukb$CH_VAF<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(mcvaf[mcvaf > 0 & mcvaf < 0.1]),]$CH_VAF<-1
tukb[tukb$eid %in% names(mcvaf[mcvaf >= 0.1]),]$CH_VAF<-2
tukb$CH_VAF<-as.factor(tukb$CH_VAF)

## Number of variants
tukb$CHIP_N<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(Tchip[Tchip == 1]),]$CHIP_N<-1
tukb[tukb$eid %in% names(Tchip[Tchip > 1]),]$CHIP_N<-2
tukb$CHIP_N<-as.factor(tukb$CHIP_N)

## mCA
tukb$mCA_N<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(Tcnv[Tcnv == 1]),]$mCA_N<-1
tukb[tukb$eid %in% names(Tcnv[Tcnv > 1]),]$mCA_N<-2
tukb$mCA_N<-as.factor(tukb$mCA_N)

## Combined
tukb$CH_N<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(Tch[Tch == 1]),]$CH_N<-1
tukb[tukb$eid %in% names(Tch[Tch > 1]),]$CH_N<-2
tukb$CH_N<-as.factor(tukb$CH_N)


## CHIP and VAF features
tukb$CHIP_var<-rep(0,nrow(tukb))
tukb[tukb$CHIP_N == 1 & tukb$CHIP_VAF == 1,]$CHIP_var<-1
tukb[tukb$CHIP_N == 1 & tukb$CHIP_VAF == 2,]$CHIP_var<-2
tukb[tukb$CHIP_N == 2,]$CHIP_var<-3
tukb$CHIP_var<-as.factor(tukb$CHIP_var)

## mCA
tukb$mCA_var<-rep(0,nrow(tukb))
tukb[tukb$mCA_N == 1 & tukb$mCA_VAF == 1,]$mCA_var<-1
tukb[tukb$mCA_N == 1 & tukb$mCA_VAF == 2,]$mCA_var<-2
tukb[tukb$mCA_N == 2,]$mCA_var<-3
tukb$mCA_var<-as.factor(tukb$mCA_var)

## Combined
tukb$CH_var<-rep(0,nrow(tukb))
tukb[tukb$CH_N == 1 & tukb$CH_VAF == 1,]$CH_var<-1
tukb[tukb$CH_N == 1 & tukb$CH_VAF == 2,]$CH_var<-2
tukb[tukb$CH_N == 2,]$CH_var<-3
tukb$CH_var<-as.factor(tukb$CH_var)


## MCI
tukb$MCI<-rep("Normal",nrow(tukb))
tukb[tukb$BCI == "Plt_high",]$MCI<-"High"
tukb[tukb$BCI == "Plt_low",]$MCI<-"Low"
tukb$MCI<-factor(tukb$MCI,levels=c("Normal","Low","High"))


#=====================================================================
## Ten fold cross validation


## Positive and negative data
t1<-unique(tukb[tukb$Myeloid == 1,])
t2<-unique(tukb[tukb$Myeloid == 0 & tukb$Lymphoid == 0,])



## sample
Pmod1<-Pmod2<-Pmod3<-Pmod4<-Pmod5<-Pmod6<-Pmod7<-Pmod8<-Pmod9<-Pmod10<-c()
PT1<-PT2<-PT3<-PT4<-PT5<-PT6<-PT7<-PT8<-PT9<-PT10<-Lab<-c()


S11<-sample(1:nrow(t1),nrow(t1),replace=FALSE)
S12<-sample(1:nrow(t2),nrow(t2),replace=FALSE)
N11<-round(nrow(t1)/10)
N12<-round(nrow(t2)/10)


for (j in 1:10){
	## Subset data
	ST11<-(j-1)*N11+1
	ST12<-(j-1)*N12+1
	EN11<-j*N11
	EN12<-j*N12
	if (j == 10){
		EN11<-nrow(t1)
		EN12<-nrow(t2)
	}
	## Pull training and testing data
	## Test data
	t11<-t1[S11[ST11:EN11],]
	t12<-t2[S12[ST12:EN12],]
	temp1<-rbind(t11,t12)
	Lab<-c(Lab,temp1$Myeloid)
	## Training data
	temp2<-rbind(t1,t2)
	temp2<-temp2[!temp2$eid %in% temp1$eid,]
	#================================================================
	# Train and predict
	## Demographic
	mod1<-glm(Myeloid ~ Age_categ + Sex + EverSmoke + Caucasian, data =  temp2,family="binomial")
	PT1<-c(PT1,predict(mod1,temp1,type="response"))
	# Add CHIP
	mod2<-glm(Myeloid ~ Age_categ + Sex + EverSmoke + Caucasian + CHIP_N + CHIP_VAF, data =  temp2,family="binomial")
	PT2<-c(PT2,predict(mod2,temp1,type="response"))
	## Add mCA only
	mod3<-glm(Myeloid ~ Age_categ + Sex + EverSmoke + Caucasian + mCA_N + mCA_VAF, data = temp2,family="binomial")
	PT3<-c(PT3,predict(mod3,temp1,type="response"))
	## Add both as separate variables
	mod4<-glm(Myeloid ~ Age_categ + Sex + EverSmoke + Caucasian + CHIP_N + CHIP_VAF + mCA_N + mCA_VAF, data = temp2,family="binomial")
	PT4<-c(PT4,predict(mod4,temp1,type="response"))
	## Add both as combined variable
	#mod5<-glm(Myeloid ~ Age_categ + Sex + EverSmoke + Caucasian + CH_var, data = temp2,family="binomial")
	#PT5<-c(PT5,predict(mod5,temp1,type="response"))
	## Add CBC	
	mod6<-glm(Myeloid ~ Age_categ + Sex + EverSmoke + Caucasian + red_blood_cell_count + platelet_count + neutrophil_count + MCI, data = temp2,family="binomial")
	PT6<-c(PT6,predict(mod6,temp1,type="response"))
	## Adding CBC + CHIP
	mod7<-glm(Myeloid ~ Age_categ + Sex + EverSmoke + Caucasian + red_blood_cell_count + platelet_count + neutrophil_count + MCI + CHIP_N + CHIP_VAF, data = temp2,family="binomial")
	PT7<-c(PT7,predict(mod7,temp1,type="response"))
	## Adding CBC + mCA
	mod8<-glm(Myeloid ~ Age_categ + Sex + EverSmoke + Caucasian + red_blood_cell_count + platelet_count + neutrophil_count + MCI + mCA_N + mCA_VAF, data = temp2,family="binomial")
	PT8<-c(PT8,predict(mod8,temp1,type="response"))
	## Adding CBC + CHIP + mCA
	mod9<-glm(Myeloid ~ Age_categ + Sex + EverSmoke + Caucasian + red_blood_cell_count + platelet_count + neutrophil_count + MCI + CHIP_N + CHIP_VAF + mCA_N + mCA_VAF, data = temp2,family="binomial")
	PT9<-c(PT9,predict(mod9,temp1,type="response"))
}


##
P1<-roc(Lab,PT1,plot=TRUE)
P2<-roc(Lab,PT2,plot=TRUE)
P3<-roc(Lab,PT3,plot=TRUE)
P4<-roc(Lab,PT4,plot=TRUE)
#P5<-roc(Lab,PT5,plot=TRUE)
P6<-roc(Lab,PT6,plot=TRUE)
P7<-roc(Lab,PT7,plot=TRUE)
P8<-roc(Lab,PT8,plot=TRUE)
P9<-roc(Lab,PT9,plot=TRUE)

## plot the AUC
pdf(paste(figDir,"scripts/Final_figs/Fig_S6a_MY.pdf",sep=""),width=9,height=7)
r1 <- plot(P1, print.auc = TRUE, lty = 1, col = "black", print.auc.y = .05, add = FALSE)
r2 <- plot(P2, print.auc = TRUE, lty = 2, col = "red", print.auc.y = .1, add = TRUE)
r3 <- plot(P3, print.auc = TRUE, lty = 2, col = "green", print.auc.y = .15, add = TRUE)
r4 <- plot(P4, print.auc = TRUE, lty = 2, col = "orange", print.auc.y = .2, add = TRUE)
#r5 <- plot(P5, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .25, add = TRUE)
r6 <- plot(P6, print.auc = TRUE, lty = 1, col = "black", print.auc.y = .25, add = TRUE)
r7 <- plot(P7, print.auc = TRUE, lty = 1, col = "red", print.auc.y = .3, add = TRUE)
r8 <- plot(P8, print.auc = TRUE, lty = 1, col = "green", print.auc.y = .35, add = TRUE)
r9 <- plot(P9, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .4, add = TRUE)
dev.off()


setEPS()
postscript(paste(figDir,"scripts/Final_figs/Fig_S6a_MY.eps",sep=""),width=9,height=7)
r1 <- plot(P1, print.auc = TRUE, lty = 1, col = "grey", print.auc.y = .05, add = FALSE)
r2 <- plot(P2, print.auc = TRUE, lty = 2, col = "red", print.auc.y = .1, add = TRUE)
r3 <- plot(P3, print.auc = TRUE, lty = 2, col = "green", print.auc.y = .15, add = TRUE)
r4 <- plot(P4, print.auc = TRUE, lty = 2, col = "orange", print.auc.y = .2, add = TRUE)
#r5 <- plot(P5, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .25, add = TRUE)
r6 <- plot(P6, print.auc = TRUE, lty = 1, col = "black", print.auc.y = .25, add = TRUE)
r7 <- plot(P7, print.auc = TRUE, lty = 1, col = "red", print.auc.y = .3, add = TRUE)
r8 <- plot(P8, print.auc = TRUE, lty = 1, col = "green", print.auc.y = .35, add = TRUE)
r9 <- plot(P9, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .4, add = TRUE)
dev.off()


pdf(paste(figDir,"scripts/Final_figs/Fig_3c_MY.pdf",sep=""),width=9,height=7)
r1 <- plot(P1, print.auc = TRUE, lty = 1, col = "grey", print.auc.y = .05, add = FALSE)
r6 <- plot(P6, print.auc = TRUE, lty = 1, col = "black", print.auc.y = .25, add = TRUE)
r7 <- plot(P7, print.auc = TRUE, lty = 1, col = "red", print.auc.y = .3, add = TRUE)
r9 <- plot(P9, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .4, add = TRUE)
dev.off()


setEPS()
postscript(paste(figDir,"scripts/Final_figs/Fig_3c_MY.eps",sep=""),width=9,height=7)
r1 <- plot(P1, print.auc = TRUE, lty = 1, col = "grey", print.auc.y = .05, add = FALSE)
r6 <- plot(P6, print.auc = TRUE, lty = 1, col = "black", print.auc.y = .3, add = TRUE)
r7 <- plot(P7, print.auc = TRUE, lty = 1, col = "red", print.auc.y = .35, add = TRUE)
r9 <- plot(P9, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .45, add = TRUE)
dev.off()

myobj<-list(Lab,PT1,PT2,PT3,PT4,PT5,PT6,PT7,PT8,PT9,PT10)
save(myobj,file=paste(figDir,"scripts/Final_figs/Fig_3c_MY_object.RData",sep=""))


## DF for AUC Fig.
rdf<-data.frame(Label = Lab,
	DEM = PT1,
	DEM_M_CHIP = PT2,
	DEM_M_mCA = PT3,
	DEM_M_CHIP_M_mCA = PT4,
	DEM_CBC = PT6,
	DEM_CBC_M_CHIP = PT7,
	DEM_CBC_M_mCA = PT8,
	DEM_CBC_M_CHIP_M_mCA = PT9)

write.table(rdf,file=paste(figDir,"scripts/Final_figs/figure_df/Fig_3c_MY.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

print (warnings())

}
#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
## Model

if (CMD == "LY"){


## Subset
ukbc1<-ukbc[(!is.na(ukbc$red_blood_cell_count)) & (!is.na(ukbc$platelet_count)) & (!is.na(ukbc$neutrophil_count)) & (!is.na(ukbc$lymphocyte_count)),]
ukbc22<-ukbc1[ukbc1$lymphocyte_count > ALC_upper,]

ukbc1<-ukbc1[!ukbc1$eid %in% c(ukbc22$eid),]


## UKBC
ukbc1$BCI<-rep("Normal",nrow(ukbc1))
ukbc22$BCI<-rep("ALC_high",nrow(ukbc22))

## UKBC
ukbc1<-rbind(ukbc1,ukbc22)
ukbc1$BCI<-factor(ukbc1$BCI,levels=c("Normal","ALC_high"))


#==================================================================

temp<-ukbc1
temp$has_disease<-temp$CLL
temp$Censor_date<-temp$CLL_date
temp<-prepCox(temp)

## Age, Sex, EverSmoke, Caucasian, lymphocyte count, RBC count, Platelet count, nCH, VAF CH
CN<-c("eid","Age","Age_categ","Sex","EverSmoke","Caucasian","lymphocyte_count","red_blood_cell_count","platelet_count","neutrophil_count","Myeloid","Lymphoid","CLL","Year","BCI")
tukb<-temp[,CN]
tukb$CLL<-as.factor(tukb$CLL)

#=====================================================================
## Add CHIP and mCA Data

Tchip<-table(chip$eid)
Tcnv<-table(cnv[cnv$LCNV == 1 & cnv$eid %in% tukb$eid,]$eid)
Tch<-table(c(chip$eid,cnv[cnv$LCNV == 1,]$eid))


## For individuals with multiple CH alterations
cvaf<-c()
mvaf<-c()
mcvaf<-c()
for (EID in names(Tch)){
	## CHIP
	if (EID %in% mychip$eid){
		v1<-chip[chip$eid == EID,]$VAF
	}else{
		v1<-c(0)
	}
	## mCA
	if (EID %in% cnv$eid){
		v2<-cnv[cnv$eid == EID & cnv$LCNV == 1,]$CELL_FRAC
	}else{
		v2<-c(0)
	}
	##
	cvaf[EID]<-max(v1)
	mvaf[EID]<-max(v2)
	mcvaf[EID]<-max(c(v1,v2))
}


## Adding VAf
tukb$CHIP_VAF<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(cvaf[cvaf > 0 & cvaf < 0.1]),]$CHIP_VAF<-1
tukb[tukb$eid %in% names(cvaf[cvaf >= 0.1]),]$CHIP_VAF<-2
tukb$CHIP_VAF<-as.factor(tukb$CHIP_VAF)

## Mca
tukb$mCA_VAF<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(mvaf[mvaf > 0 & mvaf < 0.1]),]$mCA_VAF<-1
tukb[tukb$eid %in% names(mvaf[mvaf >= 0.1]),]$mCA_VAF<-2
tukb$mCA_VAF<-as.factor(tukb$mCA_VAF)

## Combined
tukb$CH_VAF<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(mcvaf[mcvaf > 0 & mcvaf < 0.1]),]$CH_VAF<-1
tukb[tukb$eid %in% names(mcvaf[mcvaf >= 0.1]),]$CH_VAF<-2
tukb$CH_VAF<-as.factor(tukb$CH_VAF)

## Number of variants
tukb$CHIP_N<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(Tchip[Tchip == 1]),]$CHIP_N<-1
tukb[tukb$eid %in% names(Tchip[Tchip > 1]),]$CHIP_N<-2
tukb$CHIP_N<-as.factor(tukb$CHIP_N)

## mCA
tukb$mCA_N<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(Tcnv[Tcnv == 1]),]$mCA_N<-1
tukb[tukb$eid %in% names(Tcnv[Tcnv > 1]),]$mCA_N<-2
tukb$mCA_N<-as.factor(tukb$mCA_N)

## Combined
tukb$CH_N<-rep(0,nrow(tukb))
tukb[tukb$eid %in% names(Tch[Tch == 1]),]$CH_N<-1
tukb[tukb$eid %in% names(Tch[Tch > 1]),]$CH_N<-2
tukb$CH_N<-as.factor(tukb$CH_N)


## CHIP and VAF features
tukb$CHIP_var<-rep(0,nrow(tukb))
tukb[tukb$CHIP_N == 1 & tukb$CHIP_VAF == 1,]$CHIP_var<-1
tukb[tukb$CHIP_N == 1 & tukb$CHIP_VAF == 2,]$CHIP_var<-2
tukb[tukb$CHIP_N == 2,]$CHIP_var<-3
tukb$CHIP_var<-as.factor(tukb$CHIP_var)

## mCA
tukb$mCA_var<-rep(0,nrow(tukb))
tukb[tukb$mCA_N == 1 & tukb$mCA_VAF == 1,]$mCA_var<-1
tukb[tukb$mCA_N == 1 & tukb$mCA_VAF == 2,]$mCA_var<-2
tukb[tukb$mCA_N == 2,]$mCA_var<-3
tukb$mCA_var<-as.factor(tukb$mCA_var)

## Combined
tukb$CH_var<-rep(0,nrow(tukb))
tukb[tukb$CH_N == 1 & tukb$CH_VAF == 1,]$CH_var<-1
tukb[tukb$CH_N == 1 & tukb$CH_VAF == 2,]$CH_var<-2
tukb[tukb$CH_N == 2,]$CH_var<-3
tukb$CH_var<-as.factor(tukb$CH_var)


## MCI
tukb$MCI<-rep("Normal",nrow(tukb))
tukb[tukb$BCI == "ALC_high",]$MCI<-"High"
tukb$MCI<-factor(tukb$MCI,levels=c("Normal","High"))


#=====================================================================
## Ten fold cross validation


## Positive and negative data
t1<-unique(tukb[tukb$CLL == 1,])
t2<-unique(tukb[tukb$CLL == 0 & tukb$Myeloid == 0 & tukb$Lymphoid == 0,])


## sample
PT1<-PT2<-PT3<-PT4<-PT5<-PT6<-PT7<-PT8<-PT9<-PT10<-Lab<-c()


S11<-sample(1:nrow(t1),nrow(t1),replace=FALSE)
S12<-sample(1:nrow(t2),nrow(t2),replace=FALSE)
N11<-round(nrow(t1)/10)
N12<-round(nrow(t2)/10)


for (j in 1:10){
	## Subset data
	ST11<-(j-1)*N11+1
	ST12<-(j-1)*N12+1
	EN11<-j*N11
	EN12<-j*N12
	if (j == 10){
		EN11<-nrow(t1)
		EN12<-nrow(t2)
	}
	## Pull training and testing data
	## Test data
	t11<-t1[S11[ST11:EN11],]
	t12<-t2[S12[ST12:EN12],]
	temp1<-rbind(t11,t12)
	Lab<-c(Lab,temp1$CLL)
	## Training data
	temp2<-rbind(t1,t2)
	temp2<-temp2[!temp2$eid %in% temp1$eid,]
	#================================================================
	# Train and predict
	## Demographic
	mod1<-glm(CLL ~ Age_categ + Sex + EverSmoke + Caucasian, data =  temp2,family="binomial")
	PT1<-c(PT1,predict(mod1,temp1,type="response"))
	# Add CHIP
	mod2<-glm(CLL ~ Age_categ + Sex + EverSmoke + Caucasian + CHIP_N + CHIP_VAF, data =  temp2,family="binomial")
	PT2<-c(PT2,predict(mod2,temp1,type="response"))
	## Add mCA only
	mod3<-glm(CLL ~ Age_categ + Sex + EverSmoke + Caucasian + mCA_N + mCA_VAF, data = temp2,family="binomial")
	PT3<-c(PT3,predict(mod3,temp1,type="response"))
	## Add both as separate variables
	mod4<-glm(CLL ~ Age_categ + Sex + EverSmoke + Caucasian + CHIP_N + CHIP_VAF + mCA_N + mCA_VAF, data = temp2,family="binomial")
	PT4<-c(PT4,predict(mod4,temp1,type="response"))
	## Add both as combined variable
	#mod5<-glm(CLL ~ Age_categ + Sex + EverSmoke + Caucasian + CH_var, data = temp2,family="binomial")
	#PT5<-c(PT5,predict(mod5,temp1,type="response"))
	## Add CBC	
	mod6<-glm(CLL ~ Age_categ + Sex + EverSmoke + Caucasian + lymphocyte_count + MCI, data = temp2,family="binomial")
	PT6<-c(PT6,predict(mod6,temp1,type="response"))
	## Adding CBC + CHIP
	mod7<-glm(CLL ~ Age_categ + Sex + EverSmoke + Caucasian + lymphocyte_count + MCI + CHIP_N + CHIP_VAF, data = temp2,family="binomial")
	PT7<-c(PT7,predict(mod7,temp1,type="response"))
	## Adding CBC + mCA
	mod8<-glm(CLL ~ Age_categ + Sex + EverSmoke + Caucasian + lymphocyte_count + MCI + mCA_N + mCA_VAF, data = temp2,family="binomial")
	PT8<-c(PT8,predict(mod8,temp1,type="response"))
	## Adding CBC + CHIP + mCA
	mod9<-glm(CLL ~ Age_categ + Sex + EverSmoke + Caucasian + lymphocyte_count + MCI + CHIP_N + CHIP_VAF + mCA_N + mCA_VAF, data = temp2,family="binomial")
	PT9<-c(PT9,predict(mod9,temp1,type="response"))
}


##
P1<-roc(Lab,PT1,plot=TRUE)
P2<-roc(Lab,PT2,plot=TRUE)
P3<-roc(Lab,PT3,plot=TRUE)
P4<-roc(Lab,PT4,plot=TRUE)
#P5<-roc(Lab,PT5,plot=TRUE)
P6<-roc(Lab,PT6,plot=TRUE)
P7<-roc(Lab,PT7,plot=TRUE)
P8<-roc(Lab,PT8,plot=TRUE)
P9<-roc(Lab,PT9,plot=TRUE)


## plot the AUC
pdf(paste(figDir,"scripts/Final_figs/Fig_S6b_LY.pdf",sep=""),width=9,height=7)
r1 <- plot(P1, print.auc = TRUE, lty = 1, col = "grey", print.auc.y = .05, add = FALSE)
r2 <- plot(P2, print.auc = TRUE, lty = 2, col = "blue", print.auc.y = .1, add = TRUE)
r3 <- plot(P3, print.auc = TRUE, lty = 2, col = "green", print.auc.y = .15, add = TRUE)
r4 <- plot(P4, print.auc = TRUE, lty = 2, col = "orange", print.auc.y = .2, add = TRUE)
#r5 <- plot(P5, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .25, add = TRUE)
r6 <- plot(P6, print.auc = TRUE, lty = 1, col = "black", print.auc.y = .25, add = TRUE)
r7 <- plot(P7, print.auc = TRUE, lty = 1, col = "blue", print.auc.y = .3, add = TRUE)
r8 <- plot(P8, print.auc = TRUE, lty = 1, col = "green", print.auc.y = .35, add = TRUE)
r9 <- plot(P9, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .4, add = TRUE)
dev.off()


setEPS()
postscript(paste(figDir,"scripts/Final_figs/Fig_S6b_LY.eps",sep=""),width=9,height=7)
r1 <- plot(P1, print.auc = TRUE, lty = 1, col = "grey", print.auc.y = .05, add = FALSE)
r2 <- plot(P2, print.auc = TRUE, lty = 2, col = "blue", print.auc.y = .1, add = TRUE)
r3 <- plot(P3, print.auc = TRUE, lty = 2, col = "green", print.auc.y = .15, add = TRUE)
r4 <- plot(P4, print.auc = TRUE, lty = 2, col = "orange", print.auc.y = .2, add = TRUE)
#r5 <- plot(P5, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .25, add = TRUE)
r6 <- plot(P6, print.auc = TRUE, lty = 1, col = "black", print.auc.y = .25, add = TRUE)
r7 <- plot(P7, print.auc = TRUE, lty = 1, col = "blue", print.auc.y = .3, add = TRUE)
r8 <- plot(P8, print.auc = TRUE, lty = 1, col = "green", print.auc.y = .35, add = TRUE)
r9 <- plot(P9, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .4, add = TRUE)
dev.off()


pdf(paste(figDir,"scripts/Final_figs/Fig_3d_LY.pdf",sep=""),width=9,height=7)
r1 <- plot(P1, print.auc = TRUE, lty = 1, col = "grey", print.auc.y = .05, add = FALSE)
r6 <- plot(P6, print.auc = TRUE, lty = 1, col = "black", print.auc.y = .25, add = TRUE)
r8 <- plot(P8, print.auc = TRUE, lty = 1, col = "green", print.auc.y = .35, add = TRUE)
r9 <- plot(P9, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .4, add = TRUE)
dev.off()


setEPS()
postscript(paste(figDir,"scripts/Final_figs/Fig_3d_LY.eps",sep=""),width=9,height=7)
r1 <- plot(P1, print.auc = TRUE, lty = 1, col = "grey", print.auc.y = .05, add = FALSE)
r6 <- plot(P6, print.auc = TRUE, lty = 1, col = "black", print.auc.y = .3, add = TRUE)
r8 <- plot(P8, print.auc = TRUE, lty = 1, col = "green", print.auc.y = .35, add = TRUE)
r9 <- plot(P9, print.auc = TRUE, lty = 1, col = "orange", print.auc.y = .45, add = TRUE)
dev.off()

myobj<-list(Lab,PT1,PT2,PT3,PT4,PT5,PT6,PT7,PT8,PT9,PT10)
save(myobj,file=paste(figDir,"scripts/Final_figs/Fig_3c_LY_object.RData",sep=""))

## DF for AUC Fig.
rdf<-data.frame(Label = Lab,
	DEM = PT1,
	DEM_L_CHIP = PT2,
	DEM_L_mCA = PT3,
	DEM_L_CHIP_L_mCA = PT4,
	DEM_CBC = PT6,
	DEM_CBC_L_CHIP = PT7,
	DEM_CBC_L_mCA = PT8,
	DEM_CBC_L_CHIP_L_mCA = PT9)

write.table(rdf,file=paste(figDir,"scripts/Final_figs/figure_df/Fig_3d_LY.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

print (warnings())
}
