#!/usr/bin/Rscript


## Fig 3a,b

## libraries

library(ggplot2)
library(data.table)
library(colortools)
library(survival)
library(survminer)
library(ggpubr)

ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
hemeICDFile<-"metta/heme_malignancies_icd10.txt"

figDir<-""


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


#=============================================================================
## WES sample

ukbc<-data.frame(fread(paste(figDir,"Fig_1_WES_malignancy_df.txt",sep="")))

ukbc$Caucasian<-as.factor(ukbc$Caucasian)
ukbc$ever_smoked<-factor(ukbc$ever_smoked,levels=c("No","Yes"))
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))


#=====================================================================
## Subsets based on blood indices

ukbc1<-ukbc[(!is.na(ukbc$red_blood_cell_count)) & (!is.na(ukbc$platelet_count)) & (!is.na(ukbc$neutrophil_count)) & (!is.na(ukbc$lymphocyte_count)),]
ukbc11<-ukbc1[ukbc1$red_blood_cell_count < RBC_lower | ukbc1$platelet_count < PLT_lower | ukbc1$neutrophil_count < ANC_lower,]
ukbc12<-ukbc1[ukbc1$red_blood_cell_count > RBC_upper | ukbc1$platelet_count > PLT_upper | ukbc1$neutrophil_count > ANC_upper,]
ukbc12<-ukbc12[!ukbc12$eid %in% ukbc11$eid,]

ukbc1<-ukbc1[!ukbc1$eid %in% c(ukbc11$eid,ukbc12$eid),]


## UKBC
ukbc1$BCI<-rep("Normal",nrow(ukbc1))
ukbc11$BCI<-rep("Plt_low",nrow(ukbc11))
ukbc12$BCI<-rep("Plt_high",nrow(ukbc12))


ukbc1<-rbind(ukbc1,ukbc11,ukbc12)
ukbc1$BCI<-factor(ukbc1$BCI,levels=c("Normal","Plt_low","Plt_high"))

#=====================================================================
## Myeloid

temp<-ukbc1[ukbc1$BCI %in% c("Normal","Plt_low","Plt_high"),]

temp$has_disease<-temp$Myeloid
temp$Censor_date<-temp$Myeloid_date
temp<-prepCox(temp)

#=====================================================================
## Make groups for VAF, CH, etc

T<-table(c(mychip$eid,cnv[cnv$MCNV == 1,]$eid))

## CHIP
t1<-mychip[mychip$eid %in% names(T[T==1]),]
maxVAF<-t1$VAF
names(maxVAF)<-as.character(t1$eid)

## mCA
t1<-cnv[cnv$MCNV == 1 & cnv$eid %in% names(T[T==1]),]
maxVAF1<-t1$CELL_FRAC
names(maxVAF1)<-as.character(t1$eid)
maxVAF<-c(maxVAF,maxVAF1)

for (EID in names(T[T>1])){
	v1<-mychip[mychip$eid == EID,]$VAF
	v2<-cnv[cnv$eid == EID & cnv$MCNV == 1,]$CELL_FRAC
	maxVAF[EID]<-max(c(v1,v2))
}

temp$maxVAF<-maxVAF[as.character(temp$eid)]
temp[is.na(temp$maxVAF),]$maxVAF<-0

## N_hits
temp$nCH<-as.numeric(T[as.character(temp$eid)])
temp[is.na(temp$nCH),]$nCH<-0


## CH class
temp$CH<-rep("Control",nrow(temp))
temp[!temp$BCI %in% c("Plt_high","Plt_low") & temp$nCH == 1 & temp$maxVAF < 0.1,]$CH<-"N_small"
temp[!temp$BCI %in% c("Plt_high","Plt_low") & temp$nCH == 1 & temp$maxVAF >= 0.1,]$CH<-"N_large"
temp[!temp$BCI %in% c("Plt_high","Plt_low") & temp$nCH > 1,]$CH<-"N_double"
temp[temp$BCI %in% c("Plt_high","Plt_low") & temp$nCH == 0,]$CH<-"I_control"
temp[temp$BCI %in% c("Plt_high","Plt_low") & temp$nCH == 1 & temp$maxVAF < 0.1,]$CH<-"I_small"
temp[temp$BCI %in% c("Plt_high","Plt_low") & temp$nCH == 1 & temp$maxVAF >= 0.1,]$CH<-"I_large"
temp[temp$BCI %in% c("Plt_high","Plt_low") & temp$nCH > 1,]$CH<-"I_double"

temp$CH<-factor(temp$CH,levels=c("Control","N_small","N_large","N_double","I_control","I_small","I_large","I_double"))

print (nrow(temp))

tmyeloid<-temp[,c("eid","CH","Status","Year")]
colnames(tmyeloid)<-c("eid","CH","Status_MY","Year_MY")


## Myeloid malignancies
YLAB<-"Cumulative incidence of\nmyeloid malignancies"
fit<-survfit(Surv(Year,Status) ~ CH,data=temp)
p3<-ggsurvplot(fit,data=temp,fun="event",
		legend.title = "CHIP", legend=c(0.2,0.75), legend.labs=c("No CH","1 CH, VAF<0.1","1 CH, VAF>=0.1",">1 CH","I_Control","I_1 CH, VAF<0.1","I_1 CH, VAF>=0.1","I_>1 CH"),
		palette=c("grey",rep("black",3),"grey",rep(ctList[2],3)),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		linetype = c(1,1,2,3,2,1,2,3),
		break.time.by = 1, surv.scale="percent")

## cox
cox3 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = temp)

## Tables
mvdf2<-extractCox(cox3,"Myeloid",c("CHN_small","CHN_large","CHN_double"),c("CH, VAF<0.1","CH, VAF>=0.1",">1 CH"),"Any")
mvdf3<-extractCox(cox3,"Myeloid",c("CHI_control","CHI_small","CHI_large","CHI_double"),c("CBC- No CH","CBC- CH, VAF<0.1","CBC- CH, VAF>=0.1","CBC- >1 CH"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "No CH", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf2,mvdf3)
N<-c("Control","N_small","N_large","N_double","I_control","I_small","I_large","I_double")
tdf$Total<-as.numeric(table(temp$CH)[N])
tdf$Positive<-as.numeric(table(temp[temp$Status ==1,]$CH)[N])


tdf$CHIP<-factor(tdf$CHIP,levels=c("CBC- >1 CH","CBC- CH, VAF>=0.1","CBC- CH, VAF<0.1","CBC- No CH",">1 CH","CH, VAF>=0.1","CH, VAF<0.1","No CH"))
p4 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8, fatten = 3, shape=15,linetype=c(1,1,2,3,1,1,2,3)) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_bw()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[2],3),"grey",rep("black",3),"grey"))+
	scale_y_log10()


forestdf1<-tdf[,c("CHIP","HR","CI_low","CI_high")]
forestdf1$Malignancy<-rep("Myeloid")


## Fig
tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[tdf$CHIP != "CH",c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf2<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


P4<-ggarrange(ggarrange(p4,ptdf2,ncol=2,nrow=2,heights=c(2,1),widths=c(1,1)),
	p3$plot,
	ncol=2,nrow=1,widths=c(2,1.3))


#=====================================================================
#=====================================================================
#=====================================================================
## CLL

print ("Starting CLL")

#=====================================================================
## Subsets based on blood indices

ukbc1<-ukbc[(!is.na(ukbc$red_blood_cell_count)) & (!is.na(ukbc$platelet_count)) & (!is.na(ukbc$neutrophil_count)) & (!is.na(ukbc$lymphocyte_count)),]

ukbc22<-ukbc1[ukbc1$lymphocyte_count > ALC_upper,]
ukbc1<-ukbc1[!ukbc1$eid %in% c(ukbc22$eid),]

## UKBC
ukbc1$BCI<-rep("Normal",nrow(ukbc1))
ukbc22$BCI<-rep("ALC_high",nrow(ukbc22))


ukbc1<-rbind(ukbc1,ukbc22)
ukbc1$BCI<-factor(ukbc1$BCI,levels=c("Normal","ALC_high"))

#=====================================================================

temp<-ukbc1[ukbc1$BCI %in% c("Normal","ALC_high"),]

temp$has_disease<-temp$CLL
temp$Censor_date<-temp$CLL_date
temp<-prepCox(temp)

#=====================================================================
## Make groups for VAF, CH, etc

T<-table(c(chip$eid,cnv[cnv$LCNV == 1,]$eid))

## CHIP
t1<-chip[chip$eid %in% names(T[T==1]),]
maxVAF<-t1$VAF
names(maxVAF)<-as.character(t1$eid)

## mCA
t1<-cnv[cnv$LCNV == 1 & cnv$eid %in% names(T[T==1]),]
maxVAF1<-t1$CELL_FRAC
names(maxVAF1)<-as.character(t1$eid)
maxVAF<-c(maxVAF,maxVAF1)

for (EID in names(T[T>1])){
	v1<-chip[chip$eid == EID,]$VAF
	v2<-cnv[cnv$eid == EID & cnv$LCNV == 1,]$CELL_FRAC
	maxVAF[EID]<-max(c(v1,v2))
}

temp$maxVAF<-maxVAF[as.character(temp$eid)]
temp[is.na(temp$maxVAF),]$maxVAF<-0

## N_hits
temp$nCH<-as.numeric(T[as.character(temp$eid)])
temp[is.na(temp$nCH),]$nCH<-0


## CH class
temp$CH<-rep("Control",nrow(temp))
temp[temp$BCI != "ALC_high" & temp$nCH == 1 & temp$maxVAF < 0.1,]$CH<-"N_small"
temp[temp$BCI != "ALC_high" & temp$nCH == 1 & temp$maxVAF >= 0.1,]$CH<-"N_large"
temp[temp$BCI != "ALC_high" & temp$nCH > 1,]$CH<-"N_double"
temp[temp$BCI == "ALC_high" & temp$nCH == 0,]$CH<-"I_control"
temp[temp$BCI == "ALC_high" & temp$nCH == 1 & temp$maxVAF < 0.1,]$CH<-"I_small"
temp[temp$BCI == "ALC_high" & temp$nCH == 1 & temp$maxVAF >= 0.1,]$CH<-"I_large"
temp[temp$BCI == "ALC_high" & temp$nCH > 1,]$CH<-"I_double"

temp$CH<-factor(temp$CH,levels=c("Control","N_small","N_large","N_double","I_control","I_small","I_large","I_double"))

print (nrow(temp))

ST<-temp$Status
YR<-temp$Year
ALC<-temp$CH
names(ST)<-names(YR)<-names(ALC)<-as.character(temp$eid)

tmyeloid$Status_LY<-ST[as.character(tmyeloid$eid)]
tmyeloid$Year_LY<-YR[as.character(tmyeloid$eid)]
tmyeloid$CH_LY<-ALC[as.character(tmyeloid$eid)]

## reformat data
print (table(tmyeloid$CH))
tmyeloid$M_CH<-rep("No M-CH",rep(nrow(tmyeloid)))
tmyeloid[tmyeloid$CH %in% c("N_small","I_small"),]$M_CH<-"1 M-CH, VAF<0.1"
tmyeloid[tmyeloid$CH %in% c("N_large","I_large"),]$M_CH<-"1 M-CH, VAF>=0.1"
tmyeloid[tmyeloid$CH %in% c("N_double","I_double"),]$M_CH<-">1 M-CH"

## Myeloid cell parameters
tmyeloid$CBC_myeloid<-rep("Normal",rep(nrow(tmyeloid)))
tmyeloid[tmyeloid$CH %in% c("I_control","I_small","I_large","I_double"),]$CBC_myeloid<-"Abnormal"

## reformat data
tmyeloid$L_CH<-rep("No L-CH",rep(nrow(tmyeloid)))
tmyeloid[tmyeloid$CH_LY %in% c("N_small","I_small"),]$L_CH<-"1 L-CH, VAF<0.1"
tmyeloid[tmyeloid$CH_LY %in% c("N_large","I_large"),]$L_CH<-"1 L-CH, VAF>=0.1"
tmyeloid[tmyeloid$CH_LY %in% c("N_double","I_double"),]$L_CH<-">1 L-CH"

## Myeloid cell parameters
tmyeloid$Lymphocyte_count<-rep("Normal",rep(nrow(tmyeloid)))
tmyeloid[tmyeloid$CH_LY %in% c("I_control","I_small","I_large","I_double"),]$Lymphocyte_count<-"High"


t<-tmyeloid[,c("M_CH","CBC_myeloid","Status_MY","Year_MY","L_CH","Lymphocyte_count","Status_LY","Year_LY")]
write.table(t,file = paste(figDir,"scripts/Final_figs/figure_df/F3_ab_incidence.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

## Myeloid malignancies
YLAB<-"Cumulative incidence of\nCLL/SLL"
fit<-survfit(Surv(Year,Status) ~ CH,data=temp)
p5<-ggsurvplot(fit,data=temp,fun="event",
		legend.title = "CHIP", legend=c(0.2,0.75), legend.labs=c("No CH","1 CH, VAF<0.1","1 CH, VAF>=0.1",">1 CH","I_Control","I_1 CH, VAF<0.1","I_1 CH, VAF>=0.1","I_>1 CH"),
		palette=c("grey",rep("black",3),"grey",rep(ctList[1],3)),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		linetype = c(1,1,2,3,1,1,2,3),
		break.time.by = 1, surv.scale="percent")

## cox
cox6 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = temp[!temp$CH %in% c("N_large","N_double"),])

## Tables
mvdf1<-extractCox(cox6,"Lymphoid",c("CHN_small"),c("CH, VAF<0.1"),"Any")
mvdf2<-data.frame(Malignancy = "", CHIP = c("CH, VAF>=0.1",">1 CH"), HR = NA, CI_low = NA, CI_high = NA, P = NA, VAF = "Any")
mvdf3<-extractCox(cox6,"Lymphoid",c("CHI_control","CHI_small","CHI_large","CHI_double"),c("CBC- No CH","CBC- CH, VAF<0.1","CBC- CH, VAF>=0.1","CBC- >1 CH"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "No CH", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf2,mvdf3)
N<-c("Control","N_small","N_large","N_double","I_control","I_small","I_large","I_double")
tdf$Total<-as.numeric(table(temp$CH)[N])
tdf$Positive<-as.numeric(table(temp[temp$Status ==1,]$CH)[N])


tdf$CHIP<-factor(tdf$CHIP,levels=c("CBC- >1 CH","CBC- CH, VAF>=0.1","CBC- CH, VAF<0.1","CBC- No CH",">1 CH","CH, VAF>=0.1","CH, VAF<0.1","No CH"))
p6 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8, fatten = 3, shape=15,linetype=c(1,1,2,3,1,1,2,3)) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_bw()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[1],3),"grey",rep("black",3),"grey"))+
	scale_y_log10()


forestdf2<-tdf[,c("CHIP","HR","CI_low","CI_high")]
forestdf2$Malignancy<-rep("Lymphoid")

## Fig
tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[tdf$CHIP != "CH",c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf6<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


P5<-ggarrange(ggarrange(p6,ptdf6,ncol=2,nrow=2,heights=c(2,1),widths=c(1,1)),
	p5$plot,
	ncol=2,nrow=1,widths=c(2,1.3))


forestdf1$CHIP<-gsub("CH","M-CH",forestdf1$CHIP)
forestdf1$CHIP<-gsub("CBC","Abnormal CBC",forestdf1$CHIP)
forestdf2$CHIP<-gsub("CH","L-CH",forestdf2$CHIP)
forestdf2$CHIP<-gsub("CBC","High lymphocyte count",forestdf2$CHIP)


t<-rbind(forestdf1,forestdf2)

write.table(t,file = paste(figDir,"scripts/Final_figs/figure_df/F3_ab_forestPlots.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

#=======================================================================================
#========================================================================================

P<-ggarrange(P4,P5,nrow=2)

pdf(paste(figDir,"scripts/Final_figs/Fig_3ab_CHIP_mCA_CBC.pdf",sep=""),width=10,height=6)
print (P)
dev.off()


