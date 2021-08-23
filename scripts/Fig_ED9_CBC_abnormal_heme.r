#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

CMD<-args[1]


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
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"


forcedCensorDate<-as.Date("2020-03-31")


##===================================================================
##ICD10codes
hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")
hicdlist<-as.character(hicd[hicd$Group %in% c("Myeloid","Lymphoid","Unspecified"),]$ICD)

## ICDlist
mds_icdList<-c(paste("D46",0:9,sep=""),"C946")
aml_icdList<-paste("C92",0,sep="")
oml_icdList<-paste("C92",1:9,sep="")
mpn_icdList<-c("D752","D45","C945","D474","D471","D473")
cll_icdList<-c("C911","C830")
PC_icdList<-c(paste("C90",0:9,sep=""),"D472")
icdList<-as.character(hicd[hicd$Group %in% c("Lymphoid","Plasma_cell"),]$ICD)
nhl<-icdList[!icdList %in% c(PC_icdList,"C911","C830")]
myeloid_icdList<-as.character(hicd[hicd$Group == "Myeloid",]$ICD)


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

##===================================================================
## Functions


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

if (CMD %in% c("WES")){


#=====================================================================
## Reload data

ukbc<-data.frame(fread(paste(figDir,"Fig_1_WES_malignancy_df.txt",sep="")))


ukbc$Caucasian<-as.factor(ukbc$Caucasian)
ukbc$ever_smoked<-factor(ukbc$ever_smoked,levels=c("No","Yes"))
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))


#=====================================================================
## Subsets based on blood indices

ukbc1<-ukbc[(!is.na(ukbc$red_blood_cell_count)) & (!is.na(ukbc$platelet_count)) & (!is.na(ukbc$neutrophil_count)) & (!is.na(ukbc$lymphocyte_count)),]
ukbc12<-ukbc1[ukbc1$red_blood_cell_count > RBC_upper | ukbc1$platelet_count > PLT_upper | ukbc1$neutrophil_count > ANC_upper,]
ukbc11<-ukbc1[ukbc1$red_blood_cell_count < RBC_lower | ukbc1$platelet_count < PLT_lower | ukbc1$neutrophil_count < ANC_lower,]

ukbc22<-ukbc1[ukbc1$lymphocyte_count > ALC_upper,]
ukbc21<-ukbc1[ukbc1$lymphocyte_count < ALC_lower,]
ukbc1<-ukbc1[!ukbc1$eid %in% c(ukbc11$eid,ukbc12$eid,ukbc21$eid,ukbc22$eid),]

Ihigh<-intersect(c(ukbc12$eid),c(ukbc22$eid))
Ilow<-intersect(c(ukbc11$eid),c(ukbc21$eid))
I<-intersect(c(ukbc11$eid),c(ukbc22$eid))

## High lymphocyte and myeloid cells
ukbc12<-ukbc12[!ukbc12$eid %in% c(Ihigh,Ilow),]
ukbc22<-ukbc22[!ukbc22$eid %in% c(Ihigh,Ilow),]
## Low lymphocyte count and low myeloid counts
ukbc11<-ukbc11[!ukbc11$eid %in% c(Ihigh,Ilow),]
ukbc21<-ukbc21[!ukbc21$eid %in% c(Ihigh,Ilow),]
## If ALC is higher and cytopenias, include it in ALC
ukbc11<-ukbc11[!ukbc11$eid %in% I,]
ukbc11<-ukbc11[!ukbc11$eid %in% c(ukbc12$eid),]


## UKBC
ukbc1$BCI<-rep("Normal",nrow(ukbc1))
ukbc11$BCI<-rep("Plt_low",nrow(ukbc11))
ukbc12$BCI<-rep("Plt_high",nrow(ukbc12))
ukbc21$BCI<-rep("ALC_low",nrow(ukbc21))
ukbc22$BCI<-rep("ALC_high",nrow(ukbc22))


ukbc1<-rbind(ukbc1,ukbc11,ukbc12,ukbc21,ukbc22)
ukbc1$BCI<-factor(ukbc1$BCI,levels=c("Normal","Plt_low","Plt_high","ALC_low","ALC_high"))

#=====================================================================
## Myeloid leukemia

tukb<-ukbc1
tukb$has_disease<-tukb$Myeloid
tukb$Censor_date<-tukb$Myeloid_date

tukb<-prepCox(tukb)

#####
YLAB<-"Cumulative incidence of\nmyeloid malignancies"
fit<-survfit(Surv(Year,Status) ~ BCI,data=tukb)
p1<-ggsurvplot(fit,data=tukb,fun="event",
		legend.title = "CHIP", legend=c(0.2,0.75), legend.labs=c("Normal","Plt_low","Plt_high","ALC_low","ALC_high"),
		palette=c("black",rep(ctList[2],2),rep(ctList[1],2)),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		linetype = c(1,1,2,1,2),
		break.time.by = 1, surv.scale="percent")

## Cox
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + BCI, data = tukb)


mvdf1<-extractCox(cox1,"Myeloid",c("BCIPlt_low","BCIPlt_high","BCIALC_low","BCIALC_high"),c("Plt_low","Plt_high","ALC_low","ALC_high"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "Baseline", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1)
tdf$Total<-as.numeric(table(tukb$BCI)[c("Normal","Plt_low","Plt_high","ALC_low","ALC_high")])
tdf$Positive<-as.numeric(table(tukb[tukb$Status ==1,]$BCI)[c("Normal","Plt_low","Plt_high","ALC_low","ALC_high")])

tdf$CHIP<-factor(tdf$CHIP,levels=c("ALC_high","ALC_low","Plt_high","Plt_low","Baseline"))

if (nrow(tdf[tdf$Positive <= 1,]) >= 1){
	tdf[tdf$Positive <= 1,]$HR<-NA
	tdf[tdf$Positive <= 1,]$CI_high<-NA
	tdf[tdf$Positive <= 1,]$CI_low<-NA
	tdf[tdf$Positive <= 1,]$P<-NA
}

p2 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8, fatten = 3, shape=15,linetype=c(1,1,2,1,2)) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_bw()+
	theme(panel.grid=element_blank(),
		panel.grid.major.x=element_line(color="grey"),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[1],2),rep(ctList[2],2),"black"))+
	scale_y_log10()


## Fig
tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[tdf$CHIP != "CH",c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf1<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


P1<-ggarrange(ggarrange(p2,ptdf1,ncol=2,nrow=2,heights=c(2,1),widths=c(1,1.5)),
	p1$plot,
	ncol=2,nrow=1,widths=c(2,1.3))

#=====================================================================
## CLL

tukb<-ukbc1
tukb$has_disease<-tukb$CLL
tukb$Censor_date<-tukb$CLL_date

tukb<-prepCox(tukb)


#####
YLAB<-"Cumulative ncidence of\nCLL/SLL"
fit<-survfit(Surv(Year,Status) ~ BCI,data=tukb)
p11<-ggsurvplot(fit,data=tukb,fun="event",
		legend.title = "CHIP", legend=c(0.2,0.75), legend.labs=c("Normal","Plt_low","Plt_high","ALC_low","ALC_high"),
		palette=c("black",rep(ctList[2],2),rep(ctList[1],2)),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		linetype = c(1,1,2,1,2),
		break.time.by = 1, surv.scale="percent")

## Cox
cox12 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + BCI, data = tukb)


mvdf1<-extractCox(cox12,"Myeloid",c("BCIPlt_low","BCIPlt_high","BCIALC_low","BCIALC_high"),c("Plt_low","Plt_high","ALC_low","ALC_high"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "Baseline", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1)
tdf$Total<-as.numeric(table(tukb$BCI)[c("Normal","Plt_low","Plt_high","ALC_low","ALC_high")])
tdf$Positive<-as.numeric(table(tukb[tukb$Status ==1,]$BCI)[c("Normal","Plt_low","Plt_high","ALC_low","ALC_high")])

tdf$CHIP<-factor(tdf$CHIP,levels=c("ALC_high","ALC_low","Plt_high","Plt_low","Baseline"))

if (nrow(tdf[tdf$Positive <= 1,]) >= 1){
	tdf[tdf$Positive <= 1,]$HR<-NA
	tdf[tdf$Positive <= 1,]$CI_high<-NA
	tdf[tdf$Positive <= 1,]$CI_low<-NA
	tdf[tdf$Positive <= 1,]$P<-NA
}

p12 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8, fatten = 3, shape=15,linetype=c(1,1,2,1,2)) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_bw()+
	theme(panel.grid=element_blank(),
		panel.grid.major.x=element_line(color="grey"),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[1],2),rep(ctList[2],2),"black"))+
	scale_y_log10()


## Fig
tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[tdf$CHIP != "CH",c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf12<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


P2<-ggarrange(ggarrange(p12,ptdf12,ncol=2,nrow=2,heights=c(2,1),widths=c(1,1.5)),
	p11$plot,
	ncol=2,nrow=1,widths=c(2,1.3))

P<-ggarrange(P1,P2,nrow=2)

pdf(paste(figDir,"scripts/Final_figs/Fig_ED9a_CBC_heme_malignancy.pdf",sep=""),width=12,height=4)
print (P)
dev.off()

}

#==============================================================================


if (CMD %in% c("500")){

ukbc<-data.frame(fread(paste(figDir,"Fig_1_UKB500_malignancy_df.txt",sep="")))


ukbc$Caucasian<-as.factor(ukbc$Caucasian)
ukbc$ever_smoked<-factor(ukbc$ever_smoked,levels=c("No","Yes"))
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))


#=====================================================================
## Group data

ukbc1<-ukbc[(!is.na(ukbc$red_blood_cell_count)) & (!is.na(ukbc$platelet_count)) & (!is.na(ukbc$neutrophil_count)) & (!is.na(ukbc$lymphocyte_count)),]
ukbc12<-ukbc1[ukbc1$red_blood_cell_count > RBC_upper | ukbc1$platelet_count > PLT_upper | ukbc1$neutrophil_count > ANC_upper,]
ukbc11<-ukbc1[ukbc1$red_blood_cell_count < RBC_lower | ukbc1$platelet_count < PLT_lower | ukbc1$neutrophil_count < ANC_lower,]

ukbc22<-ukbc1[ukbc1$lymphocyte_count > ALC_upper,]
ukbc21<-ukbc1[ukbc1$lymphocyte_count < ALC_lower,]
ukbc1<-ukbc1[!ukbc1$eid %in% c(ukbc11$eid,ukbc12$eid,ukbc21$eid,ukbc22$eid),]


Ihigh<-intersect(c(ukbc12$eid),c(ukbc22$eid))
Ilow<-intersect(c(ukbc11$eid),c(ukbc21$eid))
I<-intersect(c(ukbc11$eid),c(ukbc22$eid))

## High lymphocyte and myeloid cells
ukbc12<-ukbc12[!ukbc12$eid %in% c(Ihigh,Ilow),]
ukbc22<-ukbc22[!ukbc22$eid %in% c(Ihigh,Ilow),]
## Low lymphocyte count and low myeloid counts
ukbc11<-ukbc11[!ukbc11$eid %in% c(Ihigh,Ilow),]
ukbc21<-ukbc21[!ukbc21$eid %in% c(Ihigh,Ilow),]
## If ALC is higher and cytopenias, include it in ALC
ukbc11<-ukbc11[!ukbc11$eid %in% I,]
ukbc11<-ukbc11[!ukbc11$eid %in% c(ukbc12$eid),]


ukbc1$BCI<-rep("Normal",nrow(ukbc1))
ukbc11$BCI<-rep("Plt_low",nrow(ukbc11))
ukbc12$BCI<-rep("Plt_high",nrow(ukbc12))
ukbc21$BCI<-rep("ALC_low",nrow(ukbc21))
ukbc22$BCI<-rep("ALC_high",nrow(ukbc22))

## Combine data
ukbc1<-rbind(ukbc1,ukbc11,ukbc12,ukbc21,ukbc22)
ukbc1$BCI<-factor(ukbc1$BCI,levels=c("Normal","Plt_low","Plt_high","ALC_low","ALC_high"))



#=====================================================================
## Myeloid leukemia

tukb<-ukbc1
tukb$has_disease<-tukb$Myeloid
tukb$Censor_date<-tukb$Myeloid_date

tukb<-prepCox(tukb)


#####
YLAB<-"Incidence of myeloid leukemia"
fit<-survfit(Surv(Year,Status) ~ BCI,data=tukb)
p1<-ggsurvplot(fit,data=tukb,fun="event",
		legend.title = "CHIP", legend=c(0.2,0.75), legend.labs=c("Normal","Plt_low","Plt_high","ALC_low","ALC_high"),
		palette=c("black",rep(ctList[2],2),rep(ctList[1],2)),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		linetype = c(1,1,2,1,2),
		break.time.by = 1, surv.scale="percent")

## Cox
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + EverSmoke + Sex + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + BCI, data = tukb)


mvdf1<-extractCox(cox1,"Myeloid",c("BCIPlt_low","BCIPlt_high","BCIALC_low","BCIALC_high"),c("Plt_low","Plt_high","ALC_low","ALC_high"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "Baseline", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1)
tdf$Total<-as.numeric(table(tukb$BCI)[c("Normal","Plt_low","Plt_high","ALC_low","ALC_high")])
tdf$Positive<-as.numeric(table(tukb[tukb$Status ==1,]$BCI)[c("Normal","Plt_low","Plt_high","ALC_low","ALC_high")])

tdf$CHIP<-factor(tdf$CHIP,levels=c("ALC_high","ALC_low","Plt_high","Plt_low","Baseline"))

p2 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8, fatten = 3, shape=15,linetype=c(1,1,2,1,2)) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_bw()+
	theme(panel.grid=element_blank(),
		panel.grid.major.x=element_line(color="grey"),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[1],2),rep(ctList[2],2),"black"))+
	scale_y_log10()


## Fig
tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[tdf$CHIP != "CH",c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf1<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


P1<-ggarrange(ggarrange(p2,ptdf1,ncol=2,nrow=2,heights=c(2,1),widths=c(1,1.5)),
	p1$plot,
	ncol=2,nrow=1,widths=c(2,1.3))


#=====================================================================
## CLL

tukb<-ukbc1
tukb$has_disease<-tukb$CLL
tukb$Censor_date<-tukb$CLL_date

tukb<-prepCox(tukb)

#####
YLAB<-"Incidence of CLL"
fit<-survfit(Surv(Year,Status) ~ BCI,data=tukb)
p11<-ggsurvplot(fit,data=tukb,fun="event",
		legend.title = "CHIP", legend=c(0.2,0.75), legend.labs=c("Normal","Plt_low","Plt_high","ALC_low","ALC_high"),
		palette=c("black",rep(ctList[2],2),rep(ctList[1],2)),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		linetype = c(1,1,2,1,2),
		break.time.by = 1, surv.scale="percent")

## Cox
cox12 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + BCI, data = tukb)


mvdf1<-extractCox(cox12,"Myeloid",c("BCIPlt_low","BCIPlt_high","BCIALC_low","BCIALC_high"),c("Plt_low","Plt_high","ALC_low","ALC_high"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "Baseline", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1)
tdf$Total<-as.numeric(table(tukb$BCI)[c("Normal","Plt_low","Plt_high","ALC_low","ALC_high")])
tdf$Positive<-as.numeric(table(tukb[tukb$Status ==1,]$BCI)[c("Normal","Plt_low","Plt_high","ALC_low","ALC_high")])

tdf$CHIP<-factor(tdf$CHIP,levels=c("ALC_high","ALC_low","Plt_high","Plt_low","Baseline"))

p12 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8, fatten = 3, shape=15,linetype=c(1,1,2,1,2)) +
	#facet_grid(LMCHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_bw()+
	theme(panel.grid=element_blank(),
		panel.grid.major.x=element_line(color="grey"),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[1],2),rep(ctList[2],2),"black"))+
	scale_y_log10()


## Fig
tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[tdf$CHIP != "CH",c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf12<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


P2<-ggarrange(ggarrange(p12,ptdf12,ncol=2,nrow=2,heights=c(2,1),widths=c(1,1.5)),
	p11$plot,
	ncol=2,nrow=1,widths=c(2,1.3))

P<-ggarrange(P1,P2,nrow=2)

pdf(paste(figDir,"scripts/Final_figs/Fig_ED9b_CBC_heme_malignancy.pdf",sep=""),width=12,height=4)
print (P)
dev.off()

}
