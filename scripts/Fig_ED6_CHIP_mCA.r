## Cases with multiple groups of events in WES cohort


#args = commandArgs(trailingOnly=TRUE)
#CMD<-args[1]


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


#=====================================================================
## Reload data

ukbc<-data.frame(fread(paste(figDir,"Fig_1_WES_malignancy_df.txt",sep="")))


ukbc$Caucasian<-as.factor(ukbc$Caucasian)
ukbc$ever_smoked<-factor(ukbc$ever_smoked,levels=c("No","Yes"))
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))


#=====================================================================
## Analyzing the effect of three groups of CH

tukb<-ukbc
tukb$has_disease<-tukb$Myeloid
tukb$Censor_date<-tukb$Myeloid_date

tukb<-prepCox(tukb)

## CH category
both<-intersect(unique(mychip$eid),unique(cnv[cnv$MCNV == 1,]$eid))

tukb$CH<-rep("Control",nrow(tukb))
tukb[tukb$eid %in% mychip$eid,]$CH<-"CHIP"
tukb[tukb$eid %in% cnv[cnv$MCNV == 1,]$eid,]$CH<-"mCA"
tukb[tukb$eid %in% both,]$CH <- "CHIP_mCA"

tukb$CH<-factor(tukb$CH,levels=c("Control","CHIP","mCA","CHIP_mCA"))


## Incidence
YLAB<-"Cumulative incidence of\nmyeloid malignancies"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb)
p98<-ggsurvplot(fit,data=tukb,fun="event",
		legend.title = "CH", legend=c(0.2,0.75), legend.labs=c("No CH","M-CHIP","M-mCA","M-CHIP + M-mCA"),
		palette=c("black",ctList[2],ctList[3],"orange"),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")


## Cox
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb)
# ## Tables
mvdf1<-extractCox(cox1,"Myeloid","CHCHIP","M-CHIP","Any")
mvdf5<-extractCox(cox1,"Myeloid","CHmCA","M-mCA","Any")
mvdf6<-extractCox(cox1,"Myeloid","CHCHIP_mCA","M-CHIP + M-mCA","Any")
BL<-data.frame(Malignancy = "", CHIP = "No CH", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5,mvdf6)
tdf$Total<-as.numeric(table(tukb$CH)[c("Control","CHIP","mCA","CHIP_mCA")])
tdf$Positive<-as.numeric(table(tukb[tukb$Status ==1,]$CH)[c("Control","CHIP","mCA","CHIP_mCA")])

tdf$CHIP<-factor(tdf$CHIP,levels=c("M-CHIP + M-mCA","M-mCA","M-CHIP","No CH"))


p99 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c("orange",ctList[3],ctList[2],"black"))+
	scale_y_log10()


tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


P2<-ggarrange(ggarrange(p99,ptdf,ncol=2,nrow=2,widths=c(1,1.5),heights=c(1.5,1)),
	p98$plot,
	ncol=2,nrow=2,widths=c(2,1.4),heights=c(1,0.1))


#===================================================================
#===================================================================
#===================================================================
## Lymphoid malignancy


tukb<-ukbc
tukb$has_disease<-tukb$CLL
tukb$Censor_date<-tukb$CLL_date

tukb<-prepCox(tukb)

## CH category
both<-intersect(unique(chip$eid),unique(cnv[cnv$LCNV == 1,]$eid))

tukb$CH<-rep("Control",nrow(tukb))
tukb[tukb$eid %in% chip$eid,]$CH<-"CHIP"
tukb[tukb$eid %in% cnv[cnv$LCNV == 1,]$eid,]$CH<-"mCA"
tukb[tukb$eid %in% both,]$CH <- "CHIP_mCA"

tukb$CH<-factor(tukb$CH,levels=c("Control","CHIP","mCA","CHIP_mCA"))


## Incidence
YLAB<-"Cumulative incidence of\nCLL/SLL"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb)
p88<-ggsurvplot(fit,data=tukb,fun="event",
		legend.title = "CH", legend=c(0.2,0.75), legend.labs=c("No CH","L-CHIP","L-mCA","L-CHIP + L-mCA"),
		palette=c("black",ctList[1],ctList[3],"orange"),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")


## Cox
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb)
#print (getX(cox1,"My"))
# ## Tables
mvdf1<-extractCox(cox1,"Lymphoid","CHCHIP","L-CHIP","Any")
mvdf5<-extractCox(cox1,"Lymphoid","CHmCA","L-mCA","Any")
mvdf6<-extractCox(cox1,"Lymphoid","CHCHIP_mCA","L-CHIP + L-mCA","Any")
#mvdf5<-data.frame(Malignancy = "", CHIP = "Lymphoid", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
BL<-data.frame(Malignancy = "", CHIP = "No CH", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5,mvdf6)
tdf$Total<-as.numeric(table(tukb$CH)[c("Control","CHIP","mCA","CHIP_mCA")])
tdf$Positive<-as.numeric(table(tukb[tukb$Status ==1,]$CH)[c("Control","CHIP","mCA","CHIP_mCA")])

tdf$CHIP<-factor(tdf$CHIP,levels=c("L-CHIP + L-mCA","L-mCA","L-CHIP","No CH"))


p89 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c("orange",ctList[3],ctList[1],"black"))+
	scale_y_log10()


tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf1<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


P3<-ggarrange(ggarrange(p89,ptdf1,ncol=2,nrow=2,widths=c(1,1.5),heights=c(1.5,1)),
	p88$plot,
	ncol=2,nrow=2,widths=c(2,1.4),heights=c(1,0.1))


P<-ggarrange(P2,P3,nrow=2)

pdf(paste(figDir,"scripts/Final_figs/Fig_ED6_CHIP_mCA.pdf",sep=""),width=10,height=6)
print (P)
dev.off()



