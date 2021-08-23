
## Check if del and gain mCAs increase risk of malignancies independet of LOH

##===================================================================
##===================================================================


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
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"



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


##===================================================================
## CNV
cnv<-data.frame(fread(paste(figDir,"data/UKB_500K_CNV.txt",sep="")))
cnv<-cnv[cnv$CHR %in% 1:22,]


LYG<-c("ATM","CIITA","ITPKB/NTRK1","KMT2C","MIR16-1","NCOR2","NOTCH1")
MYG<-c("CBL","CTCF","EP300","JAK2","TP53","MPL/GNB1")

## canonical myeloid/lymphoid
gainLy<-c("gain12q","gain15q","gain17q","gain21q","gain22q13.32","gain2p","gain3q","gain8q","gain9q","tri12","tri18","tri19")
gainMy<-c("gain1q","gain21q","gain9p","tri8")
lossLy<-c("del10p","del10q","del11q","del13q","del14q","del15q","del17p","del1p","del1q","del22q","del6q","del7q","del8p")
lossMy<-c("del12q","del20q","del5q")

##
manMG<-c("TP53","CTCF","MPL/GNB1","CBL")
manLG<-c("TCL1A","ATM","HEATR3/PLCG2/IRF8")
cnv$LCNV<-rep(0,nrow(cnv))
cnv[cnv$lyGenes %in% c(LYG,manMG) | cnv$myGenes %in% c(LYG,manMG) | cnv$Canonical_CA %in% c(gainLy,lossLy),]$LCNV<-1
cnv$MCNV<-rep(0,nrow(cnv))
cnv[cnv$myGenes %in% c(MYG,manLG) | cnv$lyGenes %in% c(MYG,manLG) | cnv$Canonical_CA %in% c(gainMy,lossMy),]$MCNV<-1

## Lymphoid
cnv$Ldel<-rep(0,nrow(cnv))
cnv$Lgain<-rep(0,nrow(cnv))
cnv$Lloh<-rep(0,nrow(cnv))
cnv[cnv$lyGenes %in% c(LYG,manMG) | cnv$myGenes %in% c(LYG,manMG),]$Lloh<-1
cnv[cnv$Canonical_CA %in% gainLy,]$Lgain<-1
cnv[cnv$Canonical_CA %in% lossLy,]$Ldel<-1

## Myeloid
cnv$Mdel<-rep(0,nrow(cnv))
cnv$Mgain<-rep(0,nrow(cnv))
cnv$Mloh<-rep(0,nrow(cnv))
cnv[cnv$myGenes %in% c(MYG,manLG) | cnv$lyGenes %in% c(MYG,manLG),]$Mloh<-1
cnv[cnv$Canonical_CA %in% gainMy,]$Mgain<-1
cnv[cnv$Canonical_CA %in% lossMy,]$Mdel<-1


#=====================================================================
## Analyzing the effect of three groups of CH
## mCA

ukbc<-data.frame(fread(paste(figDir,"Fig_1_UKB500_malignancy_df.txt",sep="")))


ukbc$Caucasian<-as.factor(ukbc$Caucasian)
ukbc$ever_smoked<-factor(ukbc$ever_smoked,levels=c("No","Yes"))
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))



## CH category
ukbc$CH<-rep("Control",nrow(ukbc))
ukbc[ukbc$eid %in% cnv[cnv$LCNV == 1 & cnv$MCNV == 0 & cnv$Lloh == 1,]$eid,]$CH<-"L-mCA, LOH"
ukbc[ukbc$eid %in% cnv[cnv$LCNV == 1 & cnv$MCNV == 0 & cnv$Lloh == 0 & cnv$Ldel == 1,]$eid,]$CH<-"L-mCA, del"
ukbc[ukbc$eid %in% cnv[cnv$LCNV == 1 & cnv$MCNV == 0 & cnv$Lloh == 0 & cnv$Ldel == 0 & cnv$Lgain == 1,]$eid,]$CH<-"L-mCA, gain"
ukbc[ukbc$eid %in% cnv[cnv$LCNV == 0 & cnv$MCNV == 1 & cnv$Mloh == 1,]$eid,]$CH<-"M-mCA, LOH"
ukbc[ukbc$eid %in% cnv[cnv$LCNV == 0 & cnv$MCNV == 1 & cnv$Mloh == 0 & cnv$Mdel == 1,]$eid,]$CH<-"M-mCA, del"
ukbc[ukbc$eid %in% cnv[cnv$LCNV == 0 & cnv$MCNV == 1 & cnv$Mloh == 0 & cnv$Mdel == 0 & cnv$Mgain == 1,]$eid,]$CH<-"M-mCA, gain"

ukbc[ukbc$eid %in% cnv[cnv$LCNV == 1 & cnv$MCNV == 1,]$eid,]$CH<-"ML"
ukbc[ukbc$eid %in% cnv[cnv$LCNV == 0 & cnv$MCNV == 0,]$eid & ukbc$CH == "Control",]$CH<-"Unk"



ukbc<-ukbc[!ukbc$CH %in% c("ML","Unk"),]
ukbc$CH<-factor(ukbc$CH,levels=c("Control","M-mCA, LOH","M-mCA, del","M-mCA, gain","L-mCA, LOH","L-mCA, del","L-mCA, gain"))



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
		legend.title = "mCA", legend=c(0.2,0.75), legend.labs=c("Control","M-mCA, LOH","M-mCA, del","M-mCA, gain","L-mCA, LOH","L-mCA, del","L-mCA, gain"),
		palette=c("black",rep(ctList[2],3),rep(ctList[1],3)),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		linetype = c(1,1,2,3,1,2,3),
		break.time.by = 1, surv.scale="percent")


## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + EverSmoke + Sex + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1)


# ## Tables
mvdf1<-extractCox(cox1,"Myeloid",c("CHM-mCA, LOH","CHM-mCA, del","CHM-mCA, gain"),c("M-mCA, LOH","M-mCA, del","M-mCA, gain"),"Any")
mvdf5<-extractCox(cox1,"Myeloid",c("CHL-mCA, LOH","CHL-mCA, del","CHL-mCA, gain"),c("L-mCA, LOH","L-mCA, del","L-mCA, gain"),"Any")
#mvdf5<-data.frame(Malignancy = "", CHIP = "Lymphoid", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
BL<-data.frame(Malignancy = "", CHIP = "No mCA", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5)

vlist<-c("Control","M-mCA, LOH","M-mCA, del","M-mCA, gain","L-mCA, LOH","L-mCA, del","L-mCA, gain")
tdf$Total<-as.numeric(table(tukb1$CH)[vlist])
tdf$Positive<-as.numeric(table(tukb1[tukb1$Status ==1,]$CH)[vlist])

tdf$CHIP<-factor(tdf$CHIP,levels=c("L-mCA, gain","L-mCA, del","L-mCA, LOH","M-mCA, gain","M-mCA, del","M-mCA, LOH","No mCA"))

# #tdf[tdf == Inf]<-0

p97 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,color="black") +
	geom_pointrange(size=0.8, fatten = 3, shape=15,linetype=c(1,1,2,3,1,2,3)) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		#panel.grid.major.x=element_line(color="gray80",size=0.2),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[1],3),rep(ctList[2],3),"black"))


tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))



P4<-ggarrange(ggarrange(p97,ptdf,ncol=2,nrow=1,widths=c(1,1)),
	p96$plot,
	ncol=2,nrow=1,widths=c(2,1.3))


#=====================================================================
#=====================================================================

## Lymphoid

tukb1<-ukbc
tukb1$has_disease<-tukb1$Lymphoid
tukb1$Censor_date<-tukb1$Lymphoid_date

tukb1<-prepCox(tukb1)

# #####
YLAB<-"Cumulative incidence of\nlymphoid malignancies"
fit<-survfit(Surv(Year,Status) ~ CH,data=tukb1)
p86<-ggsurvplot(fit,data=tukb1,fun="event",
		legend.title = "mCA", legend=c(0.2,0.75), legend.labs=c("Control","M-mCA, LOH","M-mCA, del","M-mCA, gain","L-mCA, LOH","L-mCA, del","L-mCA, gain"),
		palette=c("black",rep(ctList[2],3),rep(ctList[1],3)),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		linetype = c(1,1,2,3,1,2,3),
		break.time.by = 1, surv.scale="percent")


## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb1)


## Tables
mvdf1<-extractCox(cox1,"Myeloid",c("CHM-mCA, LOH","CHM-mCA, del","CHM-mCA, gain"),c("M-mCA, LOH","M-mCA, del","M-mCA, gain"),"Any")
mvdf5<-extractCox(cox1,"Myeloid",c("CHL-mCA, LOH","CHL-mCA, del","CHL-mCA, gain"),c("L-mCA, LOH","L-mCA, del","L-mCA, gain"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "No mCA", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf5)


vlist<-c("Control","M-mCA, LOH","M-mCA, del","M-mCA, gain","L-mCA, LOH","L-mCA, del","L-mCA, gain")
tdf$Total<-as.numeric(table(tukb1$CH)[vlist])
tdf$Positive<-as.numeric(table(tukb1[tukb1$Status ==1,]$CH)[vlist])

tdf$CHIP<-factor(tdf$CHIP,levels=c("L-mCA, gain","L-mCA, del","L-mCA, LOH","M-mCA, gain","M-mCA, del","M-mCA, LOH","No mCA"))


p87 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,color="black") +
	geom_pointrange(size=0.8, fatten = 3, shape=15,linetype=c(1,1,2,3,1,2,3)) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		#panel.grid.major.x=element_line(color="gray80",size=0.2),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep(ctList[1],3),rep(ctList[2],3),"black"))


##

tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))



P5<-ggarrange(ggarrange(p87,ptdf,ncol=2,nrow=1,widths=c(1,1)),
	p86$plot,
	ncol=2,nrow=1,widths=c(2,1.3))


#==========================================================================
## Final figure
## Max A4 size is 8x11

F1<-ggarrange(P4, P5, ncol=1, nrow=2)


pdf(paste(figDir,"scripts/Final_figs/Fig_ED2_mCA_types.pdf",sep=""),width=11,height=6)
print (F1)
dev.off()

