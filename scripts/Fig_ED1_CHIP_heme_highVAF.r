


library(ggplot2)
library(data.table)
library(ggpubr)
library(mgcv)
library(colortools)

library(survival)
library(survminer)
library(gridExtra)

ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
figDir<-""



forcedCensorDate<-as.Date("2020-03-31")


#=====================================================================
## Reload data


chip<-read.table(paste(figDir,"1_CHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
mychip<-read.table(paste(figDir,"1_MCHIP_UKB.txt",sep=""),head=TRUE,sep="\t")


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
#=====================================================================
#=====================================================================


ukbc<-data.frame(fread(paste(figDir,"Fig_1_WES_malignancy_df.txt",sep="")))


ukbc$Caucasian<-as.factor(ukbc$Caucasian)
ukbc$ever_smoked<-factor(ukbc$ever_smoked,levels=c("No","Yes"))
ukbc$EverSmoke<-factor(ukbc$ever_smoked,levels=c("No","Yes"))

ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))


#=====================================================================
## Analyzing the effect of three groups of CH

tukb<-ukbc
tukb$has_disease<-tukb$Myeloid
tukb$Censor_date<-tukb$Myeloid_date

tukb<-prepCox(tukb)

tukb$CHIP<-as.factor(tukb$CHIP)
tukb$MCHIP<-as.factor(tukb$MCHIP)

#=============================
## First analyze coxph for the three CH groups separately
## Lymphoid
tp1<-tukb[tukb$MCHIP == 0,]
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CHIP, data = tp1)
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CHIP, data = tp1[tp1$CHIP_size != 1,])
mvdf1<-extractCox(cox1,"Myeloid","CHIP1","Lymphoid","CH-LD")
#mvdf2<-extractCox(cox2,"Myeloid","CHIP_size1","Lymphoid","CH-LD, VAF<0.1")
mvdf3<-extractCox(cox2,"Myeloid","CHIP1","Lymphoid","CH-LD, VAF>=0.1")
mvdf0<-data.frame(Malignancy = "Myeloid", CHIP = "Lymphoid",HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Baseline")
temp1<-rbind(mvdf1,mvdf3)
#temp1$at_risk<-c(nrow(tp1[tp1$CHIP == 1,]),nrow(tp1[tp1$CHIP_size == 1,]),nrow(tp1[tp1$CHIP_size == 2,]))
#temp1$event<-c(nrow(tp1[tp1$CHIP == 1 & tp1$Status ==1,]),nrow(tp1[tp1$CHIP_size == 1 & tp1$Status == 1,]),nrow(tp1[tp1$CHIP_size == 2 & tp1$Status == 1,]))
temp1$at_risk<-c(nrow(tp1[tp1$CHIP == 1,]),nrow(tp1[tp1$CHIP_size == 2,]))
temp1$event<-c(nrow(tp1[tp1$CHIP == 1 & tp1$Status ==1,]),nrow(tp1[tp1$CHIP_size == 2 & tp1$Status == 1,]))
temp1$Col<-rep("Ly",nrow(temp1))

print ("LCHIP done")
###
ly <- temp1

## Myeloid
tp1<-tukb[tukb$CHIP == 0,]
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + MCHIP, data = tp1)
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + MCHIP, data = tp1[tp1$Msize != 1,])


mvdf1<-extractCox(cox1,"Myeloid","MCHIP1","Myeloid","CH-MD")
#mvdf2<-extractCox(cox2,"Myeloid","Msize1","Myeloid","CH-MD, VAF<0.1")
mvdf3<-extractCox(cox2,"Myeloid","MCHIP1","Myeloid","CH-MD, VAF>=0.1")
mvdf0<-data.frame(Malignancy = "Myeloid", CHIP = "Myeloid",HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Baseline")
temp2<-rbind(mvdf0,mvdf1,mvdf3)
temp2$at_risk<-c(nrow(tp1[tp1$MCHIP == 0,]),nrow(tp1[tp1$MCHIP == 1,]),nrow(tp1[tp1$Msize == 2,]))
temp2$event<-c(nrow(tp1[tp1$MCHIP == 0 & tp1$Status == 1,]),nrow(tp1[tp1$MCHIP == 1 & tp1$Status ==1,]),nrow(tp1[tp1$Msize == 2 & tp1$Status == 1,]))
temp2$Col<-rep("My",nrow(temp2))
temp2[temp2$VAF == "Baseline",]$Col<-"Baseline"

temp2<-rbind(temp2,ly)
temp2$Col<-as.factor(temp2$Col)

temp2$VAF<-as.character(temp2$VAF)
temp2$VAF<-gsub("CH-LD","L-CHIP",temp2$VAF)
temp2$VAF<-gsub("CH-MD","M-CHIP",temp2$VAF)
temp2$VAF<-gsub("Baseline","No CHIP",temp2$VAF)
temp2$VAF<-gsub("Any ","",temp2$VAF)

## if n event is <= 1
temp2[temp2$event <= 1,]$HR<-NA
temp2[temp2$event <= 1,]$CI_high<-NA
temp2[temp2$event <= 1,]$CI_low<-NA
temp2[temp2$event <= 1,]$P<-NA

## Table
tdf<-temp2
tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("VAF","at_risk","event","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf2<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
		padding = unit(c(2, 2), "mm")))

temp2$VAF<-factor(temp2$VAF,levels=c("L-CHIP, VAF>=0.1","L-CHIP","M-CHIP, VAF>=0.1","M-CHIP", "No CHIP"))
temp2$Col<-factor(temp2$Col,levels=c("Baseline","My","Ly"))

p2 <- ggplot(data=temp2, aes(x=VAF, y=HR, ymin=CI_low, ymax=CI_high,color=Col)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8,shape=15, fatten=3) +
	#facet_grid(CHIP ~.)+
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c("black",ctList[2],ctList[1]))

F3a<-ggarrange(p2,ptdf2,ncol=2,nrow=1,widths=c(1,1.5))


#=====================================================================
#=====================================================================
## Lymphoid malignancies


tukb<-ukbc
tukb$has_disease<-tukb$Lymphoid
tukb$Censor_date<-tukb$Lymphoid_date

tukb<-prepCox(tukb)

tukb$CHIP<-as.factor(tukb$CHIP)
tukb$MCHIP<-as.factor(tukb$MCHIP)


#=============================
## First analyze coxph for the three CH groups separately
## Lymphoid
tp1<-tukb[tukb$MCHIP == 0,]
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CHIP, data = tp1)
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CHIP, data = tp1[tp1$CHIP_size != 1,])

mvdf1<-extractCox(cox1,"Myeloid","CHIP1","Lymphoid","Any CH-LD")
#mvdf2<-extractCox(cox2,"Myeloid","CHIP_size1","Lymphoid","CH-LD, VAF<0.1")
mvdf3<-extractCox(cox2,"Myeloid","CHIP1","Lymphoid","CH-LD, VAF>=0.1")
mvdf0<-data.frame(Malignancy = "Myeloid", CHIP = "Lymphoid",HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Baseline")
# temp1<-rbind(mvdf0,mvdf1,mvdf2,mvdf3)
temp1<-rbind(mvdf1,mvdf3)
temp1$at_risk<-c(nrow(tp1[tp1$CHIP == 1,]),nrow(tp1[tp1$CHIP_size == 2,]))
temp1$event<-c(nrow(tp1[tp1$CHIP == 1 & tp1$Status ==1,]),nrow(tp1[tp1$CHIP_size == 2 & tp1$Status == 1,]))
temp1$Col<-rep("Ly",nrow(temp1))



#####
ly<-temp1

## Myeloid
tp1<-tukb[tukb$CHIP == 0,]
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + MCHIP, data = tp1)
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + MCHIP, data = tp1[tp1$Msize != 1,])


mvdf1<-extractCox(cox1,"Myeloid","MCHIP1","Myeloid","Any CH-MD")
#mvdf2<-extractCox(cox2,"Myeloid","Msize1","Myeloid","CH-MD, VAF<0.1")
mvdf3<-extractCox(cox2,"Myeloid","MCHIP1","Myeloid","CH-MD, VAF>=0.1")
mvdf0<-data.frame(Malignancy = "Myeloid", CHIP = "Myeloid",HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Baseline")
temp2<-rbind(mvdf0,mvdf1,mvdf3)
temp2$at_risk<-c(nrow(tp1[tp1$MCHIP == 0,]),nrow(tp1[tp1$MCHIP == 1,]),nrow(tp1[tp1$Msize == 2,]))
temp2$event<-c(nrow(tp1[tp1$MCHIP == 0 & tp1$Status == 1,]),nrow(tp1[tp1$MCHIP == 1 & tp1$Status ==1,]),nrow(tp1[tp1$Msize == 2 & tp1$Status == 1,]))
temp2$Col<-rep("My",nrow(temp2))
temp2[temp2$VAF == "Baseline",]$Col<-"Baseline"


temp2<-rbind(temp2,ly)
temp2$Col<-as.factor(temp2$Col)


temp2$VAF<-as.character(temp2$VAF)
temp2$VAF<-gsub("CH-LD","L-CHIP",temp2$VAF)
temp2$VAF<-gsub("CH-MD","M-CHIP",temp2$VAF)
temp2$VAF<-gsub("Baseline","No CHIP",temp2$VAF)
temp2$VAF<-gsub("Any ","",temp2$VAF)

## Table
tdf<-temp2
tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("VAF","at_risk","event","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf2<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
		padding = unit(c(2, 2), "mm")))


temp2$VAF<-factor(temp2$VAF,levels=c("L-CHIP, VAF>=0.1","L-CHIP","M-CHIP, VAF>=0.1","M-CHIP", "No CHIP"))
temp2$Col<-factor(temp2$Col,levels=c("Baseline","My","Ly"))
p2 <- ggplot(data=temp2, aes(x=VAF, y=HR, ymin=CI_low, ymax=CI_high,color=Col)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8,shape=15, fatten=3) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c("black",ctList[2],ctList[1]))

F3d<-ggarrange(p2,ptdf2,ncol=2,nrow=1,widths=c(1,1.5))


#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
## UKB 500K data

print ("Myeloid mCA")

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

cnv$LCNV_size<-cnv$LCNV
cnv[cnv$CELL_FRAC >= 0.1 & cnv$LCNV == 1,]$LCNV_size<-2
cnv$MCNV_size<-cnv$MCNV
cnv[cnv$CELL_FRAC >= 0.1 & cnv$MCNV == 1,]$MCNV_size<-2


##===================================================================
## Myeloid

ukbc<-data.frame(fread(paste(figDir,"Fig_1_UKB500_malignancy_df.txt",sep="")))


ukbc$Caucasian<-as.factor(ukbc$Caucasian)
ukbc$ever_smoked<-factor(ukbc$ever_smoked,levels=c("No","Yes"))
ukbc$Age_categ<-factor(ukbc$Age_categ,levels=c("<50","50-59",">60"))

ukbc$LCNV<-rep(0,nrow(ukbc))
ukbc[ukbc$eid %in% cnv[cnv$LCNV == 1,]$eid,]$LCNV<-1
ukbc$MCNV<-rep(0,nrow(ukbc))
ukbc[ukbc$eid %in% cnv[cnv$MCNV == 1,]$eid,]$MCNV<-1

ukbc$CNV<-rep(0,nrow(ukbc))
ukbc[ukbc$eid %in% cnv$eid,]$CNV<-1
ukbc$CNV_size<-ukbc$CNV
ukbc[ukbc$eid %in% cnv[cnv$CELL_FRAC >= 0.1,]$eid,]$CNV_size<-2

## CH category
ukbc$CH<-rep("Control",nrow(ukbc))
ukbc[ukbc$LCNV == 1 & ukbc$MCNV == 0,]$CH<-"Lymphoid"
ukbc[ukbc$MCNV == 1 & ukbc$LCNV == 0,]$CH<-"Myeloid"
ukbc[ukbc$LCNV == 1 & ukbc$MCNV == 1,]$CH <- "Amb"
ukbc[ukbc$LCNV == 0 & ukbc$MCNV == 0 & ukbc$CNV == 1,]$CH <- "Unk"

ukbc$CH<-factor(ukbc$CH,levels=c("Control","Myeloid","Lymphoid","Amb","Unk"))


#=====================================================================
## Analyzing the effect of three groups of CH

tukb<-ukbc
tukb$has_disease<-tukb$Myeloid
tukb$Censor_date<-tukb$Myeloid_date
tukb<-prepCox(tukb)


## CH
tukb$CH<-factor(tukb$CH,levels=c("Control","Myeloid","Lymphoid","Amb","Unk"))
tukb$Caucasian<-as.factor(tukb$Caucasian)


## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + EverSmoke + Sex + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb)
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + EverSmoke + Sex + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb[tukb$CNV_size != 1,])


# ## Tables
mvdf1<-extractCox(cox1,"Myeloid",c("CHMyeloid","CHLymphoid"),c("mCA-MD","mCA-LD"),"Any")
mvdf2<-extractCox(cox1,"Myeloid",c("CHAmb","CHUnk"),c("mCA-AD","mCA-UD"),"Any")
mvdf11<-extractCox(cox2,"Myeloid",c("CHMyeloid","CHLymphoid"),c("mCA-MD","mCA-LD"),"High")
mvdf12<-extractCox(cox2,"Myeloid",c("CHAmb","CHUnk"),c("mCA-AD","mCA-UD"),"High")
#mvdf5<-data.frame(Malignancy = "", CHIP = "Lymphoid", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
BL<-data.frame(Malignancy = "", CHIP = "Baseline", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf2)
tdf$Total<-as.numeric(table(tukb$CH)[c("Control","Myeloid","Lymphoid","Amb","Unk")])
tdf$Positive<-as.numeric(table(tukb[tukb$Status == 1,]$CH)[c("Control","Myeloid","Lymphoid","Amb","Unk")])
rdf<-rbind(mvdf11,mvdf12)
rdf$Total<-as.numeric(table(tukb[tukb$CNV_size != 1,]$CH)[c("Myeloid","Lymphoid","Amb","Unk")])
rdf$Positive<-as.numeric(table(tukb[tukb$Status == 1 & tukb$CNV_size != 1,]$CH)[c("Myeloid","Lymphoid","Amb","Unk")])
rdf$CHIP<-paste(rdf$CHIP,rdf$VAF,sep="_")

tdf<-rbind(tdf,rdf)
rownames(tdf)<-as.character(tdf$CHIP)
tdf<-tdf[c("Baseline","mCA-MD","mCA-MD_High","mCA-LD","mCA-LD_High","mCA-AD","mCA-AD_High","mCA-UD","mCA-UD_High"),]

tdf$CHIP<-as.character(tdf$CHIP)
tdf$CHIP<-gsub("mCA-UD_High","U-mCA, CF>=0.1",tdf$CHIP)
tdf$CHIP<-gsub("mCA-UD","U-mCA",tdf$CHIP)
tdf$CHIP<-gsub("mCA-AD_High","A-mCA, CF>=0.1",tdf$CHIP)
tdf$CHIP<-gsub("mCA-AD","A-mCA",tdf$CHIP)
tdf$CHIP<-gsub("mCA-LD_High","L-mCA, CF>=0.1",tdf$CHIP)
tdf$CHIP<-gsub("mCA-LD","L-mCA",tdf$CHIP)
tdf$CHIP<-gsub("mCA-MD_High","M-mCA, CF>=0.1",tdf$CHIP)
tdf$CHIP<-gsub("mCA-MD","M-mCA",tdf$CHIP)
tdf$CHIP<-gsub("Baseline","No mCA",tdf$CHIP)

tdf$CHIP<-factor(tdf$CHIP,levels=c("U-mCA, CF>=0.1","U-mCA","A-mCA, CF>=0.1","A-mCA","L-mCA, CF>=0.1","L-mCA","M-mCA, CF>=0.1","M-mCA","No mCA"))

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
	scale_color_manual(values=c(rep("grey",2),rep("orange",2),rep(ctList[1],2),rep(ctList[2],2),"black"))+
	scale_y_log10()+
	ggtitle("Myeloid malignancies")


tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1b<-ggarrange(p97,ptdf,ncol=2,nrow=1,widths=c(1,1.5))


#=====================================================================
#=====================================================================
## Lymphoid


tukb<-ukbc
tukb$has_disease<-tukb$Lymphoid
tukb$Censor_date<-tukb$Lymphoid_date
tukb<-prepCox(tukb)


## CH
tukb$CH<-factor(tukb$CH,levels=c("Control","Myeloid","Lymphoid","Amb","Unk"))
tukb$Caucasian<-as.factor(tukb$Caucasian)


## CHIP size
cox1 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb)
cox2 <- coxph(Surv(Year, Status) ~ Age_categ + Sex + EverSmoke + Caucasian + PC1 + PC2 + PC3 + PC4 + PC5 + CH, data = tukb[tukb$CNV_size != 1,])


# ## Tables
mvdf1<-extractCox(cox1,"Myeloid",c("CHMyeloid","CHLymphoid"),c("mCA-MD","mCA-LD"),"Any")
mvdf2<-extractCox(cox1,"Myeloid",c("CHAmb","CHUnk"),c("mCA-AD","mCA-UD"),"Any")
mvdf11<-extractCox(cox2,"Myeloid",c("CHMyeloid","CHLymphoid"),c("mCA-MD","mCA-LD"),"High")
mvdf12<-extractCox(cox2,"Myeloid",c("CHAmb","CHUnk"),c("mCA-AD","mCA-UD"),"High")
#mvdf5<-data.frame(Malignancy = "", CHIP = "Lymphoid", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
BL<-data.frame(Malignancy = "", CHIP = "Baseline", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1,mvdf2)
tdf$Total<-as.numeric(table(tukb$CH)[c("Control","Myeloid","Lymphoid","Amb","Unk")])
tdf$Positive<-as.numeric(table(tukb[tukb$Status == 1,]$CH)[c("Control","Myeloid","Lymphoid","Amb","Unk")])
rdf<-rbind(mvdf11,mvdf12)
rdf$Total<-as.numeric(table(tukb[tukb$CNV_size != 1,]$CH)[c("Myeloid","Lymphoid","Amb","Unk")])
rdf$Positive<-as.numeric(table(tukb[tukb$Status == 1 & tukb$CNV_size != 1,]$CH)[c("Myeloid","Lymphoid","Amb","Unk")])
rdf$CHIP<-paste(rdf$CHIP,rdf$VAF,sep="_")

tdf<-rbind(tdf,rdf)
rownames(tdf)<-as.character(tdf$CHIP)
tdf<-tdf[c("Baseline","mCA-MD","mCA-MD_High","mCA-LD","mCA-LD_High","mCA-AD","mCA-AD_High","mCA-UD","mCA-UD_High"),]

tdf$CHIP<-as.character(tdf$CHIP)
tdf$CHIP<-gsub("mCA-UD_High","U-mCA, CF>=0.1",tdf$CHIP)
tdf$CHIP<-gsub("mCA-UD","U-mCA",tdf$CHIP)
tdf$CHIP<-gsub("mCA-AD_High","A-mCA, CF>=0.1",tdf$CHIP)
tdf$CHIP<-gsub("mCA-AD","A-mCA",tdf$CHIP)
tdf$CHIP<-gsub("mCA-LD_High","L-mCA, CF>=0.1",tdf$CHIP)
tdf$CHIP<-gsub("mCA-LD","L-mCA",tdf$CHIP)
tdf$CHIP<-gsub("mCA-MD_High","M-mCA, CF>=0.1",tdf$CHIP)
tdf$CHIP<-gsub("mCA-MD","M-mCA",tdf$CHIP)
tdf$CHIP<-gsub("Baseline","No mCA",tdf$CHIP)

tdf$CHIP<-factor(tdf$CHIP,levels=c("U-mCA, CF>=0.1","U-mCA","A-mCA, CF>=0.1","A-mCA","L-mCA, CF>=0.1","L-mCA","M-mCA, CF>=0.1","M-mCA","No mCA"))

p87 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="grey") +
	geom_pointrange(size=0.8,shape=15,fatten=3) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c(rep("grey",2),rep("orange",2),rep(ctList[1],2),rep(ctList[2],2),"black"))+
	scale_y_log10()+
	ggtitle("Lymphoid malignancies")


##

tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("CHIP","Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("","Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1c<-ggarrange(p87,ptdf,ncol=2,nrow=1,widths=c(1,1.5))


#==========================================================================
## Final figurewidth
## Max A4 size is 8x11

F1<-ggarrange(F3a, F3d, F1b, F1c, nrow=4,ncol=1,heights=c(1,1,1.8,1.8))


pdf(paste(figDir,"scripts/Final_figs/Fig_ED1.pdf",sep=""),width=8,height=8)
print (F1)
dev.off()



