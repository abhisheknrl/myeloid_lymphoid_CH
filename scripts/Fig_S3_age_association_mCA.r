

## Forest plots for myeloid and lymphoid malignancies

library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)
library(mgcv)

library(survival)
library(survminer)
library(gridExtra)

ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")



setwd("")
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"


#=====================================================================
## Reload data

ukb1<-data.frame(fread(paste(figDir,"1_UKB_data.txt",sep="")))


##===================================================================
##===================================================================

## PV + mCA

ukb1$MCH<-rep(0,nrow(ukb1))
ukb1[ukb1$MCNV == 1 | ukb1$MCHIP == 1,]$MCH<-1
ukb1$LCH<-rep(0,nrow(ukb1))
ukb1[ukb1$LCNV == 1 | ukb1$CHIP == 1,]$LCH<-1


p31<-ggplot(ukb1,aes(age_attended_assessment_center,MCH))+
	geom_smooth(method="gam",formula = y ~ s(x,bs="cr"),color=ctList[2],fill=ctList[2],alpha=0.2)+
	geom_smooth(data=ukb1,aes(age_attended_assessment_center,LCH),
		method="gam",formula = y ~ s(x,bs="cr"),color=ctList[1],fill=ctList[1],alpha=0.2)+
		theme_linedraw()+
	theme(panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_line(color="grey"),
		panel.grid.minor=element_blank(),
		legend.position = c(0.8, 0.95),
		legend.key.size = unit(0.3, 'cm'),
		legend.title = element_blank(),
		axis.text=element_text(color="black"))+
	xlab("Age")+
	ylab("Proportion with CH")+
	ggtitle("WES cohort")



##===================================================================
##===================================================================
##===================================================================
## mCA data

ukb<-data.frame(fread(paste(figDir,"data/UKB_500K.txt",sep="")))

ukb$tempMCNV<-ukb$MCNV
ukb[ukb$LCNV == 1,]$tempMCNV<-0
ukb$tempLCNV<-ukb$LCNV
ukb[ukb$MCNV == 1,]$tempLCNV<-0

ukb$tempMLCNV<-ukb$MCNV
ukb[ukb$LCNV == 0,]$tempMLCNV<-0

ukb$tempnCNV<-ukb$CNV
ukb[ukb$MCNV == 1,]$tempnCNV<-0
ukb[ukb$LCNV == 1,]$tempnCNV<-0



p11<-ggplot(ukb,aes(Age,CNV))+
	geom_smooth(method="gam",formula = y ~ s(x,bs="cr"),color="black",fill="black",alpha=0.2)+
	geom_smooth(data=ukb,aes(Age,tempMCNV),method="gam",formula = y ~ s(x,bs="cr"),color=ctList[2],fill=ctList[2],alpha=0.2)+
	geom_smooth(data=ukb,aes(Age,tempLCNV),method="gam",formula = y ~ s(x,bs="cr"),color=ctList[1],fill=ctList[1],alpha=0.2)+
	geom_smooth(data=ukb,aes(Age,tempMLCNV),method="gam",formula = y ~ s(x,bs="cr"),color="orange",fill="orange",alpha=0.2)+
	geom_smooth(data=ukb,aes(Age,tempnCNV),method="gam",formula = y ~ s(x,bs="cr"),color="grey",fill="grey",alpha=0.2)+
	theme_linedraw()+
	theme(panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_line(color="grey"),
		panel.grid.minor=element_blank(),
		legend.position = c(0.8, 0.95),
		legend.key.size = unit(0.3, 'cm'),
		legend.title = element_blank(),
		axis.text=element_text(color="black"))+
	xlab("Age")+
	ylab("Proportion with mCA")+
	ggtitle("SNP-array cohort")



##===================================================================

F<-ggarrange(p11,p31,ncol=2,nrow=1)


pdf(paste(figDir,"scripts/Final_figs/Fig_S3.pdf",sep=""),width=8,height=4)
print (F)
dev.off()
