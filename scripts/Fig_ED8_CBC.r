

## CBC enrichment

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



#=====================================================================
## Filtering out related samples and non-caucasian individuals

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

## Factorize the variables
ukbc$CHIP<-as.factor(ukbc$CHIP)
ukbc$MCHIP<-as.factor(ukbc$MCHIP)
ukbc$CHIP_size<-as.factor(ukbc$CHIP_size)
ukbc$Msize<-as.factor(ukbc$Msize)
ukbc$Caucasian<-as.factor(ukbc$Caucasian)

#=====================================================================
#===================================================================================
## Enrichment analysis for blood counts

library(ggrepel)

cbcAssociation<-function(tukb,cnList,COL){
	dict<-c("ALC","WBC","Lymphs%","Neuts%","PLT","PCT","ANC","RDW","Monos","Eos%","Monos%","Baso")
	names(dict)<-c("lymphocyte_count","white_blood_cell_count","lymphocyte_percentage","neutrophil_percentage","platelet_count","platelet_crit","neutrophil_count","rbc_distribution_width","monocyte_count","eosinophil_percentage","monocyte_percentage","basophil_count")
	for (i in 1:length(cnList)){
		CN<-cnList[i]
		print (CN)
		temp<-data.frame(eid = tukb$eid,
			Sex = tukb$genetic_sex,
			Age = tukb$age_attended_assessment_center,
			Smoking = tukb$ever_smoked,
			Caucasian = tukb$Caucasian,
			CHIP = tukb$CHIP)#,
			#CHIP_size = tukb$CHIP_size)
		##
		temp$CHIP<-as.factor(temp$CHIP)
		#temp$CHIP_size<-as.factor(temp$CHIP_size)
		temp$Smoking<-factor(temp$Smoking,levels=c("No","Yes"))
		temp$Age_sq<-temp$Age*temp$Age
		## Adding the parameter of interest
		temp$Value<-tukb[,c(CN)]
		temp<-temp[!is.na(temp$Value),]
		## Splitting by sex before normalizing
		mt<-temp[temp$Sex == "Male",]
		ft<-temp[temp$Sex == "Female",]
		mt$Norm<-(mt$Value - mean(mt$Value))/sd(mt$Value)
		ft$Norm<-(ft$Value - mean(ft$Value))/sd(ft$Value)
		##
		#mod1<-lm(CHIP ~ Age+Age_sq+Smoking + Caucasian + Norm, data=mt,family="binomial")
		#mod2<-lm(CHIP ~ Age+Age_sq+Smoking + Caucasian + Norm, data=ft,family="binomial")
		mod1<-lm(Norm ~ Age+Age_sq+Smoking + Caucasian + CHIP, data=mt)
		mod2<-lm(Norm ~ Age+Age_sq+Smoking + Caucasian + CHIP, data=ft)
		#mod3<-lm(Norm ~ Sex+Age+Smoking+CHIP, data=temp)
		#print (summary(mod))
		Est<-summary(mod1)$coefficients[6,1]
		SE<-summary(mod1)$coefficients[6,2]
		P<-summary(mod1)$coefficients[6,4]
		temp1<-data.frame(Feature = "CHIP",
			Outcome = CN,
			Covariates = "Age, Smoking",
			Sex = "Male",
			Estimate = Est,
			se = SE,
			P = P)
		Est<-summary(mod2)$coefficients[6,1]
		SE<-summary(mod2)$coefficients[6,2]
		P<-summary(mod2)$coefficients[6,4]
		temp2<-data.frame(Feature = "CHIP",
			Outcome = CN,
			Covariates = "Age, Smoking",
			Sex = "Female",
			Estimate = Est,
			se = SE,
			P = P)
		if (i ==1){
			tdf<-rbind(temp1,temp2)
		}else{
			tdf<-rbind(tdf,temp1,temp2)
		}
	}
	## q
	#tdf1<-tdf[tdf$Sex == "Both",]
	#tdf1$q<-p.adjust(tdf1$P,method="fdr")
	tdf2<-tdf[tdf$Sex == "Male",]
	tdf2$q<-p.adjust(tdf2$P,method="fdr")
	tdf3<-tdf[tdf$Sex == "Female",]
	tdf3$q<-p.adjust(tdf3$P,method="fdr")


	tdf<-rbind(tdf2,tdf3)
	tdf$Odds<-exp(tdf$Estimate)
	tdf$Upper_ci<-exp(tdf$Estimate+1.96*tdf$se)
	tdf$Lower_ci<-exp(tdf$Estimate-1.96*tdf$se)
	tdf<-tdf[order(tdf$Odds,decreasing=TRUE),]
	tdf$Outcome2<-dict[as.character(tdf$Outcome)]

	###
	fp1 <- ggplot(data=tdf, aes(x=-log10(q), y=Estimate)) +
		geom_vline(xintercept=-log10(0.05), lty=2) +
		geom_point(color = "grey",size=1) +
		coord_flip() + 
		xlab("-log10 (p.adjust)") + ylab("Effect size") +
		theme_linedraw()+
		theme(panel.grid=element_blank(),
			legend.position="none")+
		geom_point(data = tdf[tdf$q < 0.05,],aes(x=-log10(q), y=Estimate),color=COL,size=2)+
		geom_text_repel(data = tdf[tdf$q < 0.05,],aes(x=-log10(q), y=Estimate,label=Outcome2,color = Sex), max.overlaps=20)+
		scale_color_manual(values=c("darkorange","darkorchid3"))

	return (list(tdf,fp1))
}


#cnList<-c("white_blood_cell_count","haemoglobin_concentration","red_blood_cell_count","hematocrit_percentage","mean_corpuscular_volumne_mcv","mean_corpuscular_haemoglobin","mean_corpuscular_haemoglobin_concentration","rbc_distribution_width","platelet_count","platelet_crit","mean_platelet_volume","platelet_distribution_width","lymphocyte_count","monocyte_count","neutrophil_count","eosinophil_count","basophil_count","nucleated_rbc_count","lymphocyte_percentage","monocyte_percentage","neutrophil_percentage","eosinophil_percentage","basophil_percentage","nucleated_rbc_percentage","reticulocyte_percentage","reticulocyte_count","mean_reticulocyte_volume","mean_sphered_cell_volume","immature_reticulocyte_fraction","high_light_scatter_reticulocyte_percentage","high_light_scatter_reticulocyte_count")
cnList<-c("white_blood_cell_count","haemoglobin_concentration","red_blood_cell_count","hematocrit_percentage","mean_corpuscular_volumne_mcv","mean_corpuscular_haemoglobin","mean_corpuscular_haemoglobin_concentration","rbc_distribution_width","platelet_count","platelet_crit","lymphocyte_count","monocyte_count","neutrophil_count","eosinophil_count","basophil_count","lymphocyte_percentage","monocyte_percentage","neutrophil_percentage","eosinophil_percentage","basophil_percentage")

print (cnList)
print (length(cnList))


tukb<-ukbc
p1<-cbcAssociation(tukb,cnList,ctList[1])[[2]]

### For myeloid CHIP
tukb<-ukbc
tukb$CHIP<-tukb$MCHIP
p2<-cbcAssociation(tukb,cnList,ctList[2])[[2]]

#====================================================
## UKB 500

ukbc<-data.frame(fread(paste(figDir,"data/UKB_500K.txt",sep="")))
ukbc$Caucasian<-rep(1,nrow(ukbc))
ukbc[is.na(ukbc$genetic_ethnic_grouping),]$Caucasian<-0

ukbc$LCNV<-as.factor(ukbc$LCNV)
ukbc$MCNV<-as.factor(ukbc$MCNV)
ukbc$CNV<-as.factor(ukbc$CNV)
ukbc$Caucasian<-as.factor(ukbc$Caucasian)

## LCNV
tukb<-ukbc[ukbc$LCNV == 0 | ukbc$MCNV == 0,]
tukb$CHIP<-tukb$LCNV
p3<-cbcAssociation(tukb,cnList,ctList[1])[[2]]

### For MCNV
tukb<-ukbc[ukbc$LCNV == 0 | ukbc$MCNV == 0,]
tukb$CHIP<-tukb$MCNV
p4<-cbcAssociation(tukb,cnList,ctList[2])[[2]]

### For MCNV
tukb1<-ukbc[ukbc$LCNV == 0 & ukbc$MCNV == 0,]
tukb2<-ukbc[ukbc$LCNV == 1 & ukbc$MCNV == 1,]
tukb<-rbind(tukb1,tukb2)
tukb$CHIP<-tukb$MCNV
p5<-cbcAssociation(tukb,cnList,"orange")[[2]]

## All other CNVs
tukb<-ukbc[ukbc$LCNV == 0 & ukbc$MCNV == 0,]
tukb$CHIP<-tukb$CNV
p6<-cbcAssociation(tukb,cnList,"black")[[2]]




F1a<-ggarrange(p2,p4,ncol=2,nrow=1)
F1b<-ggarrange(p1,p3,ncol=2,nrow=1)
F1c<-ggarrange(p5,p6,ncol=2,nrow=1)


F<-ggarrange(F1a,F1b,F1c,ncol=1,nrow=3)

#====================================================
#====================================================
#====================================================

pdf(paste(figDir,"scripts/Final_figs/Fig_ED8.pdf",sep=""),width=9,height=9)
print (F)
dev.off()


setEPS()
postscript(paste(figDir,"scripts/Final_figs/Fig_ED8.eps",sep=""),width=9,height=9)
print (F)
dev.off()


