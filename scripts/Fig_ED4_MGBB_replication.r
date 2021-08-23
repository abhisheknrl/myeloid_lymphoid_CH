

library(data.table)
library(ggplot2)
library(stringr)


mgbData<-"data/2021-06-02_MGB_incident_heme.onc_wCHIP.csv"


lchipFile<-"data/MGBB_LCHIP_calls.txt"
mchipFile<-"data/MGBB_CHIP_call_AN_CG_2020_11_11.txt"
mcaFile<-"data/MGBB_mCA_calls_new.txt"

## M-CHIP and L-CHIP genes
mygeneFile<-"CHIP_calls/CHIP_genes_18June2020.txt"
lygeneFile<-"CHIP_calls/Lymphoid_geneList.txt"

## ICD codes
icdFile<-"meta/heme_malignancies_icd10.txt"
icdFile9<-"meta/heme_malignancies_icd9.txt"


reviewFile<-"data/MGBB_Biobank_review_map.txt"
dnaCollectionFile<-"data/MGBB_DNA_collection.txt"

## Confirmed cases
myCases<-"data/MGBB_Myeloid_confirmed_cases.txt"
lyCases<-"data/MGBB_Lymphoid_confirmed_cases.txt"

figDir<-""

#==============================================================
## DOB File from Aswin (only these cases can be reviewed)

Review<-data.frame(fread(reviewFile))

myValid<-data.frame(fread(myCases))
lyValid<-data.frame(fread(lyCases))

## DNA collection date
dna<-data.frame(fread(dnaCollectionFile))
dna<-dna[dna$SUBJECTID %in% myValid$Biobank.Subject.ID | dna$SUBJECTID %in% lyValid$ID,]

DNA<-as.character(dna$COLLECT_DATE)
names(DNA)<-as.character(dna$SUBJECTID)

## DNA collection date
myValid$DNA<-DNA[as.character(myValid$Biobank.Subject.ID)]
lyValid$DNA<-DNA[as.character(lyValid$ID)]

## Days to diagnosis from DNA collection date
myValid$Diff<-as.numeric(difftime(as.Date(myValid$Diagnosis.Date,format="%m/%d/%Y"), as.Date(myValid$DNA, format="%m/%d/%Y"), units=c("days")))
lyValid$Diff<-as.numeric(difftime(as.Date(lyValid$Date,format="%m/%d/%Y"), as.Date(lyValid$DNA, format="%m/%d/%Y"), units=c("days")))


MYDAYS<-myValid$Diff
LYDAYS<-lyValid$Diff
names(MYDAYS)<-as.character(myValid$Biobank.Subject.ID)
names(LYDAYS)<-as.character(lyValid$ID)


#=================================================================
## Heme ICD codes
hicd<-data.frame(fread(icdFile))
hicd9<-data.frame(fread(icdFile9))


#=================================================================
## CHIP and mCAs

mchip<-data.frame(fread(mchipFile))
mchip<-mchip[mchip$CHIP_GKG_call_v2 == "CHIP",]
mchip$VAF<-round(mchip$AD_alt/(mchip$AD_alt+mchip$AD_ref),3)
lchip<-data.frame(fread(lchipFile))


## mCA
cnv<-data.frame(fread(mcaFile))
cnv<-cnv[cnv$CHR %in% paste("chr",1:22,sep=""),]


mygeneList<-readLines(mygeneFile)
mygeneList<-mygeneList[!mygeneList %in% c("MYD88","STAT5B","STAT3")]

## Annotating CNV data with the genes
### https://ashpublications.org/blood/article/115/14/2731/27247/Copy-neutral-loss-of-heterozygosity-a-novel

myg<-c()

for (i in 1:nrow(cnv)){
	copy<-as.character(cnv$type)[i]
	if (copy != "CN-LOH"){
		myg[i]<-""
		next
	}
	if (cnv$CHR[i] %in% c("X","Y")){
		myg[i]<-""
		next
	}
	gl<-unlist(strsplit(as.character(cnv$Genes)[i],"[,]"))
	I<-intersect(gl,mygeneList)
	if ("DNMT3A" %in% I){
		myg[i]<-"DNMT3A"
	}else if ("TET2" %in% I){
		myg[i]<-"TET2"
	}else if ("JAK2" %in% I){
		myg[i]<-"JAK2"
	}else if ("TP53" %in% I){
		myg[i]<-"TP53"
	}else if ("KIT" %in% I){
		myg[i]<-"KIT"
	}else if ("GNB1" %in% I){
		myg[i]<-"MPL/GNB1"
	}else if ("EP300" %in% I){
		myg[i]<-"EP300"
	}else if ("CBL" %in% I){
		myg[i]<-"CBL"
	}else if ("CTCF" %in% I){
		myg[i]<-"CTCF"
	}else if ("EZH2" %in% I){
		myg[i]<-"EZH2"
	}else{
		myg[i]<-paste(I,collapse=",")
	}
}

cnv$myGenes<-myg


## Lymphoid genes
lygeneList<-readLines(lygeneFile)
lygeneList<-lygeneList[!lygeneList %in% mygeneList]

lyg<-c()
for (i in 1:nrow(cnv)){
	copy<-as.character(cnv$type)[i]
	if (copy != "CN-LOH"){
		lyg[i]<-""
		next
	}
	if (cnv$CHR[i] %in% c("X","Y")){
		lyg[i]<-""
		next
	}
	gl<-unlist(strsplit(as.character(cnv$Genes)[i],"[,]"))
	I<-intersect(gl,lygeneList)
	if ("ATM" %in% I){
		lyg[i]<-"ATM"
	}else if ("CD79B" %in% I){
		lyg[i]<-"CD79B"
	}else if ("CD79A" %in% I){
		lyg[i]<-"CD79A"
	}else if ("NOL9" %in% I){
		lyg[i]<-"NOL9"
	}else if ("ITPKB" %in% I | "NTRK1" %in% I){
		lyg[i]<-"ITPKB/NTRK1"
	}else if ("NCOR2" %in% I){
		lyg[i]<-"NCOR2"
	}else if ("MIR16-1" %in% I){
		lyg[i]<-"MIR16-1"
	}else if ("TCL1A" %in% I){
		lyg[i]<-"TCL1A"
	}else if ("CIITA" %in% I){
		lyg[i]<-"CIITA"
	}else if ("KMT2C" %in% I){
		lyg[i]<-"KMT2C"
	}else if ("NOTCH1" %in% I){
		lyg[i]<-"NOTCH1"
	}else if ("CDKN2A" %in% I){
		lyg[i]<-"CDKN2A"
	}else if ("CHD2" %in% I){
		lyg[i]<-"CHD2"
	}else if ("ZNF217" %in% I){
		lyg[i]<-"ZNF217"
	}else{
		lyg[i]<-paste(I,collapse=",")
	}
}

cnv$lyGenes<-lyg



#=================================================================
## Myeloid

LYG<-c("ATM","CIITA","ITPKB/NTRK1","KMT2C","MIR16-1","NCOR2","NOTCH1")
MYG<-c("CBL","CTCF","EP300","JAK2","TP53","MPL/GNB1")

## canonical myeloid/lymphoid
canLy<-c("gain12q","gain15q","gain17q","gain21q","gain22q13.32","Gain22q13.32","gain2p","gain3q","gain8q","gain9q","del10p","del10q","del11q","del13q","del13q14","del14q","del15q","del17p","del1p","del1q","del22q","del6q","del7q","del8p","tri12","tri18","tri19")
canMy<-c("gain1q","gain21q","gain9p","del12q","del20q","del5q","tri8")

##

manMG<-c("TP53","CTCF","MPL/GNB1","CBL")
manLG<-c("TCL1A","ATM","HEATR3/PLCG2/IRF8")
cnv$LCNV<-rep(0,nrow(cnv))
cnv[cnv$lyGenes %in% c(LYG,manMG) | cnv$myGenes %in% c(LYG,manMG) | cnv$Canonical_mCA %in% canLy,]$LCNV<-1
cnv$MCNV<-rep(0,nrow(cnv))
cnv[cnv$myGenes %in% c(MYG,manLG) | cnv$lyGenes %in% c(MYG,manLG) | cnv$Canonical_mCA %in% canMy,]$MCNV<-1



#=================================================================
## Data file

f<-read.csv(mgbData)
f<-f[!is.na(f$inWES) | !is.na(f$array),]
f<-f[!is.na(f$age) & f$Gender %in% c("Male","Female"),]



f$L_CHIP<-rep(0,nrow(f))
f$M_CHIP<-rep(0,nrow(f))
f$CNV<-rep(0,nrow(f))
f$mCA<-rep("No mCA",nrow(f))

f[f$Biobank.Subject.ID %in% mchip$sample_id,]$M_CHIP<-1
f[f$Biobank.Subject.ID %in% lchip$sample_id,]$L_CHIP<-1
f[f$Biobank.Subject.ID %in% cnv$sample_id,]$CNV<-1

m<-cnv[cnv$MCNV == 1,]
l<-cnv[cnv$LCNV == 1,]
f[f$Biobank.Subject.ID %in% cnv[cnv$MCNV == 1,]$sample_id,]$mCA<-"M-mCA"
f[f$Biobank.Subject.ID %in% cnv[cnv$LCNV == 1,]$sample_id,]$mCA<-"L-mCA"
f[f$Biobank.Subject.ID %in% intersect(m$sample_id,l$sample_id),]$mCA<-"A-mCA"
f[f$Biobank.Subject.ID %in% cnv[!cnv$sample_id %in% c(m$sample_id,l$sample_id),]$sample_id,]$mCA<-"U-mCA"


f<-f[f$Biobank.Subject.ID %in% Review$V1,]

T<-table(f$Biobank.Subject.ID)
f<-f[f$Biobank.Subject.ID %in% names(T[T==1]),]


#==============================================================
#==============================================================
#==============================================================
#==============================================================
## Adding chart reviewed diagnoses

f$MYDAYS<-MYDAYS[as.character(f$Biobank.Subject.ID)]
f$LYDAYS<-LYDAYS[as.character(f$Biobank.Subject.ID)]

## Excluding the prevalent cases from the chart review
excl1<-f[!is.na(f$MYDAYS) & f$MYDAYS < (365.25/2),]$Biobank.Subject.ID
excl2<-f[!is.na(f$LYDAYS) & f$LYDAYS < (365.25/2),]$Biobank.Subject.ID
excl<-c(excl1,excl2)

f<-f[!f$Biobank.Subject.ID %in% excl,]

## Adding the days of follow-up for myeloid and lymphoid
## First myeloid
f1<-f[!is.na(f$MYDAYS),]
f1$Diff<-f1$followup_maxdays - f1$MYDAYS
f1<-f1[f1$Diff >= 0,]
f1$MY<-rep(1,nrow(f1))
f1$MY_t<-f1$MYDAYS/365.25

f2<-f[!f$Biobank.Subject.ID %in% f1$Biobank.Subject.ID,]
f2$MY<-rep(0,nrow(f2))
f2$MY_t<-f2$followup_maxdays/365.25
f<-rbind(f1[,colnames(f2)],f2)

## Now lymphoid
f1<-f[!is.na(f$LYDAYS),]
f1$Diff<-f1$followup_maxdays - f1$LYDAYS
f1<-f1[f1$Diff >= 0,]
f1$LY<-rep(1,nrow(f1))
f1$LY_t<-f1$LYDAYS/365.25

f2<-f[!f$Biobank.Subject.ID %in% f1$Biobank.Subject.ID,]
f2$LY<-rep(0,nrow(f2))
f2$LY_t<-f2$followup_maxdays/365.25
f<-rbind(f1[,colnames(f2)],f2)


#==============================================================
## related individual exclusion

mt<-f
mt<-mt[mt$first_or_SecondDegRelative_toRemove == 0,]
mt<-mt[mt$age >= 18,]

#==============================================================
## WES cohort

mtwes<-mt[!is.na(mt$inWES),]
mtwes<-mtwes[!is.na(mtwes$heme.malig_t),]
mtwes<-mtwes[mtwes$heme.malig_t >= 0.5 & mtwes$myeloid_t >= 0.5 & mtwes$lymphoid_t >= 0.5,]

print ("Age distribution")
print (summary(mtwes$age))

p0<-ggplot(mtwes,aes(age,L_CHIP))+
	geom_smooth(method="gam",formula = y ~ s(x,bs="cr"),color="blue",fill="blue",alpha=0.2)+
	geom_smooth(data=mtwes,aes(age,M_CHIP),method="gam",formula = y ~ s(x,bs="cr"),color="red",fill="red",alpha=0.2)+
	theme_linedraw()+
	theme(panel.grid=element_blank(),
		panel.grid.major.y=element_line(color="grey"))+
	xlab("Age")+
	ylab("Proportion with CHIP")


#==============================================================
## Gene distribution

lchip<-lchip[lchip$sample_id %in% mtwes$Biobank.Subject.ID,]
mchip<-mchip[mchip$sample_id %in% mtwes$Biobank.Subject.ID,]


getFig<-function(temp1,mycol,XLAB,isMY="F"){
	#==============================
	## genes
	T<-table(as.character(temp1$Hugo_Symbol))
	T<-T[order(T,decreasing=TRUE)]
	K<-1
	for (gene in names(T)[1:25]){
		temp<-temp1[temp1$Hugo_Symbol == gene,]
		##
		tempdf<-data.frame(Gene = gene,
			Count = as.numeric(T[gene]),
			N = nrow(temp))
		if (K==1){
			K<-2
			tdf<-tempdf
		}else{
			tdf<-rbind(tdf,tempdf)
		}
	}
	##
	tdf$Gene<-factor(tdf$Gene,levels=names(T[order(T,decreasing=TRUE)]))
	sf2d<-ggplot(tdf,aes(y=Count,x=Gene))+
		geom_bar(stat="identity",color=NA,fill=mycol)+
		theme_classic()+
		xlab(XLAB)+
		ylab("Number of mutations")+
		scale_y_log10()+
		theme(legend.position="none",
			panel.grid.major.y=element_line(color="gray80"))+
		theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=rel(0.8)),
			axis.text.y=element_text(color="black"),
			axis.ticks = element_line(color="black"))
	if (isMY == "T"){
		sf2d<-ggplot(tdf,aes(y=Count,x=Gene))+
		#geom_hline(yintercept = 30,color="grey")+
		geom_bar(stat="identity",color=NA,fill=mycol)+
		theme_classic()+
		xlab(XLAB)+
		ylab("Number of mutations")+
		scale_y_log10()+
		theme(legend.position="none",
			panel.grid.major.y=element_line(color="gray80"))+
		theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=rel(0.8)),
			axis.text.y=element_text(color="black"),
			axis.ticks = element_line(color="black"))
	}
	return (list(sf2d,tdf))
}

F1b<-getFig(lchip[lchip$sample_id %in% mtwes$Biobank.Subject.ID,],"blue","L-CHIP genes")[[1]]
F1a<-getFig(mchip[mchip$sample_id %in% mtwes$Biobank.Subject.ID,],"red","M-CHIP genes","T")[[1]]

#==============================================================
#==============================================================
### Heme malignancies

## insufficient sample size


#==============================================================
#==============================================================
#==============================================================
#==============================================================
## SNP-array cohort

mtsnp<-mt[!is.na(mt$array),]
mtsnp<-mtsnp[!is.na(mtsnp$heme.malig_t),]
mtsnp<-mtsnp[mtsnp$heme.malig_t >= 0.5 & mtsnp$myeloid_t >= 0.5 & mtsnp$lymphoid_t >= 0.5,]

## Exclude array that does not have any mCA detected (Not sure about it)
mtsnp<-mtsnp[mtsnp$array != "mega",]


## cnv
cnv<-cnv[cnv$sample_id %in% mtsnp$Biobank.Subject.ID,]


p1<-ggplot(mtsnp[mtsnp$mCA %in% c("No mCA","U-mCA"),],aes(age,CNV))+
	geom_smooth(method="gam",formula = y ~ s(x,bs="cr"),color="grey",fill="grey",alpha=0.2)+
	geom_smooth(data=mtsnp[mtsnp$mCA %in% c("No mCA","M-mCA"),],aes(age,CNV),method="gam",formula = y ~ s(x,bs="cr"),color="red",fill="red",alpha=0.2)+
	geom_smooth(data=mtsnp[mtsnp$mCA %in% c("No mCA","L-mCA"),],aes(age,CNV),method="gam",formula = y ~ s(x,bs="cr"),color="blue",fill="blue",alpha=0.2)+
	geom_smooth(data=mtsnp[mtsnp$mCA %in% c("No mCA","A-mCA"),],aes(age,CNV),method="gam",formula = y ~ s(x,bs="cr"),color="orange",fill="orange",alpha=0.2)+
	theme_linedraw()+
	theme(panel.grid.minor=element_blank(),
		panel.grid.major.x=element_blank())+
	xlab("Age")+
	ylab("Proportion with CHIP")


#=====================================================================
## mtsnp association
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

## Myeloid malignancies

library(survival)
library(survminer)


#t<-mtsnp[mtsnp$myeloid_outcome == 0 | mtsnp$Biobank.Subject.ID %in% myValid$Biobank.Subject.ID,]
t<-mtsnp
q95<-5.5
t[t$MY_t > q95,]$MY<-0
t[t$MY_t > q95,]$MY_t<-q95

t$mCA<-factor(as.character(t$mCA),levels=c("No mCA","M-mCA","L-mCA","A-mCA","U-mCA"))
## Myeloid maliganncies
YLAB<-"Cumulative incidence of\nmyeloid malignancies"
fit<-survfit(Surv(MY_t,MY) ~ mCA,data=t)
p98<-ggsurvplot(fit,data=t,fun="event",
		legend.title = "CH", legend=c(0.2,0.75), legend.labs=c("Baseline","M-mCA","L-mCA","A-mCA","U-mCA"),
		palette=c("black","red","blue","orange","grey"),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")

cox1 <- coxph(Surv(MY_t, MY) ~ age + age2 + Gender + PC1 + PC2 + PC3 + PC4 + PC5 + mCA, data = t)

mvdf1<-extractCox(cox1, "Myeloid",c("mCAM-mCA","mCAL-mCA","mCAA-mCA","mCAU-mCA"),c("M-mCA","L-mCA","A-mCA","U-mCA"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "No mCA", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf1)

tdf$Total<-as.numeric(table(t$mCA)[c("No mCA","M-mCA","L-mCA","A-mCA","U-mCA")])
tdf$Positive<-as.numeric(table(t[t$MY == 1,]$mCA)[c("No mCA","M-mCA","L-mCA","A-mCA","U-mCA")])

if (nrow(tdf[tdf$Positive <= 1,]) >= 1){
	tdf[tdf$Positive <= 1,]$HR<-NA
	tdf[tdf$Positive <= 1,]$CI_high<-NA
	tdf[tdf$Positive <= 1,]$CI_low<-NA
	tdf[tdf$Positive <= 1,]$P<-NA
}

tdf$CHIP<-factor(tdf$CHIP,levels=c("U-mCA","A-mCA","L-mCA","M-mCA","No mCA"))
p99 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8, shape=15, fatten=3) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c("grey","orange","blue","red","black"))


tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1c<-p98
F1d<-ggarrange(ggarrange(p99,ptdf,ncol=2,nrow=1,widths=c(1,1.5)),
	nrow=2,heights=c(2,1))



#=====================================================================
#=====================================================================
## Lymphoid malignancies

t<-mtsnp

MM<-c("Multiple myeloma","Myeloma","Myeloma (no path)","Smouldering myeloma")
t<-t[!t$Biobank.Subject.ID %in% lyValid[lyValid$Diagnosis %in% MM,]$ID,]

q95<-5.5
t[t$LY_t > q95,]$LY<-0
t[t$LY_t > q95,]$LY_t<-q95

t$mCA<-factor(as.character(t$mCA),levels=c("No mCA","M-mCA","L-mCA","A-mCA","U-mCA"))
## Myeloid maliganncies
YLAB<-"Cumulative incidence of\nlymphoid malignancies"
fit<-survfit(Surv(LY_t,LY) ~ mCA,data=t)
p88<-ggsurvplot(fit,data=t,fun="event",
		legend.title = "CH", legend=c(0.2,0.75), legend.labs=c("Baseline","M-mCA","L-mCA","A-mCA","U-mCA"),
		palette=c("black","red","blue","orange","grey"),
		risk.table = FALSE,
		ggtheme = theme_classic(),
		censor = FALSE,
		xlab = "Years", ylab = YLAB,
		break.time.by = 1, surv.scale="percent")

cox11 <- coxph(Surv(LY_t, LY) ~ age + age2 + Gender + PC1 + PC2 + PC3 + PC4 + PC5 + mCA, data = t)


mvdf11<-extractCox(cox11, "Lymphoid",c("mCAM-mCA","mCAL-mCA","mCAA-mCA","mCAU-mCA"),c("M-mCA","L-mCA","A-mCA","U-mCA"),"Any")
BL<-data.frame(Malignancy = "", CHIP = "No mCA", HR = 1, CI_low = 1, CI_high = 1, P = 1, VAF = "Any")
tdf<-rbind(BL,mvdf11)

tdf$Total<-as.numeric(table(t$mCA)[c("No mCA","M-mCA","L-mCA","A-mCA","U-mCA")])
tdf$Positive<-as.numeric(table(t[t$LY == 1,]$mCA)[c("No mCA","M-mCA","L-mCA","A-mCA","U-mCA")])

if (nrow(tdf[tdf$Positive <= 1,]) >= 1){
	tdf[tdf$Positive <= 1,]$HR<-NA
	tdf[tdf$Positive <= 1,]$CI_high<-NA
	tdf[tdf$Positive <= 1,]$CI_low<-NA
	tdf[tdf$Positive <= 1,]$P<-NA
}

tdf$CHIP<-factor(tdf$CHIP,levels=c("U-mCA","A-mCA","L-mCA","M-mCA","No mCA"))
p89 <- ggplot(data=tdf, aes(x=CHIP, y=HR, ymin=CI_low, ymax=CI_high,color=CHIP)) +
	geom_hline(yintercept=1, lty=2,col="black") +
	geom_pointrange(size=0.8, shape=15, fatten=3) +
	coord_flip() + 
	xlab("") + ylab("Hazard ratio (95% CI)") +
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.line.y=element_blank(),
		axis.text=element_text(color="black"),
		axis.ticks.y=element_blank(),
		legend.position="none")+
	scale_color_manual(values=c("grey","orange","blue","red","black"))


tdf$HR_ci<-paste(round(tdf$HR,2)," (",round(tdf$CI_low,2),", ",round(tdf$CI_high,2),")",sep="")
tdf$P<-format(tdf$P, format = "e", digits = 3)
tdf1<-tdf[,c("Total","Positive","HR_ci","P")]
colnames(tdf1)<-c("Total","Event","HR (95% CI)","p-value")
ptdf1<-ggtexttable(tdf1,
	rows=NULL,
	theme=ttheme(base_style="classic",base_size = 9,
	padding = unit(c(2, 2), "mm")))


F1e<-p88
F1f<-ggarrange(ggarrange(p89,ptdf1,ncol=2,nrow=1,widths=c(1,1.5)),
	nrow=2,heights=c(2,1))


## output figure

P<-ggarrange(ggarrange(F1a,F1b,ncol=1,nrow=2),
	ggarrange(p0, p1, nrow=1, ncol= 2),
	ggarrange(F1d, F1c$plot, nrow=1, ncol = 2, widths=c(1.5,1)),
	ggarrange(F1f, F1e$plot, nrow=1, ncol = 2, widths=c(1.5,1)),
	nrow=4)



pdf(paste(figDir,"scripts/Final_figs/ED_Fig_4_MGBB.pdf",sep=""),width=8,height=10)
print (P)
dev.off()


mt1<-mt[mt$myeloid_outcome == 1 | mt$lymphoid_outcome == 1,]

eidList<-unique(c(mchip$sample_id,lchip$sample_id,cnv$sample_id,mt1$Biobank.Subject.ID))


EL<-paste("S",1:length(eidList),sep="")
names(EL)<-as.character(eidList)

mchip$SID<-EL[as.character(mchip$sample_id)]
lchip$SID<-EL[as.character(lchip$sample_id)]
cnv$SID<-EL[as.character(cnv$sample_id)]
mt1$SID<-EL[as.character(mt1$Biobank.Subject.ID)]


write.table(mchip,file=paste(figDir,"scripts/data/MGBB_M_CHIP.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
write.table(lchip,file=paste(figDir,"scripts/data/MGBB_L_CHIP.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
write.table(cnv,file=paste(figDir,"scripts/data/MGBB_mCA_classified.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")



