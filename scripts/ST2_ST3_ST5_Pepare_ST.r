
##===================================================================
##===================================================================


library(ggplot2)
library(data.table)
library(ggpubr)


setwd("")
hemeICDFile<-"meta/heme_malignancies_icd10.txt"
figDir<-""

ukbc1<-data.frame(fread(paste(figDir,"data/UKB_500K.txt",sep="")))



chip<-read.table(paste(figDir,"1_CHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
mychip<-read.table(paste(figDir,"1_MCHIP_UKB.txt",sep=""),head=TRUE,sep="\t")

chip<-read.table(paste(figDir,"1_CHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
mychip<-read.table(paste(figDir,"1_MCHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
cnv<-data.frame(fread(paste(figDir,"data/UKB_500K_CNV.txt",sep="")))
cnv[is.na(cnv$CELL_FRAC),]$CELL_FRAC<-0.001
cnv<-cnv[cnv$CHR %in% 1:22,]
cnv<-cnv[cnv$eid %in% ukbc1$eid,]


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
#=====================================================================
## FIG

ukbc<-data.frame(fread(paste(figDir,"scripts/data/First_heme_diagnosis.txt",sep="")))


#=====================================================================
## eidlist
eidlist<-unique(c(chip$eid,mychip$eid,cnv$eid,ukbc$eid))
EL<-paste("S",1:length(eidlist),sep="")
names(EL)<-as.character(eidlist)

chip$ID<-EL[as.character(chip$eid)]
mychip$ID<-EL[as.character(mychip$eid)]
cnv$ID<-EL[as.character(cnv$eid)]
ukbc$ID<-EL[as.character(ukbc$eid)]

#=====================================================================
## Write tables
## ST2 -  L-CHIP
chip$MCHIP<-rep(0,nrow(chip))
chip[chip$eid %in% mychip$eid,]$MCHIP<-1
chip$mCA<-rep(0,nrow(chip))
chip[chip$eid %in% cnv$eid,]$mCA<-1
CN<-c("ID","eid","CHR","POS","ref","alt","AD_ref","AD_alt","Hugo_Symbol","Variant_Classification","Variant_Type","Annotation_Transcript","cDNA_Change","Protein_Change","Category","MCHIP","mCA","Overlapping_mCA")

overlapping<-c()
for (i in 1:nrow(chip)){
	EID=chip$eid[i]
	if (EID %in% cnv$eid){
		CHR<-chip$CHR[i]
		Gene<-as.character(chip$Hugo_Symbol[i])
		t<-cnv[cnv$eid == EID,]
		genelist<-unlist(strsplit(paste(t[t$CHR == CHR,]$Genes,collapse=","),","))
		if (Gene %in% genelist){
			#print (EID)
			overlapping<-c(overlapping,EID)
		}
	}
}

chip$Overlapping_mCA<-rep(0,nrow(chip))
chip[chip$eid %in% overlapping,]$Overlapping_mCA<-1

t<-chip[,CN]
write.table(t,file = paste(figDir,"scripts/Final_figs/ST_2_L-CHIP.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

#=====================================================================
## ST3 - MCHIP
mychip$CHIP<-rep(0,nrow(mychip))
mychip[mychip$eid %in% chip$eid,]$CHIP<-1
mychip$mCA<-rep(0,nrow(mychip))
mychip[mychip$eid %in% cnv$eid,]$mCA<-1
CN<-c("ID","eid","CHR","POS","ref","alt","AD_ref","AD_alt","Hugo_Symbol","Variant_Classification","Variant_Type","Annotation_Transcript","cDNA_Change","Protein_Change","CHIP","mCA","Overlapping_mCA")

overlapping<-c()
for (i in 1:nrow(mychip)){
	EID=mychip$eid[i]
	if (EID %in% cnv$eid){
		CHR<-mychip$CHR[i]
		Gene<-as.character(mychip$Hugo_Symbol[i])
		t<-cnv[cnv$eid == EID,]
		genelist<-unlist(strsplit(paste(t[t$CHR == CHR,]$Genes,collapse=","),","))
		if (Gene %in% genelist){
			#print (EID)
			overlapping<-c(overlapping,EID)
		}
	}
}

mychip$Overlapping_mCA<-rep(0,nrow(mychip))
mychip[mychip$eid %in% overlapping,]$Overlapping_mCA<-1

t<-mychip[,CN]
write.table(t,file = paste(figDir,"scripts/Final_figs/ST_3_M-CHIP.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

#=====================================================================
## ST4 - Heme malignancy cases

CN<-c("ID","eid","CHIP","MCHIP","LCNV","MCNV","CNV","LCHIP_var","MCHIP_var","CNV_var","ICD")

chip$VAR<-paste(chip$Hugo_Symbol,chip$Protein_Change,sep="_")
mychip$VAR<-paste(mychip$Hugo_Symbol,mychip$Protein_Change,sep="_")
cnv$VAR<-paste(cnv$Canonical_CA,cnv$myGenes,cnv$lyGenes,sep="_")

LV<-c()
MV<-c()
CV<-c()
for (i in 1:nrow(ukbc)){
	EID<-ukbc$eid[i]
	LV[as.character(EID)]<-paste(chip[chip$eid == EID,]$VAR,collapse=";")
	MV[as.character(EID)]<-paste(mychip[mychip$eid == EID,]$VAR,collapse=";")
	CV[as.character(EID)]<-paste(cnv[cnv$eid == EID,]$VAR,collapse=";")
}

print (head(LV))

ukbc$LCHIP_var<-LV[as.character(ukbc$eid)]
ukbc$MCHIP_var<-MV[as.character(ukbc$eid)]
ukbc$CNV_var<-CV[as.character(ukbc$eid)]


t<-ukbc[,CN]
write.table(t,file = paste(figDir,"scripts/Final_figs/ST_4_Heme_malignancies.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

#=====================================================================
## ST 6
CN<-c("ID","eid","CHR","START_MB","END_MB","SIZE_MB","BAF..SE.","LRR..SE.","COPY_CHANGE","CELL_FRAC","Q_VALUE","START_MB_RANGE","END_MB_RANGE","Canonical_CA","myGenes","lyGenes","LCNV","MCNV")

t<-cnv[,CN]
write.table(t,file = paste(figDir,"scripts/Final_figs/ST_6_mCA.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

