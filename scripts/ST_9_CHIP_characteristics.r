


library(ggplot2)
library(data.table)


setwd("")
hemeICDFile<-"meta/heme_malignancies_icd10.txt"
figDir<-""


##===================================================================
##ICD10codes
hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")
hicdlist<-as.character(hicd[hicd$Group %in% c("Myeloid","Lymphoid","Unspecified"),]$ICD)


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

#=====================================================================
## Reload data

ukbc<-data.frame(fread(paste(figDir,"1_UKB_data.txt",sep="")))

getDF<-function(ukbc,Group){
	N<-nrow(ukbc)
	Age<-as.numeric(quantile(ukbc$age_attended_assessment_center,c(0.25,0.5,0.75)))
	Sex<-nrow(ukbc[ukbc$Sex == "Female",])
	Smoked <-nrow(ukbc[ukbc$ever_smoked == "Yes",])
	Caucasian<-nrow(ukbc[!is.na(ukbc$genetic_ethnic_grouping),])
	## Normal CBC
	ukbc<-ukbc[(!is.na(ukbc$red_blood_cell_count)) & (!is.na(ukbc$platelet_count)) & (!is.na(ukbc$neutrophil_count)) & (!is.na(ukbc$lymphocyte_count)),]
	CBC<-nrow(ukbc[ukbc$lymphocyte_count < ALC_upper & ukbc$lymphocyte_count >ALC_lower & ukbc$neutrophil_count < ANC_upper & ukbc$neutrophil_count > ANC_lower & ukbc$red_blood_cell_count < RBC_upper & ukbc$red_blood_cell_count > RBC_lower & ukbc$platelet_count < PLT_upper & ukbc$platelet_count > PLT_lower,])
	Anemia<-nrow(ukbc[ukbc$neutrophil_count <= ANC_lower | ukbc$red_blood_cell_count <= RBC_lower | ukbc$platelet_count <= PLT_lower,])
	Cytosis<-nrow(ukbc[ukbc$neutrophil_count >= ANC_upper | ukbc$red_blood_cell_count >= RBC_upper | ukbc$platelet_count >= PLT_upper,])
	HALC<-nrow(ukbc[ukbc$lymphocyte_count >= ALC_upper,])
	LALC<-nrow(ukbc[ukbc$lymphocyte_count <= ALC_lower,])
	## Create a data frame
	tdf<-data.frame(Group = Group,
		N = N,
		Age = paste(Age[1]," (",Age[2],", ",Age[3],")",sep=""),
		Female = paste(Sex," (",round(Sex*100/N,1),"%)",sep=""),
		Smoked = paste(Smoked," (",round(Smoked*100/N,1),"%)",sep=""),
		Caucasian = paste(Caucasian," (",round(Caucasian*100/N,1),"%)",sep=""),
		CBC = paste(CBC," (",round(CBC*100/N,1),"%)",sep=""),
		Anemia = paste(Anemia," (",round(Anemia*100/N,1),"%)",sep=""),
		Cytosis = paste(Cytosis," (",round(Cytosis*100/N,1),"%)",sep=""),
		LALC = paste(LALC," (",round(LALC*100/N,1),"%)",sep=""),
		HALC = paste(HALC," (",round(HALC*100/N,1),"%)",sep=""))
	return (tdf)
}

tdf1<-getDF(ukbc[ukbc$CHIP == 0 & ukbc$MCHIP == 0,],"No CHIP")
tdf2<-getDF(ukbc[ukbc$CHIP == 1 & ukbc$MCHIP == 0,],"L-CHIP")
tdf3<-getDF(ukbc[ukbc$CHIP == 0 & ukbc$MCHIP == 1,],"M-CHIP")
tdf4<-getDF(ukbc[ukbc$CHIP == 1 & ukbc$MCHIP == 1,],"L-CHIP + M-CHIP")

## UKB 500
ukbc1<-data.frame(fread(paste(figDir,"data/UKB_500K.txt",sep="")))
tdf5<-getDF(ukbc1[ukbc1$CNV == 0,],"No mCA")
tdf6<-getDF(ukbc1[ukbc1$LCNV == 1 & ukbc1$MCNV == 0,],"L-mCA")
tdf7<-getDF(ukbc1[ukbc1$LCNV == 0 & ukbc1$MCNV == 1,],"M-mCA")
tdf8<-getDF(ukbc1[ukbc1$LCNV == 1 & ukbc1$MCNV == 1,],"A-mCA")
tdf9<-getDF(ukbc1[ukbc1$LCNV == 0 & ukbc1$MCNV == 0 & ukbc1$CNV == 1,],"U-mCA")

tdf<-rbind(tdf1,tdf2,tdf3,tdf4,tdf5,tdf6,tdf7,tdf8,tdf9)

write.table(tdf,file=paste(figDir,"scripts/Final_figs/ST_8_CH_characteristics.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
