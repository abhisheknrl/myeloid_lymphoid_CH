

## Fig 2


library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)


library(survival)
library(survminer)


ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
hemeICDFile<-"meta/heme_malignancies_icd10.txt"


figDir<-""


#=====================================================================
##
ukb<-data.frame(fread(paste(figDir,"1_UKB_data.txt",sep="")))

#=====================================================================
##

hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")
myeloid_icdList<-as.character(hicd[hicd$Group == "Myeloid",]$ICD)
lymphoid_icdList<-as.character(hicd[hicd$Group == "Lymphoid",]$ICD)
unspe_icdList<-as.character(hicd[hicd$Group == "Unspecified",]$ICD)


## Diagnosis with dates
getDiag<-function(ukb,icdlist,columnNames,dateCols){
	K<-1
	tdf<-c()
	for (i in 1:length(columnNames)){
		CN<-columnNames[i]
		DCN<-dateCols[i]
		#print (CN)
		temp<-ukb[!is.na(ukb[,c(CN)]),]
		if (nrow(temp) == 0){
			next
		}
		temp<-temp[temp[,c(CN)] %in% icdlist,]
		if (nrow(temp)==0){
			next
		}
		## Get the date and eids
		tempdf<-data.frame(eid = temp$eid,
			ICD = temp[,c(CN)],
			Date = temp[,c(DCN)],
			DOAAC = temp$date_attending_assessment_center)
		if (K == 1){
			K<-2
			tdf<-tempdf
		}else{
			tdf<-rbind(tdf,tempdf)
		}
	}
	return (tdf)
}

## Cancer codes
getCancer<-function(ukb,icdlist,columnNames,dateCols){
	K<-1
	tdf<-c()
	for (i in 1:length(columnNames)){
		CN<-columnNames[i]
		DCN<-dateCols[i]
		#print (CN)
		temp<-ukb[!is.na(ukb[,c(CN)]),]
		if (nrow(temp) == 0){
			next
		}
		if (length(icdlist) > 0){
			temp<-temp[temp[,c(CN)] %in% icdlist,]
		}
		if (nrow(temp)==0){
			next
		}
		## Get the date and eids
		tempdf<-data.frame(eid = temp$eid,
			ICD = temp[,c(CN)],
			Date = temp[,c(DCN)],
			DOAAC = temp$date_attending_assessment_center)
		if (K == 1){
			K<-2
			tdf<-tempdf
		}else{
			tdf<-rbind(tdf,tempdf)
		}
	}
	return (tdf)
}


## Get unique samples
getUnique<-function(tdf){
	tdf$Diff<-as.numeric(difftime(tdf$Date, tdf$DOAAC,units=c("days")))
	## Finding the first diagnosis
	T<-table(tdf$eid)
	t1<-tdf[tdf$eid %in% names(T[T==1]),]
	t2<-tdf[tdf$eid %in% names(T[T>1]),]
	for (ID in names(T[T>1])){
		temp<-tdf[tdf$eid == ID,]
		temp<-temp[order(temp$Diff,decreasing=FALSE),]
		tnk<-data.frame(eid = ID,
			ICD = temp$ICD[1],
			Date = temp$Date[1],
			DOAAC = temp$DOAAC[1],
			Diff = temp$Diff[1])
		t1<-rbind(t1,tnk)
	}
	return (t1)
}


cancer10<-c("type_of_cancer_icd10",paste("type_of_cancer_icd10",1:16,sep="_"))
cancerDate10<-c("date_cancer_diagnosed",paste("date_cancer_diagnosed",1:16,sep="_"))

any10<-paste("f.41270.0.",0:212,sep="")
any9<-paste("f.41271.0.",0:46,sep="")
date10<-paste("f.41280.0.",0:212,sep="")
date9<-paste("f.41281.0.",0:46,sep="")


forcedCensorDate<-"2020-03-31"

###
## MY/LY
tdf1<-getDiag(ukb,c(myeloid_icdList,lymphoid_icdList),any10,date10)
tdf2<-getCancer(ukb,c(myeloid_icdList,lymphoid_icdList),cancer10,cancerDate10)
tempdf<-rbind(tdf1,tdf2)
tempdf$Diff2<-as.numeric(difftime(forcedCensorDate,tempdf$Date,units=c("days")))
tempdf<-tempdf[tempdf$Diff2 >= 0,]
tdf<-getUnique(tempdf[,colnames(tdf1)])

## For the unspecified types
# t<-getUnique(tempdf[!tempdf$ICD %in% c("C859","C851","C857"),colnames(tdf1)])
# t1<-tdf[!tdf$eid %in% t$eid,]
# tdf<-rbind(t,t1)


ICD<-as.character(tdf$ICD)
names(ICD)<-as.character(tdf$eid)
ukb$ICD<-ICD[as.character(ukb$eid)]

write.table(ukb[!is.na(ukb$ICD),],file=paste(figDir,"scripts/data/First_heme_diagnosis.txt",sep=""),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


#============================================================================
## ICD list


icdList<-as.character(hicd[hicd$Group %in% c("Lymphoid","Plasma_cell"),]$ICD)
CLL<-c("C911","C830")
PC<-c("C900","D472")
HL<-c("C81",paste("C81",0:9,sep=""))
FL<-c("C82",paste("C82",0:9,sep=""))
DLBCL<-c("C833")
NHL<-c("C85",paste("C85",0:9,sep=""))
WS<-c("C880")
CIRC<-c("C901","C910","C913","C914","C831","C841","C915","C916")
nhl<-icdList[!icdList %in% c(PC,CLL,HL,FL,DLBCL,NHL,WS,CIRC)]


mds_icdList<-c(paste("D46",0:9,sep=""),"C946")
aml_icdList<-paste("C92",0,sep="")
mpn_icdList<-c("D752","D45","C945","D471","D473")
OM_icdList<-hicd[hicd$Group == "Myeloid",]$ICD
OM_icdList<-OM_icdList[OM_icdList %in% c(mds_icdList,aml_icdList,mpn_icdList)]


df<-data.frame(eid = ukb[!is.na(ukb$ICD),]$eid,
	Malignancy = "Other_myeloid")


print (table(ukb$ICD))
print (table(ukb[ukb$CHIP == 1 | ukb$CNV == 1,]$ICD))

df$Malignancy<-as.character(df$Malignancy)
df[df$eid %in% ukb[ukb$ICD %in% CLL,]$eid,]$Malignancy<-"CLL"
df[df$eid %in% ukb[ukb$ICD %in% PC,]$eid,]$Malignancy<-"PC"
df[df$eid %in% ukb[ukb$ICD %in% FL,]$eid,]$Malignancy<-"FL"
df[df$eid %in% ukb[ukb$ICD %in% HL,]$eid,]$Malignancy<-"HL"
df[df$eid %in% ukb[ukb$ICD %in% DLBCL,]$eid,]$Malignancy<-"DLBCL"
df[df$eid %in% ukb[ukb$ICD %in% NHL,]$eid,]$Malignancy<-"NHL"
df[df$eid %in% ukb[ukb$ICD %in% WS,]$eid,]$Malignancy<-"Walden"
df[df$eid %in% ukb[ukb$ICD %in% CIRC,]$eid,]$Malignancy<-"CIRC"
df[df$eid %in% ukb[ukb$ICD %in% nhl,]$eid,]$Malignancy<-"Other_lymphoid"

df[df$eid %in% ukb[ukb$ICD %in% mds_icdList,]$eid,]$Malignancy<-"MDS"
df[df$eid %in% ukb[ukb$ICD %in% mpn_icdList,]$eid,]$Malignancy<-"MPN"
df[df$eid %in% ukb[ukb$ICD %in% aml_icdList,]$eid,]$Malignancy<-"AML"


T<-table(df[df$Malignancy %in% c("MDS","MPN","AML"),]$Malignancy)
T<-T[order(T,decreasing=TRUE)]
T<-T[T!=0]
T1<-table(df[!df$Malignancy %in% c("MDS","MPN","AML","Other_lymphoid","Other_myeloid","PC","CLL","CIRC","NHL"),]$Malignancy)
T1<-T1[order(T1,decreasing=TRUE)]
T1<-T1[T1!=0]


df$Malignancy<-factor(df$Malignancy,levels=c(names(T),"Other_myeloid","PC","CLL","CIRC","NHL",names(T1),"Other_lymphoid"))
rownames(df)<-paste("S",df$eid,sep="")
df<-df[order(df$Malignancy,decreasing=FALSE),]

#=====================================================================
## pheatmap

library(ComplexHeatmap)
library(circlize)


chip<-read.table(paste(figDir,"1_CHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
mychip<-read.table(paste(figDir,"1_MCHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
cnv<-read.table(paste(figDir,"data/UKB_500K_CNV.txt",sep=""),head=TRUE,sep="\t")
cnv<-cnv[cnv$CHR %in% 1:22,]


## canonical myeloid/lymphoid
CAN<-unique(c("gain12q","gain15q","gain17q","gain21q","gain22q13.32","gain2p","gain3q","gain8q","gain9q","del10p","del10q","del11q","del13q","del14q","del15q","del17p","del1p","del1q","del22q","del6q","del7q","del8p","tri12","tri18","tri19","gain1q","gain21q","gain9p","del12q","del20q","del5q","tri8"))
LG<-c("ATM","CIITA","ITPKB/NTRK1","KMT2C","MIR16-1","NCOR2","NOTCH1")
MG<-c("CBL","CTCF","EP300","JAK2","TP53","MPL/GNB1")

## CNV
cnv$VT<-paste(cnv$eid,cnv$CHR,cnv$Start,cnv$End,sep="_")
cnv1<-cnv[cnv$Canonical_CA %in% CAN,]
cnv2<-cnv[!cnv$VT %in% cnv1$VT & cnv$myGenes %in% c(MG,LG),]
cnv3<-cnv[!cnv$VT %in% c(cnv1$VT,cnv2$VT) & cnv$lyGenes %in% c(MG,LG),]
cnv4<-cnv[!cnv$VT %in% c(cnv1$VT,cnv2$VT,cnv3$VT),]

## CNV
cnv1$Hugo<-cnv1$Canonical_CA
cnv2$Hugo<-paste("LOH",cnv2$myGenes,sep="_")
cnv3$Hugo<-paste("LOH",cnv3$lyGenes,sep="_")
cnv4$Hugo<-rep("Other_mCA",nrow(cnv4))


cnv<-rbind(cnv1,cnv2,cnv3,cnv4)
mlCNV<-c("gain21q","LOH_MPL/GNB1","LOH_CBL","LOH_TP53","LOH_CTCF")
myCNV<-c("gain1q","gain9p","del12q","del20q","del5q","tri8","LOH_EP300","LOH_JAK2","LOH_TCL1A")

#myCNV<-c("gain9p","gain9q34","gain22q","gain1q","del5q","del14q23","LOH_TET2","LOH_JAK2","LOH_EP300","LOH_TCL1A")
lyCNV<-unique(as.character(cnv$Hugo))
lyCNV<-lyCNV[!lyCNV %in% c(myCNV,mlCNV,"Other_mCA")]
##
cnv1<-cnv
cnv<-cnv[cnv$eid %in% ukb$eid,]
## Make a list of Individuals first
eidList<-paste("S",ukb[!is.na(ukb$ICD),]$eid,sep="")
Genelist<-unique(c(as.character(mychip$Hugo_Symbol),myCNV,"Other_mCA",mlCNV,lyCNV,as.character(chip$Hugo_Symbol)))



#Create a matrix with the counts of each of these
m<-matrix(0,length(Genelist),length(eidList))
rownames(m)<-Genelist
colnames(m)<-eidList


## MyCHIP
for (EID in mychip$eid){
	if (!paste("S",EID,sep="") %in% eidList){
		next
	}
	T<-table(as.character(mychip[mychip$eid == EID,]$Hugo_Symbol))
	for (gene in names(T)){
		m[gene,paste("S",EID,sep="")]<-as.numeric(T[gene])
	}
}
## CHIP
for (EID in chip$eid){
	if (!paste("S",EID,sep="") %in% eidList){
		next
	}
	T<-table(as.character(chip[chip$eid == EID,]$Hugo_Symbol))
	for (gene in names(T)){
		m[gene,paste("S",EID,sep="")]<-as.numeric(T[gene])
	}
}
## CNV
for (EID in cnv$eid){
	if (!paste("S",EID,sep="") %in% eidList){
		next
	}
	T<-table(as.character(cnv[cnv$eid == EID,]$Hugo))
	for (gene in names(T)){
		## Myeloid
		if (gene %in% myCNV){
			m[gene,paste("S",EID,sep="")]<-as.numeric(T[gene])
		}
		## Lymphoid
		if (gene %in% lyCNV){
			m[gene,paste("S",EID,sep="")]<-as.numeric(T[gene])
		}
		## Ambiguous
		if (gene %in% mlCNV){
			m[gene,paste("S",EID,sep="")]<-as.numeric(T[gene])
		}
		## Unknowm
		if (gene %in% c("Other_mCA")){
			m[gene,paste("S",EID,sep="")]<-as.numeric(T[gene])
		}
	}
}

nCase<-rowSums(m)
nCase<-nCase[nCase >0]
m<-m[names(nCase),]
m[m>=2]<-2
hm1<-data.frame(m)

#=====================================================================
#=====================================================================
## ORder samples and genes

colOrder<-c("MPN","MDS","AML","Other_myeloid","PC","CLL","CIRC","NHL","DLBCL","FL","Walden","HL","Other_lymphoid")
getOrder<-function(hm1,sampdf,colOrder,myCNV,lyCNV,mlCNV,mchip,chip){
	mg<-c()
	lg<-c()
	mmca<-c()
	lmca<-c()
	amca<-c()
	SORTCOL<-c()
	for (CN in colOrder){
		print (CN)
		temp<-data.frame(hm1)[,paste("S",sampdf[sampdf$Malignancy == CN,]$eid,sep=""),]
		nCase<-rowSums(temp)
		nCase<-nCase[nCase>0]
		nCase<-nCase[order(nCase,decreasing=TRUE)]
		## Myeloid chip/mca
		tmc<-nCase[intersect(unique(as.character(mychip$Hugo_Symbol)),names(nCase))]
		tmc<-tmc[order(tmc,decreasing=TRUE)]
		tmc<-names(tmc)
		mg<-c(mg,tmc[!tmc %in% mg])
		##
		tmm<-nCase[intersect(myCNV,names(nCase))]
		tmm<-tmm[order(tmm,decreasing=TRUE)]
		tmm<-names(tmm)
		mmca<-c(mmca,tmm[!tmm %in% mmca])
		## Lymphoid chip/mca
		tlc<-nCase[intersect(unique(as.character(chip$Hugo_Symbol)),names(nCase))]
		tlc<-tlc[order(tlc,decreasing=TRUE)]
		tlc<-names(tlc)
		lg<-c(lg,tlc[!tlc %in% lg])
		##
		tlm<-nCase[intersect(lyCNV,names(nCase))]
		tlm<-tlm[order(tlm,decreasing=TRUE)]
		tlm<-names(tlm)
		lmca<-c(lmca,tlm[!tlm %in% lmca])
		## Ambiguous
		tam<-nCase[intersect(mlCNV,names(nCase))]
		tam<-tam[order(tam,decreasing=TRUE)]
		tam<-names(tam)
		amca<-c(amca,tam[!tam %in% amca])
		## Now order samples based on the genes mutated
		RN<-c(mg,mmca,lg,lmca,amca,"Other_mCA")
		RN<-RN[RN %in% names(nCase)]
		temp1<-data.frame(t(temp[RN,]))
		temp1$eid<-rownames(temp1)
		rownames(temp1)<-NULL
		setorderv(temp1, colnames(temp1)[-ncol(temp1)], rep(-1,ncol(temp1)-1))
		## SORT COLUMNS
		SORTCOL<-c(SORTCOL,temp1$eid)
	}
	SORTGENES<-c(mg,mmca,lg,lmca,amca,"Other_mCA")
	return (list(SORTGENES,SORTCOL))
}

GO<-getOrder(hm1,df,colOrder,myCNV,lyCNV,mlCNV,mchip,chip)

hm1<-hm1[GO[[1]],GO[[2]]]



RN<-rownames(hm1)
RN<-gsub("LOH_MPL.GNB1","LOH_MPL/GNB1",RN)
RN<-gsub("LOH_MIR16.1","LOH_MIR16-1",RN)
rownames(hm1)<-RN

#=====================================================================
## Add other data
#alc<-ukb$lymphocyte_count
#anc<-ukb$neutrophil_count
#rbc<-ukb$red_blood_cell_count
#plt<-ukb$platelet_count
#names(alc)<-names(anc)<-names(rbc)<-names(plt)<-paste("S",ukb$eid,sep="")
#df$ALC<-alc[rownames(df)]
#df$ANC<-anc[rownames(df)]
#df$RBC<-rbc[rownames(df)]
#df$PLT<-plt[rownames(df)]
#dfn<-colnames(df)
#dfn<-dfn[dfn!="eid"]
#df<-df[,dfn]
#df[!is.na(df$ALC) & df$ALC > 5,]$ALC<-5
#df[!is.na(df$PLT) & df$PLT > 500,]$PLT<-500

tempdf<-data.frame(Malignancy = df$Malignancy)
rownames(tempdf)<-rownames(df)

## N CH events among samples
T<-table(paste("S",c(mychip$eid,chip$eid, cnv$eid),sep=""))
nCH<-as.numeric(T[colnames(hm1)])


## Annotations
ha1 = HeatmapAnnotation(df = tempdf,
	bars=anno_barplot(nCH),
	col = list(Malignancy = c("MPN" = "red", "MDS" = "orange","AML" = "pink","Other_myeloid" = "gray31", "PC"="purple","CLL" = ctList[1], "CIRC" = "khaki", "NHL" = "brown", "DLBCL" = "green", "FL" = "violet", "Walden" = "skyblue", "HL" = "gold", "Other_lymphoid"="black")),
	na_col="white"
	)

#=====================================================================
## Row annotations



## df1
GeneList = rownames(hm1)

MG<-unique(as.character(mychip$Hugo_Symbol))
CG<-rep("M-CHIP",length(MG))
LG<-unique(as.character(chip$Hugo_Symbol))
CG<-c(CG,rep("L-CHIP",length(LG)))
MM<-unique(as.character(myCNV))
CG<-c(CG,rep("M-mCA",length(MM)))
LM<-unique(as.character(lyCNV))
CG<-c(CG,rep("L-mCA",length(LM)))
MLM<-unique(as.character(mlCNV))
CG<-c(CG,rep("A-mCA",length(MLM)))
OM<-c("Other_mCA")
CG<-c(CG,rep("U-mCA",length(OM)))
names(CG)<-c(MG,LG,MM,LM,MLM,OM)

df1<-data.frame(
	CH_categ = as.character(CG[GeneList])
)
rownames(df1)<-as.character(GeneList)

ha_row = rowAnnotation(df = df1,
	col = list(CH_categ = c("M-CHIP" = "red", "L-CHIP" = "blue","M-mCA" = "red","L-mCA" = "blue","A-mCA"="orange","U-mCA"="grey")))


RS<-as.character(df1$CH_categ)
names(RS)<-rownames(df1)

## Replacing 1s by something else for different color

t1<-t(hm1)
for (gene in colnames(t1)){
	if (gene %in% c(lyCNV,as.character(chip$Hugo_Symbol))){
		t1[,c(gene)]<-gsub(1,3,t1[,c(gene)])
	}
	## Ambiguous
	if (gene %in% mlCNV){
		t1[,c(gene)]<-gsub(1,4,t1[,c(gene)])
	}
	if (gene == "Other_mCA"){
		t1[,c(gene)]<-gsub(1,5,t1[,c(gene)])
	}
}

hm1<-t(t1)

## Heatmap
ht1<-Heatmap(hm1,
	col=c("white","red","black","blue","orange","grey"),
	split = factor(as.character(df1$CH_categ),levels=c("M-CHIP","M-mCA","L-CHIP","L-mCA","A-mCA","U-mCA")),
	cluster_rows = FALSE,
	cluster_columns = FALSE,
	row_names_side = "left",
	show_column_names = FALSE,
	top_annotation = ha1)

ht2<-ht1+ha_row


write.table(hm1,file = paste(figDir,"scripts/Final_figs/figure_df/F2_comut_matrix.txt",sep=""),
	col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")


#=====================================================================
#=====================================================================

pdf(paste(figDir,"scripts/Final_figs/Fig_ED2_comut.pdf",sep=""),width=14,height=14)
print (ht2)
dev.off()
