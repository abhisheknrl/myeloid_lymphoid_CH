

library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)
library(gg3D)


library(survival)
library(survminer)


ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"


## list of samples thtat carry somatic variants and mCAs at the same locus
m<-c("##")
l<-c("###")

#=====================================================================
## Reload data

chip<-read.table(paste(figDir,"1_CHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
mychip<-read.table(paste(figDir,"1_MCHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
cnv<-read.table(paste(figDir,"1_CNV_UKB.txt",sep=""),head=TRUE,sep="\t")

cnv[is.na(cnv$CELL_FRAC),]$CELL_FRAC<-0.001

cnv<-cnv[cnv$CHR %in% 1:22,]


#=====================================================================
## Get max VAF for each sample

Tot<-length(unique(c(mychip$eid,chip$eid,cnv$eid)))
print (paste("Total: ",Tot,sep=""))


Ic<-unique(intersect(chip$eid,cnv$eid))
Im<-unique(intersect(mychip$eid,cnv$eid))
Icm<-unique(intersect(chip$eid,mychip$eid))

I<-unique(c(Ic,Im,Icm))
print (paste("Multiple events:",length(I),sep=""))

Tot1<-table(c(mychip$eid,chip$eid,cnv$eid))
Tot1<-names(Tot1[Tot1 >1])
Tot1<-Tot1[!Tot1 %in% I]
print (paste("Multiple events same group:",length(Tot1),sep=""))


for (EID in I){
	if (EID %in% chip$eid){
		tc<-max(chip[chip$eid == EID,]$VAF)
	}else{
		tc<-0
	}
	if (EID %in% mychip$eid){
		tm<-max(mychip[mychip$eid == EID,]$VAF)
	}else{
		tm<-0
	}
	if (EID %in% cnv$eid){
		tv<-max(cnv[cnv$eid == EID,]$CELL_FRAC)
	}else{
		tv<-0
	}
	## temp
	temp<-data.frame(eid = EID,
		CHIP = tc,
		MCHIP = tm,
		CNV = tv)
	if (EID == I[1]){
		tdf<-temp
	}else{
		tdf<-rbind(tdf,temp)
	}
}

tdf$L<-rep(1,nrow(tdf))
tdf[tdf$eid %in% l,]$L<-1.3
tdf$M<-rep(1,nrow(tdf))
tdf[tdf$eid %in% m,]$M<-1.3


## Color the plots
tdf$COL<-rep("black",nrow(tdf))
tdf[tdf$CHIP == 0,]$COL<-ctList[1]
tdf[tdf$CNV == 0,]$COL<-ctList[3]
tdf[tdf$MCHIP == 0,]$COL<-ctList[2]
tdf$COL<-factor(tdf$COL,levels=c("black",ctList))

	

temp<-tdf[tdf$CHIP != 0 & tdf$CNV !=0,]
CM<-cor.test(temp[temp$CNV >= 0.1,]$CHIP,temp[temp$CNV >= 0.1,]$CNV)
p3<-ggplot(temp,aes(CHIP,CNV,size=L))+
	geom_point(color="black")+
	theme_linedraw()+
	theme(panel.grid=element_blank(),
		legend.position="none")+
	xlab("L-CHIP VAF")+
	ylab("mCA cell fraction")+
	scale_size(range=c(1,3))



temp<-tdf[tdf$MCHIP != 0 & tdf$CNV !=0,]
CM<-cor.test(temp[temp$CNV >= 0.1,]$MCHIP,temp[temp$CNV >= 0.1,]$CNV)
p4<-ggplot(temp,aes(MCHIP,CNV,size=M))+
	geom_point(color="black")+
	theme_linedraw()+
	theme(panel.grid=element_blank(),
		legend.position="none")+
	xlab("M-CHIP VAF")+
	ylab("mCA cell fraction")+
	scale_size(range=c(1,3))

P1<-ggarrange(p4,p3,ncol=2,nrow=1)


setEPS()
postscript(paste(figDir,"scripts/Final_figs/Fig_ED5_cd_VAF_CF.eps",sep=""),width=7,height=3)
print (P1)
dev.off()

pdf(paste(figDir,"scripts/Final_figs/Fig_ED5_cd_VAF_CF.pdf",sep=""),width=7,height=3)
print (P1)
dev.off()

