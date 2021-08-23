

library(stringr)
library(ggplot2)
library(data.table)
library(ggpubr)
library(colortools)

ctList<-splitComp("#0000CD")
ctList[1:2]<-c("blue","red")


setwd("")
figDir<-""
hemeICDFile<-"meta/heme_malignancies_icd10.txt"

##===================================================================

##ICD10codes
hicd<-read.delim(hemeICDFile,head=TRUE,sep="\t")

#=====================================================================
## Reload data

ukb<-data.frame(fread(paste(figDir,"1_UKB_data.txt",sep="")))

chip<-read.delim(paste(figDir,"1_CHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
mychip<-read.delim(paste(figDir,"1_MCHIP_UKB.txt",sep=""),head=TRUE,sep="\t")
cnv<-data.frame(fread(paste(figDir,"data/UKB_500K_CNV.txt",sep="")))


## 
getFig<-function(temp1,mycol){
	temp1$VC<-as.character(temp1$Variant_Classification)
	temp1[temp1$VC == "MISSENSE",]$VC<-"Missense"
	temp1[temp1$VC == "NONSENSE",]$VC<-"Nonsense"
	temp1[temp1$VC %in% c("FRAME_SHIFT_DEL","FRAME_SHIFT_INS"),]$VC<-"Frame shift"
	temp1[temp1$VC == "SPLICE_SITE",]$VC<-"Splice site"
	if ("IN_FRAME_DEL" %in% temp1$VC | "IN_FRAME_INS" %in% temp1$VC){
		temp1[temp1$VC %in% c("IN_FRAME_DEL","IN_FRAME_INS"),]$VC<-"In frame"
	}

	#==============================
	## Variant type
	T<-table(temp1$Variant_Type)
	tempdf<-data.frame(Type = names(T),
		Count = as.numeric(T))
	tempdf$Type <- gsub("SNP","SNV",as.character(tempdf$Type))
	tempdf$Type <- gsub("DNP","DNV",as.character(tempdf$Type))
	tempdf<-tempdf[order(tempdf$Count,decreasing=TRUE),]
	tempdf$Type<-factor(tempdf$Type,levels=tempdf$Type)
	sf2a<-ggplot(tempdf,aes(y=Count,x=Type))+
		geom_bar(stat="identity",fill=mycol)+
		theme_classic()+
		xlab("Variant type")+
		ylab("Number of variants")+
		theme(legend.position="none",
			axis.text.y = element_text(color="black"),
			axis.text.x=element_text(angle=45,hjust=1,color="black"),
			axis.ticks = element_line(color="black"))
	#==============================
	### Substitutions
	nuc<-c("A","C","G","T")
	names(nuc)<-c("T","G","C","A")
	temp<-temp1
	temp$ref1<-temp$ref
	temp$alt1<-temp$alt
	temp<-temp[temp$ref %in% c("A","C","G","T") & temp$alt %in% c("A","C","G","T"),]
	t1<-temp[temp$ref %in% c("G","A"),]
	t2<-temp[!temp$ref %in% c("G","A"),]
	t1$ref<-nuc[as.character(t1$ref)]
	t1$alt<-nuc[as.character(t1$alt)]
	temp<-rbind(t1,t2)
	temp$Subs<-paste(temp$ref,temp$alt,sep=" > ")
	T<-table(temp$Subs)
	tempdf<-data.frame(Type = names(T),
		Count = as.numeric(T))
	tempdf$Prop<-round(tempdf$Count/sum(tempdf$Count),3)
	sf2b<-ggplot(tempdf,aes(y=Prop,x=Type))+
		geom_bar(stat="identity",fill=mycol)+
		theme_classic()+
		xlab("Nucleotide substitutions")+
		ylab("Number of variants")+
		theme(legend.position="none",
			axis.text.y = element_text(color="black"),
			axis.text.x=element_text(angle=45,hjust=1,color="black"),
			axis.ticks = element_line(color="black"))
	#==============================
	## Variant classification
	COL<-c("green3","brown1","dodgerblue3","darkorange","seagreen3")
	names(COL)<-c("Missense","Nonsense","Frame shift","Splice site","In frame")
	LEV<-c("In frame","Splice site","Frame shift","Nonsense","Missense")
	T<-table(temp1$VC)
	tempdf<-data.frame(Variant_class = names(T),
		Count = as.numeric(T))
	tempdf<-tempdf[order(tempdf$Count,decreasing=TRUE),]
	tempdf$Color<-COL[tempdf$Variant_class]
	tempdf<-tempdf[order(tempdf$Count,decreasing=TRUE),]
	tempdf$Variant_class<-factor(tempdf$Variant_class,levels=LEV)
	sf2c<-ggplot(tempdf,aes(y=Count,x=Variant_class,label=Count))+
		geom_bar(stat="identity", fill=mycol)+
		theme_classic()+
		xlab("Variant classification")+
		ylab("Number of variants")+
		#scale_fill_manual(values=COL[LEV])+
		theme(legend.position="none",
			axis.text.y = element_text(color="black"),
			axis.text.x=element_text(angle=45,hjust=1,color="black"),
			axis.ticks = element_line(color="black"))
	#==============================
	## genes
	T<-table(as.character(temp1$Hugo_Symbol))
	T<-T[order(T,decreasing=TRUE)]
	K<-1
	for (gene in names(T)[1:25]){
		temp<-temp1[temp1$Hugo_Symbol == gene,]
		##
		nmis<-nrow(temp[temp$VC == "Missense",])
		nnon<-nrow(temp[temp$VC == "Nonsense",])
		nfs<-nrow(temp[temp$VC == "Frame shift",])
		nss<-nrow(temp[temp$VC == "Splice site",])
		nif<-nrow(temp[temp$Variant_Classification %in% c("IN_FRAME_DEL","IN_FRAME_INS"),])
		tempdf<-data.frame(Gene = gene,
			Count = c(nmis,nnon,nfs,nss,nif),
			Group= c("Missense","Nonsense","Frame shift","Splice site","In frame"),
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
	tdf$Group<-factor(tdf$Group,levels=c("Missense","Nonsense","Frame shift","Splice site","In frame"))
	#sf2d<-ggplot(tdf,aes(y=Count,x=Gene,fill=Group))+
	sf2d<-ggplot(tdf,aes(y=Count,x=Gene))+
		geom_bar(stat="identity",color=mycol,fill=mycol)+
		theme_classic()+
		xlab("Gene")+
		ylab("Number of variants")+
		#geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25)+
		#ylim(0,35)+
		scale_fill_manual(values=COL)+
		theme(legend.position="none")+
		#coord_flip()+
		theme(axis.text.x=element_text(color="black",angle=45,hjust=1))
	return (list(sf2a,sf2b,sf2c,sf2d))
}



#=====================================================================
## VAF distribution

p51<-ggplot(chip,aes(VAF))+
	geom_histogram(bins=25,fill=ctList[1])+
	theme_linedraw()+
	theme(panel.grid=element_blank())+
	ylab("Number of variants")+
	xlab("Variant allele fraction")+
	xlim(0,1)
p52<-ggplot(mychip,aes(VAF))+
	geom_histogram(bins=25,fill=ctList[2])+
	theme_linedraw()+
	theme(panel.grid=element_blank())+
	ylab("Number of variants")+
	xlab("Variant allele fraction")+
	xlim(0,1)



## N mutation carriers plot
tdf<-data.frame(eid = chip$eid,
	VT = chip$VT)
tdf<-unique(tdf)

T<-table(tdf$eid)
tempdf<-data.frame(N = c(1,">1"),
	Count = c(length(T[T==1]),length(T[T>1])))

tempdf$N<-factor(tempdf$N,levels=c(1,">1"))

p61<-ggplot(tempdf,aes(x=N,y=Count,label=Count))+
	geom_bar(color=ctList[1],fill=ctList[1],stat="identity")+
	theme_classic()+
	xlab("Number of variants")+
	ylab("Number of samples")+
	ylim(0,700)
	#geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25)


## Myeloid CHIP
tdf<-data.frame(eid = mychip$eid,
	VT = mychip$VT)
tdf<-unique(tdf)

T<-table(tdf$eid)
tempdf<-data.frame(N = c(1,">1"),
	Count = c(length(T[T==1]),length(T[T>1])))

tempdf$N<-factor(tempdf$N,levels=c(1,">1"))

p62<-ggplot(tempdf,aes(x=N,y=Count,label=Count))+
	geom_bar(color=ctList[2],fill=ctList[2],stat="identity")+
	theme_classic()+
	xlab("Number of variants")+
	ylab("Number of samples")+
	ylim(0,3000)
	#geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25)


sfl<-getFig(chip,ctList[1])
sf<-getFig(mychip,ctList[2])


F1<-ggarrange(ggarrange(ggarrange(sfl[[1]],sf[[1]],ncol=3,nrow=1,widths=c(1,1,0.1)),
		ggarrange(sfl[[2]],sf[[2]],ncol=2,nrow=1),
		ncol=2,nrow=2,heights=c(1,0.1)),
	ggarrange(ggarrange(sfl[[3]],sf[[3]],ncol=3,nrow=1,widths=c(1,1,0.1)),
		ggarrange(p61,p62,ncol=2,nrow=1),
		ncol=2,nrow=2,widths=c(1.2,1),heights=c(1,0.1)),
	ggarrange(p51,p52,ncol=2,nrow=2,heights=c(1,0.1)),
	nrow=3,ncol=1)



#===========================================================================
## Individuals with both M-CHIP and L-CHIP

I<-intersect(as.character(mychip$eid),as.character(chip$eid))
mukb<-ukb[ukb$MCHIP == 1,]
lukb<-ukb[ukb$CHIP == 1,]

K<-1
for (S in I){
	lvaf<-unique(chip[chip$eid == S,c("eid","VT","VAF"),])$VAF
	mvaf<-unique(mychip[mychip$eid == S,c("eid","VT","VAF"),])$VAF
	tempdf<-data.frame(Sample = S,
		My_VAF = max(mvaf),
		Ly_VAF = max(lvaf))
	if (K == 1){
		tdf<-tempdf
		K<-2
	}else{
		tdf<-rbind(tdf,tempdf)
	}
}


mod<-lm(Ly_VAF ~ My_VAF,data=tdf)

sf1a<-ggplot(tdf,aes(My_VAF,Ly_VAF))+
	geom_point(size=1)+
	geom_smooth(method="lm",formula = y ~ x, se=FALSE,color="black",size=0.8)+
	theme_linedraw()+
	theme(panel.grid=element_blank())+
	xlab("M-CHIP VAF")+
	ylab("L-CHIP VAF")+
	xlim(c(0,0.45))+
	ylim(c(0,0.45))



## Genes overlapping in myeloid and lymphoid chip
T1<-table(as.character(mychip[mychip$eid %in% I,]$Hugo_Symbol))
T2<-table(as.character(chip[chip$eid %in% I,]$Hugo_Symbol))
tdf1<-data.frame(Gene = names(T1),
	Count = as.numeric(T1))
tdf2<-data.frame(Gene = names(T2),
	Count = as.numeric(T2))
tdf1<-tdf1[order(tdf1$Count,decreasing=TRUE),]
tdf1$Gene<-factor(tdf1$Gene,levels=tdf1$Gene)
tdf2<-tdf2[order(tdf2$Count,decreasing=TRUE),]
tdf2$Gene<-factor(tdf2$Gene,levels=tdf2$Gene)


sf1b<-ggplot(tdf1[tdf1$Count > 1,],aes(x=Gene,y=Count,label=Count))+
	geom_bar(color=ctList[2],fill=ctList[2],stat="identity")+
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.text.x = element_text(color="black",angle=45,hjust=1))+
	xlab("M-CHIP genes")+
	ylab("Number of variants")+
	ylim(0,55)


sf1c<-ggplot(tdf2[tdf2$Count > 1,],aes(x=Gene,y=Count,label=Count))+
	geom_bar(color=ctList[1],fill=ctList[1],stat="identity")+
	theme_classic()+
	theme(panel.grid=element_blank(),
		axis.text.x = element_text(color="black",angle=45,hjust=1))+
	xlab("L-CHIP genes")+
	ylab("Number of variants")+
	ylim(0,7)


F2<-ggarrange(ggarrange(sf1a,ncol=2,widths=c(1,0.1)),
	ggarrange(sf1c,sf1b,nrow=2,ncol=1),ncol=2,nrow=1)

F<-ggarrange(F1,F2,nrow=2,ncol=1,heights=c(1.5,1))


setEPS()
postscript(paste(figDir,"scripts/Final_figs/Fig_S1_CHIP_characteristics.eps",sep=""),width=10,height=12)
print (F)
dev.off()

pdf(paste(figDir,"scripts/Final_figs/Fig_S1_CHIP_characteristics.pdf",sep=""),width=10,height=12)
print (F)
dev.off()

