
geneIntervalFile<-"meta/CHIP_genes_hg38.txt"
inDir<-"sequencing_depth/"


library(data.table)
library(ggplot2)




## CHIP intervals
chip<-data.frame(fread(geneIntervalFile))


K<-1

for (ST in 10:60){
	print (ST)
inFile<-paste(inDir,"UKB_depth_samtools_out_",ST,".txt",sep="")
## samtools output file
f<-data.frame(fread(inFile))
colnames(f)<-c("CHR","POS",paste("S",ST,"_",3:ncol(f),sep=""))
## Loop over the genes
for (i in 1:nrow(chip)){
	CHR<-chip$Chr[i]
	Start<-chip$Start[i]
	End<-chip$End[i]
	Gene<-chip$Gene[i]
	## samtools
	t<-f[f$CHR == CHR & f$POS >= Start & f$POS <= End,]
	if (nrow(t) < 1){
		#print (Gene)
		next
	}
	CM<-colMeans(t[,3:ncol(t)])
	#print (Gene)
	#print (summary(t$POS))
	tempdf<-data.frame(Gene = Gene,
		CHR = CHR,
		Start = min(t$POS),
		End = max(t$POS),
		Sample = names(CM),
		nBases = nrow(t),
		meanDepth = CM)
	##
	if (K == 1){
		tdf<-tempdf
		K<-2
	}else{
		tdf<-rbind(tdf,tempdf)
	}
}
}

write.table(tdf,file=outFile,
	col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

outFile<-"scripts/data/UKB_depth_per_gene.txt"
write.table(tdf,file=outFile,
	col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


#======================================================================
## ST10

outFile2<-"scripts/data/UKB_depth_per_gene_median.txt"


mygeneFile<-"meta/CHIP_genelist.txt"
## gene list
gl<-data.frame(fread(mygeneFile,head=FALSE))

tdf$CHIP<-rep("L-CHIP",nrow(tdf))
tdf[tdf$Gene %in% gl[gl$V2 == "M-CHIP",]$V1,]$CHIP<-"M-CHIP"

tdf$CHIP<-factor(tdf$CHIP,levels=c("M-CHIP","L-CHIP"))
print (unique(tdf[!tdf$Gene %in% gl$V1,]$Gene))

tdf<-tdf[tdf$Gene %in% gl$V1,]


## Separate table for supplement
K<-1
for (g in unique(tdf$Gene)){
	t<-tdf[tdf$Gene == g,]
	q<-quantile(t$meanDepth,c(0.25,0.5,0.75))
	tempdf<-data.frame(Gene = g,
		nBases = t$nBases[1],
		CHIP = t$CHIP[1],
		Q1 = q[1],
		Q2 = q[2],
		Q3 = q[3])
	if (K == 1){
		mydf<-tempdf
		K<-2
	}else{
		mydf<-rbind(mydf,tempdf)
	}
}


write.table(mydf,file=outFile2,
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

