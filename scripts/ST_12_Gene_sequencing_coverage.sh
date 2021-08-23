

wd=$1


## Coverage at the target intervals

chipInterval=${wd}/CHIP_genes_hg38_intervals.txt
ukbInterval=${wd}/../xgen_plus_spikein.b38.bed
wesInterval=${wd}/CHIP_genes_UKB_WES_interval.txt
bamFile_list=${wd}/UKB_depth_bam_list.txt

## Intersect between the CHIP genes and the UKB target intervals
bedtools intersect -a ${chipInterval} -b ${ukbInterval} > ${wesInterval}



for i in $(seq 10 60)
do
echo ${i}
## Make list of bam files
ls ${wd}/bams/${i}*.bam > ${bamFile_list}
## Samtools depth at these sites
outFile=${wd}/depth/UKB_depth_samtools_out_${i}.txt
samtools depth -a -b ${wesInterval} -f ${bamFile_list} -q 30 -l 30 > ${outFile} 
done



### Rscript to compute per gene average depth

Rscript scripts/ST_12_UKB_depth_seq_genes.r

