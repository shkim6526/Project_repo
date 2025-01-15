#!/bin/bash

REFERENCE=/linux/database/REFERENCE/GRCh37/UCSC
BAM_DIR=/linux/zz_NIPT/zz_Exam/zz_data/zz_sam/zz_bam
#BAM_DIR2=/mnt/zz_NIPT/zz_Exam/zz_data/zz_sam/zz_bam
OUTDIR=/linux/zz_NIPT/zz_Exam/zz_data/zz_sam/zz_bam/dir_cnvkit
#BAM_FILES=$(ls ${BAM_DIR}/*.sorted.bam)

#docker run -it -v /mnt/d/linux/:/linux biocontainers/cnvkit:v0.9.5-3-deb_cv1 cnvkit access -o $OUTDIR/zz_access.bed $REFERENCE/GRCh37.fa

## CNVkit reference ##
# Alternatively, you can create a "flat" reference of neutral copy number (i.e. log2 0.0) for each probe from the target and antitarget interval files. This still computes the GC content of each region if the reference genome is given.

#docker run -it -v /mnt/d/linux/:/linux biocontainers/cnvkit:v0.9.5-3-deb_cv1 cnvkit reference -o ${BAM_DIR}/FlatReference.cnn -f ${REFERENCE}/GRCh37.fa -t ${OUTDIR}/zz_access.bed

## CNVkit target ##
# Prepare a BED file of baited regions for use with CNVkit.
#docker run -it -v /mnt/d/linux/:/linux biocontainers/cnvkit:v0.9.5-3-deb_cv1 cnvkit target ${OUTDIR}/zz_access.bed --annotate ${OUTDIR}/refFlat.txt --split -o {OUTDIR}/my_targets.bed

## CNVkit batch ##
# The [batch --method wgs] option uses the given reference genome's sequencing-accessible regions ("access" BED) as the "targets" - these will be calculated on the fly if not provided. No "antitarget" regions are used. Since the input does not contain useful per-target gene labels, a gene annotation database is required and used to label genes in the outputs

# ** Final object: Copy number reference profile (.cnn) **
# In addtion to the columns present in the “target” and “antitarget” .cnn files, the reference .cnn file has the columns:
# - GC content of the sequence region (gc)
# - RepeatMasker-masked proportion of the sequence region (rmask)
# - Statistical spread or dispersion (spread)
# The log2 coverage depth is the robust average of coverage depths, excluding extreme outliers, observed at the corresponding bin in each the sample .cnn files used to construct the reference. The spread is a similarly robust estimate of the standard deviation of normalized log2 coverages in the bin. The depth column is the robust average of absolute-scale coverage depths from the input .cnn files, but without any bias corrections.
# To manually review potentially problematic targets in the built reference, you can sort the file by the spread column; bins with higher values are the noisy ones.
# It is important to keep the copy number reference file consistent for the duration of a project, reusing the same reference for bias correction of all tumor samples in a cohort. If your library preparation protocol changes, it’s usually best to build a new reference file and use the new file to analyze the samples prepared under the new protocol.

#docker run -it -v /mnt/d/linux/:/linux biocontainers/cnvkit:v0.9.5-3-deb_cv1 cnvkit batch -m wgs --normal ${BAM_DIR}/SRR19224104.sorted.bam ${BAM_DIR}/SRR19224105.sorted.bam ${BAM_DIR}/SRR19224106.sorted.bam ${BAM_DIR}/SRR19224107.sorted.bam ${BAM_DIR}/SRR19224108.sorted.bam ${BAM_DIR}/SRR19224109.sorted.bam ${BAM_DIR}/SRR19224110.sorted.bam ${BAM_DIR}/SRR19224111.sorted.bam ${BAM_DIR}/SRR19224112.sorted.bam ${BAM_DIR}/SRR19224113.sorted.bam ${BAM_DIR}/SRR19224114.sorted.bam ${BAM_DIR}/SRR19224115.sorted.bam ${BAM_DIR}/SRR19224116.sorted.bam ${BAM_DIR}/SRR19224117.sorted.bam ${BAM_DIR}/SRR19224118.sorted.bam ${BAM_DIR}/SRR19224119.sorted.bam ${BAM_DIR}/SRR19224120.sorted.bam --targets ${OUTDIR}/zz_access.bed --fasta ${REFERENCE}/GRCh37.fa --output-reference ${OUTDIR}/OutReference.cnn --output-dir ${OUTDIR} --diagram --scatter

## CNVkit fix ##
# Combine the uncorrected target and antitarget coverage tables (.cnn) and correct for biases in regional coverage and GC content, according to the given reference. Output a table of copy number ratios (.cnr).

reference_cnn=${OUTDIR}/OutReference.cnn 
#for i in `ls /mnt/zz_NIPT/zz_Exam/zz_data/zz_sam/zz_bam/*.bam`
for i in `ls | grep _1.fastq`
do

echo $i

	sample=`echo $i | sed -e 's/_1.*//g'`
	target_cnn=${sample}".sorted.targetcoverage.cnn"
	antitarget_cnn=${sample}".sorted.antitargetcoverage.cnn"

## Normalization
#docker run -it -v /mnt/d/linux/:/linux biocontainers/cnvkit:v0.9.5-3-deb_cv1 cnvkit fix ${OUTDIR}/${target_cnn} ${OUTDIR}/${antitarget_cnn} ${reference_cnn} -o ${OUTDIR}/${sample}.cnr

# ** Final object: Bin-level log2 ratios (.cnr) **
# In addtion to the "chromosome", "start", "end", "gene", "log2" and "depth" columns present in .cnn files, the .cnr file includes each bin's proportional weight or reliability (weight).
# The weight value is dereived from several sources:
# - The size of the bin relative to the average bin size (for targets or antitargets, separately)
# - For a paired or pooled reference, the deviation of the reference log2 value from neutral coverage (i.e. distance from 0.0)
# - For a pooled reference, the inverse of the variance (i.e. square of "spread" in the reference) of normalized log2 coverage values seen among all normal samples at that bin.
# This calculated value is used to weight the bin log2 ratio values during segmentation. Also, when a genomic region is plotted with CNVkit's "scatter" command, the size of the plotted datapoints is proportional to each bin's weight - a relatively small point indicates a less reliable bin.

## Segment
#docker run -it -v /mnt/d/linux/:/linux biocontainers/cnvkit:v0.9.5-3-deb_cv1 cnvkit segment ${OUTDIR}/${sample}.cnr -m cbs -o ${OUTDIR}/${sample}.cns

# ** Final object: Segmented log2 rations (.cns) **
# In addition to the "chromosome", "start", "end", "gene", "log2" and "depth" columns present in .cnr files, the .cns file format has the additional column "probes" indicating the number of bins covered by the segment.
# The gene column concatenates the gene names of all the bins that the segment covers. The weight column sums the bin-level weights, and the depth and log2 is the weighted mean of the input bin-level values corresponding to the segment.

## Call
#docker run -it -v /mnt/d/linux/:/linux biocontainers/cnvkit:v0.9.5-3-deb_cv1 cnvkit call ${OUTDIR}/${sample}.cns -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o ${OUTDIR}/${sample}.call.cns

## Scatter plot
#docker run -it -v /mnt/d/linux/:/linux biocontainers/cnvkit:v0.9.5-3-deb_cv1 cnvkit scatter ${OUTDIR}/${sample}.cnr --y-min -2.0 --y-max 2.0  -t -s ${OUTDIR}/${sample}.call.cns -o ${OUTDIR}/${sample}.scatter.png
#for j in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
#do
#	docker run -it -v /mnt/d/linux/:/linux biocontainers/cnvkit:v0.9.5-3-deb_cv1 cnvkit scatter ${OUTDIR}/${sample}.cnr --y-min -2.0 --y-max 2.0 -c $j -s ${OUTDIR}/${sample}.call.cns -o ${OUTDIR}/${sample}.$j.scatter.png
#done

done


