#!/bin sh

#iterating though references
for f in *.fasta

do
#iterating through read datasets
for fq in *.fastq

do
#mapping reference datasets
tmap index -f $f
tmap mapall -f $f -r $fq -n 2 -v -Y -u -o 0 stage1 map4 >$f.$fq.TMAP.sam

#modifying sequence alignment
#samtools view -bSh -o $f.$fq.TMAP.bam $f.$fq.TMAP.sam
#java -Xmx4g -jar ~/bin/AddOrReplaceReadGroups.jar I= $f.$fq.TMAP.bam O= $f.$fq.TMAP.h.bam RGID= $fq RGLB="200bp" RGPU="316chip" RGPL="ion" RGSM=$f 
#java -Xmx4g -jar ~/bin/SortSam.jar I= $f.$fq.TMAP.h.bam O= $f.$fq.TMAP.h.sorted.bam SO=coordinate
#java -Xmx4g -jar ~/bin/MarkDuplicates.jar REMOVE_DUPLICATES=TRUE I=$f.$fq.TMAP.h.sorted.bam O=$f.$fq.TMAP.h.sorted.dedup.bam M=$f.$fq.TMAP.metric.txt
#samtools index $f.$fq.TMAP.h.sorted.dedup.bam
#java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -R $f -I $f.$fq.TMAP.h.sorted.dedup.bam -T RealignerTargetCreator -o $f.$fq.TMAP.h.sorted.dedup.bam.intervals
#java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -R $f -I $f.$fq.TMAP.h.sorted.dedup.bam -T IndelRealigner -targetIntervals $f.$fq.TMAP.h.sorted.dedup.bam.intervals -o  $f.$fq.TMAP.h.sorted.dedup.realinged.bam
 
#GATK variant calling
#samtools faidx $f
#java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -R $f -I $f.$fq.TMAP.h.sorted.dedup.realinged.bam -T UnifiedGenotyper -glm SNP -out_mode EMIT_ALL_SITES -o $f.$fq.TMAP.h.sorted.dedup.realinged.vcf

done
done