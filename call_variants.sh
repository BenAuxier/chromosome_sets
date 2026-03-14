#!/bin/bash
set -euo pipefail

#The above command is called a shebang (or hashbang), and it's used at the very top of a script file to tell the sys>

REFERENCE=/lustre/BIF/nobackup/auxie001/one_house/genome/Botrytis_cinerea.ASM83294v1.dna_sm.toplevel.fa
#/mnt/LTR_userdata/auxie001/programs/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index $REFERENCE


#mkdir -p bams
# This calls the bwa-mem2 program and runs the index command on the file specified by the REFERENCE variable. The in>
#for i in B0s_1 B0s_2 B0s_3 B0s_4 B0s_5 B20s_1 B20s_2 B20s_3 B20s_4 B20s_5 B20s_6 B20s_7 B20s_8 B100s_1 B100s_2 B100s_3
#do
#    echo $i
#    r1=/lustre/BIF/nobackup/auxie001/one_house/data/X204SC25015728-Z01-F002/01.RawData/$i/B*_1.fq.gz
#    echo $i
#    r2=/lustre/BIF/nobackup/auxie001/one_house/data/X204SC25015728-Z01-F002/01.RawData/$i/B*_2.fq.gz
#    echo "Processing $i..."
#    /mnt/LTR_userdata/auxie001/programs/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 8 -R "@RG\tID:${i}\tSM:${i}\tLB:fumigatus" "$REFERENCE" $r1 $r2 | \
#    samtools view -b | \
#    samtools fixmate -@ 6 -m - - | \
#    samtools sort -@ 8 -m 3G - |
#    samtools markdup -@ 6 - bams/"${i}.sorted.bam"
#done

# Index the reference genome for samtools and GATK

#samtools faidx $REFERENCE
#/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk CreateSequenceDictionary -R $REFERENCE

# This script indexes sorted BAM files and calls variants for each sample using GATK HaplotypeCaller.
# It processes each .sorted.bam file in the specified directory, generates a BAM index (.bai), and outputs a per-sam>

mkdir -p vcfs
for bam in bams/*.sorted.bam
do
    sample_name=$(basename "$bam" .sorted.bam)
    # Index the BAM file
    samtools index $bam
    echo $sample_name
    # Run GATK HaplotypeCaller in GVCF mode
    /mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk HaplotypeCaller \
        -I $bam \
        -R $REFERENCE \
        -O vcfs/$sample_name\.gvcf.gz \
        -ERC GVCF \
        -ploidy 2
done

/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk CombineGVCFs \
  -R "$REFERENCE" \
  -O vcfs/combined.gvcf.gz \
  $(printf -- "-V %s " vcfs/*.gvcf.gz)


/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk GenotypeGVCFs -R $REFERENCE -O vcfs/combined.vcf.gz \
					-V vcfs/combined.gvcf.gz
