/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk SplitVcfs -I vcfs/combined.vcf.gz --SNP_OUTPUT vcfs/combined_SNPs.vcf.gz --INDEL_OUTPUT vcfs/combined_INDELs.vcf.gz --STRICT false

echo "starting filtering SNPs"
bcftools filter -O u -s "QD"			-m + -e "QD < 2.0" vcfs/combined_SNPs.vcf.gz | \
	bcftools filter -O u -s "MQ"		-m + -e "MQ < 40.0" | \
	bcftools filter -O u -s "FS"		-m + -e "FS > 60.0" | \
	bcftools filter -O u -s "MQRank"	-m + -e "MQRankSum < -12.5" | \
	bcftools filter -O u -s "ReadPos"	-m + -e "ReadPosRankSum < -8.0" | \
	bcftools view -O z -f "PASS" > vcfs/combined_SNPs_filtered.vcf.gz
echo "finished filtering SNPs"

bcftools index vcfs/combined_SNPs_filtered.vcf.gz

echo "starting filtering INDELs"
bcftools filter -O u -s "QD"                    -m + -e "QD < 2.0" vcfs/combined_INDELs.vcf.gz | \
        bcftools filter -O u -s "FS"         -m + -e "FS > 200.0" | \
        bcftools filter -O u -s "ReadPos"    -m + -e "ReadPosRankSum < -20.0" | \
        bcftools view -O z -f "PASS" > vcfs/combined_INDELs_filtered.vcf.gz
echo "finished filtering INDELs"

bcftools index vcfs/combined_INDELs_filtered.vcf.gz

bcftools concat -O z -a vcfs/combined_INDELs_filtered.vcf.gz vcfs/combined_SNPs_filtered.vcf.gz > vcfs/combined_filtered.vcf.gz

bcftools +fill-tags vcfs/combined_filtered.vcf.gz -Oz -o vcfs/combined_filtered.VAF.vcf.gz -- -t FORMAT/VAF

/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk IndexFeatureFile -I vcfs/combined_filtered.vcf.gz


