library(ggplot2)
library(vcfR)
library(cowplot)
data <- read.vcfR("/Users/ben/Library/CloudStorage/OneDrive-WageningenUniversity&Research/one_house/combined_filtered.vcf.gz")
fix <- as.data.frame(getFIX(data))
#extract Allele depths, and clean up to make ref vs alt
ad <- as.data.frame(extract.gt(data,element="AD"))
#710 multiallelic variants, which can be ignored
ad <- as.data.frame(ad[is.biallelic(data)==TRUE,])
fix <- as.data.frame(fix[is.biallelic(data)==TRUE,])

ad_alt <- as.data.frame(lapply(ad, function(x) as.numeric(sub(".+,","",x))))
ad_ref <- as.data.frame(lapply(ad, function(x) as.numeric(sub(",.+","",x))))


#add the fix information to the AD
ad <- cbind(fix,ad)




#one filter is to keep only those where in the 5 replicate no UV lines there are at most 1 read for teh alternate allele
#pure_starting <- ad_alt$B0s_1<2 & ad_alt$B0s_2<2 & ad_alt$B0s_3<2 & ad_alt$B0s_4<2 & ad_alt$B0s_5<2
#now filter the three data frames to keep only sites with pure ancestral lines
#pure_ad <- ad[pure_starting==T,]
#pure_ad_alt <- ad_alt[pure_starting==T,]
#pure_ad_ref <- ad_ref[pure_starting==T,]
#now we want to calculate ratio of the reference to alternate
ad_ratio <- ad_alt/(ad_alt+ad_ref)
ad_total <- ad_alt+ad_ref
#since some sites have zero ref and alt, this gives NA, which causes errors and we need to fix
ad_ratio[is.na(ad_ratio)] <- 0
#we want to figure out which variants are unique to a given sample, so we make a criterion that counts any ratio > 0.1 as "alternate", we want variants that only appear in a single sample
unique <- rowSums(ad_ratio > 0) == 1

unique_filtered_mutations <- rbind(
ad[ad_ratio$B20s_1 > 0.15 & ad_total$B20s_1 > 60 & unique,],
ad[ad_ratio$B20s_2 > 0.15 & ad_total$B20s_2 > 44 & unique,],
ad[ad_ratio$B20s_3 > 0.15 & ad_total$B20s_3 > 15 & unique,],
ad[ad_ratio$B20s_4 > 0.15 & ad_total$B20s_4 > 15 & unique,],
ad[ad_ratio$B20s_5 > 0.15 & ad_total$B20s_5 > 15 & unique,],
ad[ad_ratio$B20s_6 > 0.15 & ad_total$B20s_6 > 15 & unique,],
ad[ad_ratio$B20s_7 > 0.15 & ad_total$B20s_7 > 15 & unique,],
ad[ad_ratio$B20s_8 > 0.15 & ad_total$B20s_8 > 15 & unique,],
ad[ad_ratio$B100s_1 > 0.15 & ad_total$B100s_1 > 61 & unique,],
ad[ad_ratio$B100s_2 > 0.15 & ad_total$B100s_2 > 24 & unique,],
ad[ad_ratio$B100s_3 > 0.15 & ad_total$B100s_3 > 16 & unique,]
)

unique_filtered_mutations$type <- paste0(unique_filtered_mutations$REF,"_",unique_filtered_mutations$ALT)
table(unique_filtered_mutations$type)
sum(unique_filtered_mutations$type == "T_C")
sum(unique_filtered_mutations$type == "A_G")
sum(nchar(unique_filtered_mutations$type) == 3)

ratios_20s_variants <- c(
  ad_ratio$B20s_1[ad_ratio$B20s_1 > 0.15 & ad_total$B20s_1 > 60 & unique],
  ad_ratio$B20s_2[ad_ratio$B20s_2 > 0.15 & ad_total$B20s_2 > 44 & unique],
  ad_ratio$B20s_3[ad_ratio$B20s_3 > 0.15 & ad_total$B20s_3 > 44 & unique],
  ad_ratio$B20s_4[ad_ratio$B20s_4 > 0.15 & ad_total$B20s_4 > 17 & unique],
  ad_ratio$B20s_5[ad_ratio$B20s_5 > 0.15 & ad_total$B20s_5 > 30 & unique],
  ad_ratio$B20s_6[ad_ratio$B20s_6 > 0.15 & ad_total$B20s_6 > 30 & unique],
  ad_ratio$B20s_7[ad_ratio$B20s_7 > 0.15 & ad_total$B20s_7 > 61 & unique],
  ad_ratio$B20s_8[ad_ratio$B20s_8 > 0.15 & ad_total$B20s_8 > 68 & unique] 
)

ratios_100s_variants <- c(
  ad_ratio$B100s_1[ad_ratio$B100s_1 > 0.15 & ad_total$B100s_1 > 61 & unique],
  ad_ratio$B100s_2[ad_ratio$B100s_2 > 0.15 & ad_total$B100s_2 > 24 & unique],
  ad_ratio$B100s_3[ad_ratio$B100s_3 > 0.15 & ad_total$B100s_3 > 16 & unique]
)

plot_20s <- ggplot() + theme_classic(base_size=16) +
  geom_histogram(aes(x=ratios_20s_variants),bins=100) +
  scale_x_continuous(limits=c(0,1.02),expand=c(0.01,0.01)) +
  scale_y_continuous(limits=c(0,45),expand=c(0,0)) +
  geom_vline(xintercept=0.147,lty=2) +
  labs(x=expression(paste(italic("de novo "),"allele ratio in colonies from 20 s UV exposed spores")))

plot_100s <- ggplot() + theme_classic(base_size=16) +
  geom_histogram(aes(x=ratios_100s_variants),bins=100) +
  scale_x_continuous(limits=c(0,1.02),expand=c(0.01,0.01)) +
  scale_y_continuous(limits=c(0,45),expand=c(0,0)) +
  geom_vline(xintercept=0.147,lty=2) +
  labs(x=expression(paste(italic("de novo "),"allele ratio in colonies from 100 s UV exposed spores")))

plot_grid(plot_20s,plot_100s,nrow=2)


