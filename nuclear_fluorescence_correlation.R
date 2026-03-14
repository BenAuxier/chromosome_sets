library(ggplot2)
library(tidyr)
library(readxl)
library(scales)

setwd("/Users/ben/Library/CloudStorage/OneDrive-WageningenUniversity&Research/one_house/Quantification_of_acriflavin_staining_(excel file)")

botrytis_feb2 <- read.csv("2Feb_Mixture_of_micro_and_conidiospores_Bc_Bbig_frame.csv")
botrytis_jan26 <- read.csv("26Jan_Mixture_of_micro_and_conidiospores_Bc_Bbig_frame.csv")
botrytis <- rbind(botrytis_feb2,botrytis_jan26)
head(botrytis)
nrow(botrytis)

#we want to calculate the model fit without the microspores
no_micro_botrytis <- botrytis[botrytis$Comment != "microscpores",]
botrytis_nomicro_model <- lm(no_micro_botrytis$Total_intensity ~ no_micro_botrytis$Number_of_nuclei)
summary(botrytis_nomicro_model)

#we want to make a best fit line along part of the graph
xtest <- 1
xmin <- 2
xmax <- 12
ytest <- coef(botrytis_nomicro_model)[1] + coef(botrytis_nomicro_model)[2]*xtest
ymin <- coef(botrytis_nomicro_model)[1] + coef(botrytis_nomicro_model)[2]*xmin
ymax <- coef(botrytis_nomicro_model)[1] + coef(botrytis_nomicro_model)[2]*xmax

ggplot(botrytis, aes(x=Number_of_nuclei,y=Total_intensity)) +
  #geom_abline(intercept=coef(botrytis_nomicro_model)[1],slope=coef(botrytis_nomicro_model)[2],color="blue")+
  geom_segment(aes(x=xmin, xend=xmax, y=ymin, yend=ymax),lty=2)+
  geom_jitter(width=0.2,size=0.3,alpha=0.5) +
  theme_classic() +
  stat_summary(data=subset(botrytis,Number_of_nuclei != 1),fun=mean,geom="point",size=3,shape=1)+
  stat_summary(data=subset(botrytis,Number_of_nuclei != 1),fun.data=mean_se,geom="errorbar",width=0.2)+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12))+
  labs(x="Number of nuclei per spore",y=expression("Total fluorescence arbitrary units (" %*% 10^6 * ")")) +
  scale_y_continuous(labels=label_number(scale = 1e-6))+
  annotate("text",
           x = 11,
           y = 1e4,
           label = sprintf("italic(R)^2 == %.3f", summary(botrytis_nomicro_model)$r.squared),
           parse = TRUE,
           size = 4)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave("botrytis_model.svg",width=3.5,height=3.5)


neuro <- read.csv("3Feb_conidiospores_Nc_Bbig frame.csv")
head(neuro)

#made linear model
neuro_model <- lm(neuro$Total_intensity ~ neuro$Number_of_nuclei)
summary(neuro_model)

#we want to make a best fit line along part of the graph
nc_xmin <- 1
nc_xmax <- 9
nc_ymin <- coef(neuro_model)[1] + coef(neuro_model)[2]*nc_xmin
nc_ymax <- coef(neuro_model)[1] + coef(neuro_model)[2]*nc_xmax

ggplot(neuro, aes(x=Number_of_nuclei,y=Total_intensity)) +
  #geom_abline(intercept=coef(botrytis_nomicro_model)[1],slope=coef(botrytis_nomicro_model)[2],color="blue")+
  geom_segment(aes(x=nc_xmin, xend=nc_xmax, y=nc_ymin, yend=nc_ymax),lty=2) +
  geom_jitter(width=0.2,size=0.3,alpha=0.5) +
  theme_classic() +
  stat_summary(data=neuro,fun=mean,geom="point",size=3,shape=1)+
  stat_summary(data=neuro,fun.data=mean_se,geom="errorbar",width=0.2)+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9))+
  labs(x="Number of nuclei per spore",y=expression("Total fluorescence arbitrary units (" %*% 10^6 * ")")) +
  scale_y_continuous(labels=label_number(scale = 1e-6))+
  annotate("text",
           x = 8,
           y = 1e4,
           label = sprintf("italic(R)^2 == %.3f", summary(neuro_model)$r.squared),
           parse = TRUE,
           size = 4)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave("neuro_model.svg",width=3.5,height=3.5)
