setwd("./")
library("ggplot2")
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(viridis)
library(ggpubr)



#########################################add isoform quantification results for samples from KakutaniLab and Sazelab
l1=read.table("./list_mutant/KakutaniLab.ibm.isoform.txt",header = T)
l2=read.table("./list_mutant/SazeLab.ibm.isoform.txt",header = T)
l3=read.table("./list_mutant/IBM2mutant.ibm1.isoform.txt",header = T)
l4=read.table("./list_mutant/xinjianLab.ibm.isoform.txt",header = T)
l5=read.table("./list_mutant/ZilbermanLab.ibm.isoform.txt",header = T)


c1=l1 %>%
  group_by(Sample) %>%
  summarise(TPM_mean=mean(geneTPM),TPM_sd=sd(geneTPM),long_mean=mean(long),
            long_sd=sd(long),short_mean=mean(short),short_sd=sd(short)) %>%
  pivot_longer(
    cols=-Sample,
    names_to = c("type",".value"),
    names_sep = "_"
  )

c2=l2 %>%
  group_by(Sample) %>%
  summarise(TPM_mean=mean(geneTPM),TPM_sd=sd(geneTPM),long_mean=mean(long),
            long_sd=sd(long),short_mean=mean(short),short_sd=sd(short)) %>%
  pivot_longer(
    cols=-Sample,
    names_to = c("type",".value"),
    names_sep = "_"
  )

c3=l3 %>%
  group_by(Sample) %>%
  summarise(TPM_mean=mean(geneTPM),TPM_sd=sd(geneTPM),long_mean=mean(long),
            long_sd=sd(long),short_mean=mean(short),short_sd=sd(short)) %>%
  pivot_longer(
    cols=-Sample,
    names_to = c("type",".value"),
    names_sep = "_"
  )

c4=l4 %>%
  group_by(Sample) %>%
  summarise(TPM_mean=mean(geneTPM),TPM_sd=sd(geneTPM),long_mean=mean(long),
            long_sd=sd(long),short_mean=mean(short),short_sd=sd(short)) %>%
  pivot_longer(
    cols=-Sample,
    names_to = c("type",".value"),
    names_sep = "_"
  )

c5=l5 %>%
  group_by(Sample) %>%
  summarise(TPM_mean=mean(geneTPM),TPM_sd=sd(geneTPM),long_mean=mean(long),
            long_sd=sd(long),short_mean=mean(short),short_sd=sd(short)) %>%
  pivot_longer(
    cols=-Sample,
    names_to = c("type",".value"),
    names_sep = "_"
  )


p1=ggplot(c1,aes(x=Sample,y=mean,fill=type,ymin=mean-sd, ymax=mean+sd))+
  geom_bar(stat="identity", width=.5, position = "dodge")+
  geom_errorbar( width=0.2, alpha=1, size=.7,
                 position = position_dodge(0.5),color="black")+
  scale_fill_manual(values = c("#E77577", "#0065A2", "#FFD872"),
                    labels=c("IBM1-L","IBM1-S","total"))+
  xlab("KakutaniLab")+
  ylab("IBM1 Gene expression (TPM)")+
  theme_classic(base_size = 15)
p1

p2=ggplot(c2,aes(x=Sample,y=mean,fill=type,ymin=mean-sd, ymax=mean+sd))+
  geom_bar(stat="identity", width=.5, position = "dodge")+
  geom_errorbar( width=0.2, alpha=1, size=.7,
                 position = position_dodge(0.5),color="black")+
  scale_fill_manual(values = c("#E77577", "#0065A2", "#FFD872"),
                    labels=c("IBM1-L","IBM1-S","total"))+
  xlab("Sazelab")+
  ylab("IBM1 Gene expression (TPM)")+
  theme_classic(base_size = 15)
p2

p3=ggplot(c3,aes(x=Sample,y=mean,fill=type,ymin=mean-sd, ymax=mean+sd))+
  geom_bar(stat="identity", width=.5, position = "dodge")+
  geom_errorbar( width=0.2, alpha=1, size=.7,
                 position = position_dodge(0.5),color="black")+
  scale_fill_manual(values = c("#E77577", "#0065A2", "#FFD872"),
                    labels=c("IBM1-L","IBM1-S","total"))+
  xlab("IBM2 mutant-Plant Epigenetics")+
  ylab("IBM1 Gene expression (TPM)")+
  theme_classic(base_size = 15)
p3

p4=ggplot(c4,aes(x=Sample,y=mean,fill=type,ymin=mean-sd, ymax=mean+sd))+
  geom_bar(stat="identity", width=.5, position = "dodge")+
  geom_errorbar( width=0.2, alpha=1, size=.7,
                 position = position_dodge(0.5),color="black")+
  scale_fill_manual(values = c("#E77577", "#0065A2", "#FFD872"),
                    labels=c("IBM1-L","IBM1-S","total"))+
  xlab("jianxi lab")+
  ylab("IBM1 Gene expression (TPM)")+
  theme_classic(base_size = 15)
p4

p5=ggplot(c5,aes(x=Sample,y=mean,fill=type,ymin=mean-sd, ymax=mean+sd))+
  geom_bar(stat="identity", width=.5, position = "dodge")+
  geom_errorbar( width=0.2, alpha=1, size=.7,
                 position = position_dodge(0.5),color="black")+
  scale_fill_manual(values = c("#E77577", "#0065A2", "#FFD872"),
                    labels=c("IBM1-L","IBM1-S","total"))+
  xlab("Zilberman lab")+
  ylab("IBM1 Gene expression (TPM)")+
  theme_classic(base_size = 15)
p5

pdf('lab5_isoformratio.pdf',
    width=12,
    height=6)
plot_grid(p1,p2,p3, p4, p5,nrow=2)
dev.off()


