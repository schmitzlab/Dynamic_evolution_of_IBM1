#setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/3_RA_SMAD9_Project/temp_scan1kb/")
setwd("C:/Users/zywlm/OneDrive - University of Georgia/02_my project/12_ibm1/v5.0-refine/4_expression/")
library("ggplot2")
library(tidyverse)
library(dplyr)
library(cowplot)

# the script combined IBM1 gene expression with intron methylatin status information to check if there
# there is a significant different for expression between species with or without intronic methylation

# read tables
d_ex=read.table("./34.species.methylgeneTPM.txt",header = T)

d_mi=read.csv("./intron.status.txt",header = T)

# combine express and methylated intron information

ds_mi=d_mi %>%
  group_by(label) %>%
  count(status) %>%
  count(label) %>%
  separate(col = label, into = c("geneid", "species"), sep = "_(?=[^_]+_[^_]+$)", remove = TRUE, convert = FALSE) %>%
  mutate(n = if_else(n > 1, "with", "without"))


d_ex_IBM=d_ex[d_ex$orthgene=='IBM1'&d_ex$replicate==1,]

# output
write.csv(ds_mi,"./sp_intron_metylstatus.csv",quote=F, row.names  = F)
write.csv(d_ex_IBM,"./sp_IBM1_expression.csv",quote=F, row.names  = F)
# input
df_m=read.csv("./sp_intron_metylstatus.csv",header = T)
df_e=read.csv("./sp_IBM1_expression.csv",header=T)
# combine information
df_com=inner_join(df_m,df_e,by=c("geneid","species"))

# make analysis comparing expression between species with or without methylated introns
p1=ggplot(df_com,aes(x=n,y=TPM,fill=n))+
  geom_boxplot()+
  theme_classic(base_size = 20)+
  ylab("gene expression (TPM)")+
  xlab("")+
  theme(legend.position = "“none")

t.test(TPM~n,data=df_com)

# comarison between Capsella_rubella and Boechera_stricta

df_sub=df_com[df_com$species %in% c("Capsella_rubella","Boechera_stricta"),]

p2=ggplot(df_sub,aes(x=n,y=TPM+0.1,fill=n))+
  geom_bar(stat = "identity")+
  facet_grid(.~species)+
  theme_classic(base_size = 20)+
  ylab("gene expression (TPM)")+
  xlab("")+
  theme(legend.position = "“none")


pdf('exp.with.without.mintron.pdf',
    width=10,
    height=6)

plot_grid(p1, p2, nrow = 1,rel_widths = c(1, 1))

dev.off()