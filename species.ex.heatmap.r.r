#setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/3_RA_SMAD9_Project/temp_scan1kb/")
setwd("C:/Users/zywlm/OneDrive - University of Georgia/02_my project/12_ibm1/v3.0/4_expression/")
library("ggplot2")
library(tidyverse)
library(dplyr)
library(cowplot)

##########################################workspace
#########read tables

d_exp=read.table("./34.species.methylgeneTPM.txt",header=T)
d_gbM=read.table("./34.species.gbMcount.txt",header=T)%>%
  mutate(ratio=gbMgene/allgene)
sp_order=read.table("./species.order.txt",header=T)

################main structure
####1. for target genes with multiple gene copies choose the one with maximum expression. for data point with multiple replicates average them.
##Consider using tradition way to get the table.
sp_str=split(d_exp,list(d_exp$species,d_exp$orthgene,d_exp$geneid),drop = T)
d_mean=data.frame()

### get the mean TPM for each geneid's TPM value, also consider distinguish 

for(i in seq_along(sp_str)){
  
  df <- data.frame(matrix(ncol = 4, nrow = 0))

  
  if(sp_str[[i]]$TPM[1]=='notmatch'){
    
    df=c(species=sp_str[[i]]$species[1],orthgene=sp_str[[i]]$orthgene[1],geneid=sp_str[[i]]$geneid[1],TPM=sp_str[[i]]$TPM[1])
  }
  else{
    tpmmean=mean(as.numeric(sp_str[[i]]$TPM))
    df=c(species=sp_str[[i]]$species[1],orthgene=sp_str[[i]]$orthgene[1],geneid=sp_str[[i]]$geneid[1],TPM=tpmmean)
    
  }
  d_mean=rbind(d_mean,df)
 # break
}
colnames(d_mean) <-c("species", "orthgene","geneid","TPM")
### get maximum expression from multiple gene copies

sp_str2=split(d_mean,list(d_mean$species,d_mean$orthgene),drop = T)
d_max=data.frame()

for(i in seq_along(sp_str2)){
  
  df <- data.frame(matrix(ncol = 4, nrow = 0))
  
  
  if(sp_str2[[i]]$TPM[1]=='notmatch'){
    
    df=c(species=sp_str2[[i]]$species[1],orthgene=sp_str2[[i]]$orthgene[1],geneid="NA",TPM="NA")
  }
  
  else{
    
    tpmmax=max(as.numeric(sp_str2[[i]]$TPM))
    index=which.max(as.numeric(sp_str2[[i]]$TPM))
    df=c(species=sp_str2[[i]]$species[1],orthgene=sp_str2[[i]]$orthgene[1],geneid=sp_str2[[i]]$geneid[index],TPM=tpmmax)
    #break  
  }

  d_max=rbind(d_max,df)
  #break
}
View(d_max)
colnames(d_max) <-c("species", "orthgene","geneid","TPM")

####2. make heatmap to express expression of target genes, and for the lost/unmatched genes, just show it by black color for example.

## 
d_max$TPM=as.numeric(d_max$TPM)
d_max$species=factor(d_max$species,levels = sp_order$orderedspecies) 
d_max$orthgene=factor(d_max$orthgene,levels=c("CMT2","ZMET1","CMT1","CMT3","IBM1","SUVH4"))

p1=ggplot(d_max, aes(x=orthgene,y=species,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 2.5,limits=c(0, 5))+
  theme_classic()
p1

####3. make scatter point for the gbM gene ratio

d_gbM$Species=factor(d_gbM$Species,levels = sp_order$orderedspecies)

p2= ggplot(d_gbM, aes(x=ratio,y=Species))+
  geom_point(aes(size=ratio),alpha=0.7,color="#3a53a4")+
  theme_classic(base_size = 15)+
  ylab("")+
  xlab("gbM genes ratio")+
  theme(legend.position = "none",axis.text.y=element_blank())
p2

####4. combine 2 figure together
pdf('34sp.ex.gbm.pdf',
    width=12,
    height=9)

plot_grid(p1, p2, nrow = 1,rel_widths = c(2, 1))

dev.off()

##############################including only result from Bra family
bra_sp=c("Arabidopsis_lyrata","Arabidopsis_thaliana","Eutrema_salsugineum","Boechera_stricta","Brassica_oleracea",
         "Brassica_rapa","Capsella_rubella","Thellungiella_parvula","Thlaspi_arvense")

d_max_bra=d_max[d_max$species %in% bra_sp,]%>%
  filter(orthgene!= 'ZMET1')

d_gbM_bra=d_gbM[d_gbM$Species %in% bra_sp,]

p1=ggplot(d_max_bra, aes(x=orthgene,y=species,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 2.5,limits=c(0, 5))+
  theme_classic()
p1

p2= ggplot(d_gbM_bra, aes(x=ratio,y=Species))+
  geom_point(aes(size=ratio),alpha=0.7,color="#3a53a4")+
  theme_classic(base_size = 15)+
  ylab("")+
  xlab("gbM genes ratio")+
  theme(legend.position = "none",axis.text.y=element_blank())
p2

pdf('Bra.9.ex.gbm.pdf',
    width=8,
    height=6)

plot_grid(p1, p2, nrow = 1,rel_widths = c(2, 1))

dev.off()