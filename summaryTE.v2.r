setwd("C:/Users/zywlm/OneDrive - University of Georgia/02_my project/12_ibm1/v5.0-refine//3_genetics/")
library("ggplot2")
library(tidyverse)
library(cowplot)
library(ggrepel)
#################################################################################################

mt=read.csv("./list/intron.methy.info.csv",header = T)
t=read.table("./list/34.addTE.anno.txt",header=T)
comb=full_join(mt,t,by = c("label" = "Geneid_species","start"="instart","end"="inend"))
tc=length(unique(paste(comb$label,comb$start,comb$end,sep="_")))###1total intron number
mc=length(mt[mt$status!='no',1])###2total mCHG methylated intron number
##3Total intron with TE
te=comb %>% filter(!is.na(motif))
tec=length(unique(paste(te$label,te$start,te$end,sep="_")))
##4total intron with TE and mCHG
mte=comb %>% filter(!is.na(motif))%>%filter(status !='no')
mtc=length(unique(paste(mte$label,mte$start,mte$end,sep="_")))

write.csv(comb, "./TE.intron.methyl.info.csv",quote = F, row.names = F)
####################################built a table for enrichment test
##in a order of 4 3 2 1

c1=c(mtc,mc,tec,tc)
c2=c(mc-mtc,mc,tc-tec,tc)
d=as.data.frame(rbind(c1,c2))
colnames(d)=c("k","n","m","N")
d$class=c("TE","no-TE")

a1=d$k
a2=d$n-a1
a3=d$m-a1
a4=d$N-d$n-a3
p=c()
for(i in 1:length(a1)){
  TeaTasting <-
    matrix(c(a1[i], a2[i], a3[i], a4[i]),
           nrow = 2,
           dimnames = list(Guess = c("Milk", "Tea"),
                           Truth = c("Milk", "Tea")))
  f=fisher.test(TeaTasting, alternative = "greater",conf.int = 0.95)
  
  p[i]=f$p.value
}
p

write.csv(d, "./TE.enrichment.table.csv",quote = F, row.names = F)

cp <- c("#F6D55C","grey")
d$kp=d$k/d$n
d$mp=d$m/d$N
d2=data.frame(class=d$class,kp=d$kp,mp=d$mp)
d3=reshape(d2,varying=c("kp","mp"),#combine wide format to long format by CG counts
           v.names="percentage",
           timevar="group",
           times=c("kp","mp"),
           direction = "long")
d3$group=factor(d3$group)
d3$class=factor(d3$class,levels=c("TE","no-TE"))
p1=ggplot(d3,aes(x=class,y=percentage,fill=group))+
  geom_bar(stat="identity",position=position_dodge())+
  scale_fill_manual(values = cp,name='',labels=c("mCHG methylated intron","Total intron"))+
  scale_x_discrete(labels=c("Intron-TE","no-TE"))+
  theme_classic(base_size=20) +
  ylab(" ")+
  xlab(" ")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))
p1
##################################output about TE type and subtype
type=as.data.frame(table(te$type))
type$pro=type$Freq/sum(type$Freq)
type=type[order(type$pro,decreasing = T),]
type$Var1=factor(type$Var1,levels = type$Var1)

df2 <- type %>% 
  mutate(csum = rev(cumsum(rev(pro))), 
         pos = pro/2 + lead(csum, 1),
         pos = if_else(is.na(pos), pro/2, pos))

p2=ggplot(type, aes(x = "" , y = pro, fill = fct_inorder(Var1))) +
  geom_bar(width = 1, stat = "identity")+
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set2") +
  geom_label(data = df2,
                   aes(y = pos, label = paste0(round(pro*100,digits=1),"%")),
                   size = 5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Group"))+
  theme_void(base_size=20)
p2


###try another way to show the repeat information
df=te %>%
  separate(col = label, into = c("geneid", "species"), sep = "_(?=[^_]+_[^_]+$)", remove = TRUE, convert = FALSE) %>%
  mutate(status = if_else(status != "no", "methyl", status))

df_sub=df[,c(1,2,17,21)]

df_summary <- df_sub %>%
  group_by(type) %>%
  summarise(num_species = n_distinct(species),
            num_species_methyl = n_distinct(species[status == "methyl"]))

df_summary2 <- df_sub %>%
  group_by(type) %>%
  summarise(num_gene = n_distinct(geneid),
            num_gene_methyl = n_distinct(geneid[status == "methyl"]))

df_long <- df_summary %>%
  pivot_longer(
    cols = c(num_species, num_species_methyl), # Columns to make long
    names_to = "measurement_type",            # New column for the original column names
    values_to = "count",                      # New column for the values
    names_prefix = "num_species_"             # Remove prefix from original column names for clarity
  )%>%
  mutate(
    source = case_when(
      measurement_type == "methyl" ~ "Methyl",
      TRUE ~ "Total"
    )
  )

df_long$type=factor(df_long$type,levels = c("Simple_repeat","DNA","LTR","LINE","RC","SINE","rRNA"))

p_dis=ggplot(df_long,aes(x=type,y=count,fill=source))+
  geom_bar(stat = "identity",position="dodge")+
  theme_classic(base_size=20) +
  xlab("")+
  ylab("the number of species")
  #theme(legend.position = "none")


pdf('TE.dis.species.pdf',
    width=8,
    height=6)
print(p_dis)

dev.off()
#########################################gene duplication distribution
mt$species=sub(".*\\_(\\w+\\_\\w+)$", "\\1", mt$label)
mt$gene=sub("(.*)\\_\\w+\\_\\w+$", "\\1", mt$label)
class=read.csv("./list/34.sp.classinfo.csv",header=T)
add_class=full_join(mt,class,by = c("species" = "name"))

dist=add_class %>%
  filter(class!='Basalmost_angiosperms')%>%
  distinct(gene,species,.keep_all = T)%>%
  group_by(species) %>%
  mutate(dup=n())%>%
  distinct(species,.keep_all = T)%>%
  group_by(dup) %>%
  mutate(count=n())%>%
  group_by(dup,class) %>%
  mutate(subcount=n())%>%
  mutate(freq =  100 *subcount/count) %>% 
  select(dup,class,count,subcount,freq)%>%
  distinct(dup,class,.keep_all = T)%>%
  data.frame()
cp=c("#F2789F","#36AE7C")
dist$class=factor(dist$class)
p3=ggplot(dist,aes(x=dup,y=subcount,fill=class))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = cp) +
  theme_classic(base_size=20) +
  xlab("the number of IBM1 gene copies")+
  ylab("the number of species")
p3
############################################variation of methylated intron pattern between gene copies.
dup=add_class %>%
  distinct(gene,species,.keep_all = T)%>%
  group_by(species) %>%
  mutate(dup=n())%>%
  distinct(species,dup)%>%
  filter(dup>=2)%>%
  data.frame()
View(dup)

dup_class=inner_join(dup,add_class,by = c("species"))
sp_g=split(dup_class, f = ~ dup_class$species+ dup_class$gene,drop = T)
mi_ls=list() ##record whether or not a gene got methylated intron or not

for(i in seq_along(sp_g)){
  index="no"
  for(j in seq_along(sp_g[[i]]$status)){
    if(sp_g[[i]][j,]$status!="no"){###  
      index="yes"
    }
  }
  df=data.frame(species=unique(sp_g[[i]]$species),gene=unique(sp_g[[i]]$gene),mintron=index)
  mi_ls=append(mi_ls,list(df)) 
  #df=rbind()
}
d_mi=do.call(rbind,mi_ls)
####################################variation from species perspective
dm2=split(d_mi, f = ~ d_mi$species,drop = T)
ms_ls=list() ##record whether duplicated IBM1 genes from a species have different methylation pattern on introns.
for(i in seq_along(dm2)){
  index="no"
  for(j in seq_along(dm2[[i]]$mintron)){
    if(length(levels(as.factor(dm2[[i]]$mintron)))==1){###  
      if(levels(as.factor(dm2[[i]]$mintron))=='yes'){
        index="all"
      }
    }
    else {
      index="var"
    }
  }
  df=data.frame(species=unique(dm2[[i]]$species),var=index)
  ms_ls=append(ms_ls,list(df)) 
  #df=rbind()
}
ms=do.call(rbind,ms_ls)
dvar=as.data.frame(table(ms$var))

p4=ggplot(dvar,aes(x=Var1,y=Freq))+
  geom_bar(stat="identity",fill="#1a53ff")+
  #scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size=20) +
  xlab("")+
  ylab("the number of species")
p4



pdf('TE.summary.pdf',
    width=12,
    height=9)
plot_grid(p1, p2,p3,p4, labels = c('A', 'B','C', 'D'),ncol = 2)

dev.off()