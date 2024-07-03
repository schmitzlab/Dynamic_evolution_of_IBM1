setwd("C:/Users/zywlm/OneDrive - University of Georgia/02_my project/12_ibm1/v4.0/2_intron_methy//")
library("ggplot2")
library(grid)
library(scales)
library(tidyverse)
library(ggtree)
library(tidytree)
library(treeio)
library(patchwork)
library(cowplot)
library(ggpubr)
#################################################################################################
#########Section 1
#####Building IBM1 gene tree with methylated gene structure information
#####v2.0, adding indicator information in the intron methylation plot.
################################################################################################

#####################################read tree
tr <- read.tree("./data/convert.besttree.34.100.tre")
##########read newick format tree with bootstrap value at branches
#bstre=read.raxml("./data/convert.besttree.34.addbsbranch.tre")
#x=as_tibble(bstre)
##########set Amborella as outlier
rtre=root(tr,outgroup = "evm_27.model.AmTr_v1.0_scaffold00072.47_Amborella_trichopoda",edgelabel = TRUE)
###extract tip labels on the tree
tips=as_tibble(rtre)%>%
  filter(!is.na(label)&!is.na(branch.length)) %>%
  pull(label)
###########read species class and family information
info=read.csv("./data/34.sp.classinfo.csv",header = T)
#######################################find overlap between tree and information table
sp=c()
t=str_split(tips,"_")
for (i in 1:length(t)){
  pre=t[[i]][length(t[[i]])-1]
  last=t[[i]][length(t[[i]])]
  string=paste0(pre,"_",last)
  sp[i]=string
}
###combined information
sp34=data.frame(label=tips,name=sp)
inset=inner_join(sp34,info,by="name")
rtre_tb=as_tibble(rtre)
rtre_info=full_join(rtre_tb,inset,by="label")
rtre_final=as.treedata(rtre_info)
#View(rtre_info)
################################drop some tips
to_drop=c("Manes.11G102400.2_Manihot_esculenta","Lj0g3v0226069.1_Lotus_japonicus","Lj5g3v2258440.1_Lotus_japonicus",
          "Glyma.10G284500.1_Glycine_max","Glyma.20G104900.3_Glycine_max","Glyma.19G068800.2_Glycine_max","Glyma.19G064000.1_Glycine_max",
          "Cucsa.067550.1_Cucumis_sativus","FvH4_5g29990.t3_Fragaria_vesca","FvH4_7g06860.t1_Fragaria_vesca","Gorai.008G148100.1_Gossypium_raimondii",
          "Thecc1EG033453t2_Theobroma_cacao","Migut.A00947.1_Mimulus_guttatus","Migut.K01225.1_Mimulus_guttatus","Solyc04g049140.3.1_Solanum_lycopersicum",
          "Thecc1EG033531t1_Theobroma_cacao","Sobic.001G175000.1_Sorghum_bicolor","Sevir.9G173600.1_Setaria_viridis","Pahal.9G172500.1_Panicum_hallii","Pahal.1G457700.1_Panicum_hallii",
          "Sobic.004G355200.1_Sorghum_bicolor","Sobic.004G006600.7_Sorghum_bicolor","Sevir.9G171300.4_Setaria_viridis","Pahal.1G006000.1_Panicum_hallii",
          "Sevir.1G006300.1_Setaria_viridis","Medtr5g047620.1_Medicago_truncatula","Medtr5g065200.1_Medicago_truncatula","Medtr1g114070.1_Medicago_truncatula")
tre_reduced=drop.tip(rtre_final,to_drop)
x=as_tibble(tre_reduced)
#View(x)
###########build tree figure
###############complete
# gg_tr=ggtree(rtre_final) +
#   #geom_text(aes(label=node),size=3)+
#   #geom_label(aes(x=branch,label=bootstrap),size=3)+
#   geom_tiplab(align = T,color="black",size=3)+
#   geom_highlight(node=118,fill='red',color='white', alpha=0.1, extend=2)+
#   xlim(0,3)

###############reduced
gg_tr=ggtree(tre_reduced) +
  #geom_text(aes(label=node),size=3)+
  #geom_text(aes(x=branch,label=bootstrap),size=3)+
  geom_tiplab(align = T,color="black",size=3)+
  geom_highlight(node=52,fill='red',color='white', alpha=0.1, extend=2.5)+
  xlim(0,3)

################preparing format for gene structure and methylation infomation

d=read.table("./data/34.structure",header = F)
d2=read.table("./data/34.domain",header = F)
me=read.table("./data/34.methyl",header = F)

colnames(d)=c('label','start','end','element')
colnames(d2)=c('label','start','end','domain')
colnames(me)=c('label','pos','strand','context','mreads','treads','index')
###################apply 
m=select(d, label)
g=select(gg_tr$data, label,y)
tt=seq_along(me$y)
d$y=left_join(select(d, label), select(gg_tr$data, label,y)) %>% pull(y)
d2$y=left_join(select(d2, label), select(gg_tr$data, label,y)) %>% pull(y)
me$y=left_join(select(me, label), select(gg_tr$data, label,y)) %>% pull(y)

me$y1=unlist(lapply(seq_along(me$y), function(i) {
     
     if (me$index[i]==1){
       add=me$mreads[i]/me$treads[i]
       if (me$strand[i] == '+') {
         me$y[i]+add*0.35
       }
       else {
         me$y[i]-add*0.35
       }
     }
     else {
       me$y[i]
     }
    }
  )
)

me$group=unlist(lapply(seq_along(me$y), function(i) {
    # 
      if (me$context[i]=='CAC'| me$context[i]=='CTC'| me$context[i]=='CTT'| me$context[i]=='CTA'| me$context[i]=='CAA'| me$context[i]=='CAT'| me$context[i]=='CCA'| me$context[i]=='CCT'| me$context[i]=='CCC'){
    
        'CHH'
    
      }
      else if (me$context[i]=='CGA'||me$context[i]=='CGT'||me$context[i]=='CGC'||me$context[i]=='CGG'){
        'CG'
      }
      else if (me$context[i]=='CAG'||me$context[i]=='CTG'||me$context[i]=='CCG'){
        'CHG'
      }
      
      else {
        'CHG'
      }
      
    }
  )
)

#View(me)
me2=me[me$index==1,]
###############################################built geom rect
geom_rect5 <- function(mapping = NULL, data = NULL,
                      stat = "identity", position = "identity",
                      ...,
                      linejoin = "mitre",
                      na.rm = FALSE,
                      show.legend = NA,
                      inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRect,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}


GeomRect <- ggproto("GeomRect", Geom,
                    default_aes = aes(colour = NA, fill = "grey35", size = 0.5, linetype = 1,
                                      alpha = NA),
                    
                    required_aes = c("xmin", "xmax", "ymin", "ymax"),
                    
                    draw_panel = function(self, data, panel_params, coord, linejoin = "mitre") {
               
                        coords <- coord$transform(data, panel_params)
                        assign("coords", coords, envir = .GlobalEnv)
                         rect=rectGrob(
                          coords$xmin, coords$ymax,
                          width = coords$xmax - coords$xmin,
                          height = coords$ymax - coords$ymin,
                          default.units = "native",
                          just = c("left", "top"),
                          gp = gpar(
                            col = coords$colour,
                            fill = alpha(coords$fill, coords$alpha),
                            lwd = coords$size * .pt,
                            lty = coords$linetype,
                            linejoin = linejoin,
                            # `lineend` is a workaround for Windows and intentionally kept unexposed
                            # as an argument. (c.f. https://github.com/tidyverse/ggplot2/issues/3037#issuecomment-457504667)
                            lineend = if (identical(linejoin, "round")) "round" else "square"
                          )
                        )
                        
                        ##################calcuate the location of centro line
                         cl=dplyr::group_by(coords, group) %>%
                           dplyr::summarize(lmin = min(xmin)-0.02,lmax=max(xmax)+0.035,
                                            loc=ymin[1]+(ymax[1]-ymin[1])/2
                           )
                         
                         
                        line= segmentsGrob(x0 = cl$lmin, y0 = cl$loc,
                                       x1 = cl$lmax, y1 = cl$loc,
                                       default.units = "native",
                                       arrow = arrow(angle = 30, 
                                                     ends = "last", type = "open",length=unit(0.02,'native')),
                                       name = NULL, gp = gpar(),  vp = NULL)
                        
                        gTree(children = gList(rect,line))

                      }
                    ,
                    
                    draw_key = draw_key_polygon
)

limits <- range(gg_tr$data$y)
expand_limits=expansion(0,.6)
limits[1] <- limits[1] + (limits[1] * expand_limits[1]) - expand_limits[2]
limits[2] <- limits[2] + (limits[2] * expand_limits[3]) + expand_limits[4]


#browseVignettes("ggtree")
#################################################################################################
#########Section 2
#####Analysis for the correlation between intron length and methylation level
#####Analysis for the location of CHG methylation block relatetive to JMJC domain
################################################################################################
library(dplyr)

####get the total /average intron length for each species
####get the methylation level on the intron region for each species

######group the structure dataframe table into lists by label
sp_str=split(d,f=d$label)
intron_ls=list()
for(i in seq_along(sp_str)){
  if(length(sp_str[[i]]$start)>1){
    df <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(df) <-c("label", "start","end","element")
    c=1
    for(j in 2:length(sp_str[[i]]$start)){
      #print(names(sp_str[[i]]))
      end=sp_str[[i]][j,2]
      start=sp_str[[i]][j-1,3]
      if(end-start>1){
        df[c,]=c(label=sp_str[[i]][j,1],start=start,end=end,element="intron")
        c=c+1
        #print(a)
      }
      #break
    }
    intron_ls=append(intron_ls,list(df)) 
  }
}
#################table for location of intron
d_intron=do.call(rbind,intron_ls)
d_intron$start=as.numeric(d_intron$start)
d_intron$end=as.numeric(d_intron$end)
# #############################################average and total intron length
in_aver=d_intron %>%
  group_by(label)%>%
  summarise(aver=mean(end-start),total=sum(end-start))

in_frag=split(d_intron,f=d_intron$label)#############location of each intron fragments
######################reorder methylation file, 
me2=me%>%
  group_by(label)%>%
  arrange(pos)
me_sp=split(me2,f=me2$label)##############group methylation tsv files by label
mreads=list()##########record total methylatied reads and total mapped reads number on intron region for each label
treads=list()
for(i in seq_along(in_frag)){
  #print(names(in_frag)[i])
  #break
  ####find the coresponding group for methylation list
  g=1
  for(k in seq_along(me_sp)){
    if(names(in_frag)[i] == names(me_sp)[k] ){
      g=k
      break
    }
  }
  ############record the context list CG/CHG/CHH
  #mr=list("CHG"=0,"CG"=0,"CHH"=0)
  #tr=list("CHG"=0,"CG"=0,"CHH"=0)
  md=data.frame("start"=rep(0,length(in_frag[[i]]$start)),"end"=rep(0,length(in_frag[[i]]$start)),"CHG"=rep(0,length(in_frag[[i]]$start)),"CG"=rep(0,length(in_frag[[i]]$start)),"CHH"=rep(0,length(in_frag[[i]]$start)))
  td=data.frame("start"=rep(0,length(in_frag[[i]]$start)),"end"=rep(0,length(in_frag[[i]]$start)),"CHG"=rep(0,length(in_frag[[i]]$start)),"CG"=rep(0,length(in_frag[[i]]$start)),"CHH"=rep(0,length(in_frag[[i]]$start)))
  shift=1
  n=length(me_sp[[g]]$pos)
  for(j in 1:length(in_frag[[i]]$start)){
    if(shift<=n){
      from=shift
      to=n
      for (p in from:to){
        if (me_sp[[g]]$pos[p]>=in_frag[[i]]$start[j] & me_sp[[g]]$pos[p]<=in_frag[[i]]$end[j]) {
          #print(me_sp[[g]]$group[p])
          #mr[[me_sp[[g]]$group[p]]]=mr[[me_sp[[g]]$group[p]]]+me_sp[[g]]$mreads[p]
          #tr[[me_sp[[g]]$group[p]]]=tr[[me_sp[[g]]$group[p]]]+me_sp[[g]]$treads[p]
          md[j,me_sp[[g]]$group[p]]=md[j,me_sp[[g]]$group[p]]+me_sp[[g]]$mreads[p]
          td[j,me_sp[[g]]$group[p]]=td[j,me_sp[[g]]$group[p]]+me_sp[[g]]$treads[p]
        }
        else if (me_sp[[g]]$pos[p] > in_frag[[i]]$end[j]){
          md[j,"start"]=in_frag[[i]]$start[j]
          md[j,"end"]=in_frag[[i]]$end[j]
          td[j,"start"]=in_frag[[i]]$start[j]
          td[j,"end"]=in_frag[[i]]$end[j]
          shift=p+1
          break
        }
      }
    }
  }
  mreads[[names(in_frag)[i]]]=md
  treads[[names(in_frag)[i]]]=td
  #break
}

db_mr=bind_rows(mreads,.id = "label")
colnames(db_mr)[4:6]=c("mCHG","mCG","mCHH")

db_tr=bind_rows(treads,.id = "label")

d_full=full_join(db_mr,db_tr,by=c("label","start","end"))

#################get the averaged methylation value as probabilty of methylation for each gene

pm=d_full %>% 
  group_by(label) %>%
  summarise(rCHG=sum(mCHG)/sum(CHG),rCG=sum(mCG)/sum(CG),rCHH=sum(mCHH)/sum(CHH))

#########combine ratio together
df2=full_join(d_full,pm,by="label")

df=full_join(in_aver,pm,by="label")
####add pvalue from binominal test


add_pvalue=function(df,cols,limit=50){
  for(i in 1:length(df[,1])){
    if(df[i,cols[1]]<limit){
      #print(df2[i,])
      df[i,cols[4]]=1
    }
    else{
      df[i,cols[4]]=binom.test(df[i,cols[1]], df[i,cols[2]], df[i,cols[3]],alternative = "greater")$p.value
      #print(df2[i,'pCHG'])
    }
  }
  return(df)
}

l1=c("mCHG","CHG","rCHG","pCHG")
l2=c("mCG","CG","rCG","pCG")
l3=c("mCHH","CHH","rCHH","pCHH")
df3=add_pvalue(df2,l1)
df3=add_pvalue(df3,l2)
df3=add_pvalue(df3,l3)

########################################calculate the relative location to jmjc

#df3[df3$label=="16425412_Carica_papaya",]
#View(d_intron)


#########decided its location relative to JmJC domain
#####using binorminal test to decide higher CG/CHG methylation intron and their relative locations to JmJC domain for each gene, 
#########generate the list
dm=d2 %>% 
  group_by(label) %>% 
  summarise(s=min(start),e=max(end))
  
for(i in 1:length(df3[,1])){
  if(df3[i,"pCHG"]<0.01){
    #break
    if(df3[i,"end"]<dm[dm$label==df3[i,1],]$s){
      df3[i,"status"]="before"

    }
    else if (df3[i,"start"]>dm[dm$label==df3[i,1],]$e){
      df3[i,"status"]="after"

    }
    else {
      df3[i,"status"]="within"
      }

  }
  else {
    df3[i,"status"]="no"
  }
}

write.csv(df3,"./intron.status.csv",quote = F,row.names = F)
##########
#####################################
df4=df3 %>%
  group_by(label) %>%
  count(status) %>%
  count(label)%>%
  ungroup()%>%
  summarize(no=length(label[n==1]),mintron=length(label[n>1]))%>%
  pivot_longer(cols=no:mintron,names_to = "group",values_to = "count")

df5_sp=df3 %>%
  group_by(label) %>%
  count(status) %>%
  filter(status!="no")%>%
  summarize(num=sum(n))
  

df5=df3 %>%
  group_by(label) %>%
  count(status) %>%
  filter(status!="no") %>%
  ungroup()%>%
  summarise(before=sum(n[status=='before']),within=sum(n[status=='within']),after=sum(n[status=='after']))%>%
  pivot_longer(cols=1:3,names_to = "group",values_to = "count")


#######################plot
d_y=d %>%
  group_by(label)%>%
  distinct(label,y)

df3y=left_join(select(df3,label,start,end,status),d_y,by="label")%>%
  filter(status!="no")

cp=c("#adcf9f","#37ae7d")
met=c("#FFC4C4","#0057e7","#ffa700")
me2$group=factor(me2$group,levels=c("CG","CHG","CHH"))

gg_st=ggplot()+
  geom_segment(data=me2,aes(x=pos,y=y,xend=pos,yend=y1,color=group),size=0.1,alpha=0.9)+
  geom_rect5(data=d,aes(xmin=start,xmax=end,ymin=y-0.12,ymax=y+0.12,group=label,fill=element),alpha=1)+
  geom_rect(data=d2,aes(xmin=start,xmax=end,ymin=y-0.12,ymax=y+0.12),fill='#F94892',alpha=1)+
  geom_segment(data=df3y,aes(x=start,y=y-0.45,xend=end,yend=y-0.45),linewidth=1.5,color="#19A7CE")+
  scale_y_continuous(limits=limits, expand=c(0,0))+
  scale_fill_manual(values=cp)+
  scale_color_manual(values=met)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab("")+
  ylab("")


p1=gg_tr + 
  gg_st +
  plot_annotation(tag_levels="A")+
  plot_layout(ncol=2,widths=c(1,1.5))



pdf('test.bombine_reduced.pdf',
    width=12,
    height=12)
p1

dev.off()

df3_group=df3 %>%
  mutate(group=if_else(status=='no',"no","mintron"),.keep="all")

p2=ggplot(df3_group,aes(x=group,y=(end-start)/1000,fill=group))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  scale_y_continuous(limits=c(0,10))+
  theme_classic()+
  xlab("mCHG intron")+
  ylab("intron length (kb)")+
  scale_fill_brewer(palette="Set1")+
  theme(legend.position = "none")  


p4=ggplot(df4,aes(x=group,y=count,fill=group))+
  geom_bar(stat="identity")+
  geom_text(aes(label=count),  color = "white",
             show.legend = FALSE,size=8,vjust = 1.5)+
  theme_classic() +
  scale_fill_brewer(palette="Set1")+
  xlab("mCHG intron")+
  ylab("the number of IBM1 genes")+
  theme(legend.position = "none")

p5=ggplot(df5,aes(x=group,y=count,fill=group))+
  geom_bar(stat="identity")+
  geom_text(aes(label=count),  color = "black", show.legend = FALSE,size=8,vjust = 1.5)+
  theme_classic() + 
  scale_fill_manual(values = c("#FFD93D","#FF8400","#159895"),name="position to JmJC")+
  theme(legend.position = "none")+
  xlab("")+
  ylab("the number of methylated introns")+
  


pdf('result.summary.pdf',
    width=8,
    height=4)
plot_grid(p2, p4,p5,nrow=1)

dev.off()
