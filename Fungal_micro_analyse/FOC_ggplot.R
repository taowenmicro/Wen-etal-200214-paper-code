
mythemealpha <- theme_bw()+
  #theme_classic()+
  # scale_color_manual(values = mi, guide = guide_legend(title = NULL))+
  # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
  theme(
    
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    
    plot.title = element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =element_text(size = 20,face = "bold",colour = "black"),
    axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
    axis.text = element_text(size = 20,face = "bold"),
    axis.text.x = element_text(colour = "black",size = 20),
    axis.text.y = element_text(colour = "black",size = 20),
    legend.text = element_text(size = 20,face = "bold"))+
  # theme(legend.position = "top")+
  
  theme(strip.text.x = element_text(size=15, angle=0),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="blue", fill="#CCCCFF"))
library(RColorBrewer)#调色板调用包

#调用所有这个包中的调色板
display.brewer.all()
#提取特定个数的调色板颜色，会出图显示
display.brewer.pal(8,"Set1")
##仅仅只显示色号,我们要在颜色上并显示色号才好
mi = brewer.pal(12,"Set1")
library("scales")
show_col(mi)



path = getwd()

dir.create("./Fon_diversity")

path = paste(path,"/Fon_diversity/",sep = "")

path





library(tidyverse)


ps1 = readRDS("./ps_final_model.rds")



ps1
ps1_rela  = transform_sample_counts(ps1, function(x) x / sum(x) );ps1_rela 



ps2 <- ps1_rela %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    Genus  == "Fusarium"
    #Species == "Fusarium_oxysporum"
  )
ps2
plot_barwt <- function (physeq, x = "Sample", y = "Abundance", fill = NULL,
                        title = NULL, facet_grid = NULL)
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

p <- plot_barwt(ps2, x="SampleType", fill="Species") +
  scale_x_discrete(limits = c(axis_order))
p <- p+mythemealpha
p

colnames(tax_table(ps2))

p <- plot_barwt(ps2, x="SampleType", fill="Genus") +
  scale_x_discrete(limits = c(axis_order))
p <- p+mythemealpha
p


FileName2 <- paste(path,"镰刀菌属",".pdf", sep = "")
library("Cairo")
ggsave(FileName2, p, width = 3, height =6, device = cairo_pdf, family = "Times New Roman" )





#可视化镰刀菌

ps2 <- ps1_rela %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    Species == "Fusarium_oxysporum"
  )
ps2

p <- plot_barwt(ps2, x="SampleType", fill="Species") +
  scale_x_discrete(limits = c(axis_order),labels = xlabel)
p <- p+mythemealpha +scale_fill_manual(values = mi)
p




FileName2 <- paste(path,"镰刀菌病原菌",".pdf", sep = "")
library("Cairo")
ggsave(FileName2, p, width = 6, height =6, device = cairo_pdf, family = "Times New Roman" )



#可视化镰刀菌和模型中的林外一株菌
xlabel = c("Disease","Health")
axis_order = c( "D_bulk","H_bulk")

ps2 <- ps1_rela %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    Species %in%c("Fusarium_oxysporum","Fusarium_keratoplasticum") 
  )
ps2

p <- plot_barwt(ps2, x="SampleType", fill="Species") +
  scale_x_discrete(limits = c(axis_order),labels = xlabel)
p <- p+mythemealpha+scale_fill_manual(values = c(mi[2],mi[1]))
p
FileName2 <- paste(path,"镰刀菌-林外一个",".pdf", sep = "")
library("Cairo")
ggsave(FileName2, p, width = 6.5, height =6, device = cairo_pdf, family = "Times New Roman" )

#可视化话镰刀菌和随机森林中另外一株菌
ps2 <- ps1_rela %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    # Species %in%c("Fusarium_oxysporum","Fusarium_keratoplasticum") 
    row.names(tax_table(ps1_rela ))%in%c("SH010924.07FU_KF986690_reps_singleton","SH020983.07FU_JN235282_refs")
  )
ps2

mapping = as.data.frame(sample_data(ps2))
otu_table = as.data.frame((vegan_otu(ps2)))
head(otu_table)
index = merge(otu_table,mapping,by = "row.names",all = F)
dim(index)

tax_table = as.data.frame(vegan_tax(ps2))
head(tax_table)

p = ggplot(index,aes(x =otu_table$SH010924.07FU_KF986690_reps_singleton,y=  otu_table$SH020983.07FU_JN235282_refs ,fill = index$SampleType,color= index$SampleType))+ 
geom_point(position = position_jitter(width = 0.3,height =0.03),pch = 21,size = 2) + 
  # stat_ellipse(type = "t", linetype = 2,level = 0.99) +
  scale_color_manual(values = c("#CC6666","#7777DD")) +
  labs(title="", 
       x='Fusarium_oxysporum', y='Fusarium_keratoplasticum') +
theme_classic()
p

FileName2 <- paste(path,"病原菌-另外镰刀菌",".pdf", sep = "")
library("Cairo")
ggsave(FileName2, p, width = 6, height =6, device = cairo_pdf, family = "Times New Roman" )


  # +geom_jitter(position=position_jitter(0.15))#+ stat_smooth(method = lm)+





p <- plot_barwt(ps2, x="SampleType", fill="Species") +
  scale_x_discrete(limits = c(axis_order))
p <- p+mythemealpha
p


#可视化非镰刀菌


ps2 <- ps1_rela %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    Genus  == "Fusarium"&
      Species != "Fusarium_oxysporum"
  )
ps2

p <- plot_barwt(ps2, x="SampleType", fill="Species") +
  scale_x_discrete(limits = c(axis_order))
p <- p+mythemealpha
p


FileName2 <- paste(path,"病原菌-其他镰刀菌",".pdf", sep = "")
library("Cairo")
ggsave(FileName2, p, width = 6, height =6, device = cairo_pdf, family = "Times New Roman" )





ps1 <- ps1_rela %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    Genus  == "Fusarium"
    
  )
ps1


vegan_otu <-  function(physeq){
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}
ps1 = tax_glom(ps1,taxrank ="Species")
ps1

otu_table = as.data.frame(t(vegan_otu(ps1)))
head(otu_table)
vegan_tax <-  function(physeq){
  tax <-  tax_table(physeq)
  
  return(as(tax,"matrix"))
}
tax_table = as.data.frame(vegan_tax(ps1))
head(tax_table)
dim(tax_table)
index = merge(otu_table,tax_table,by = "row.names",all = F)
row.names(index) = index$Species
index$Row.names = NULL
dim(index)
index = index[1:(dim(index)[2]-dim(tax_table)[2])]
index = t(index)
head(index)


wt = index

# library("corrplot")
library(ggcorrplot)
# install.packages("ggcorrplot")
library(igraph)
library(psych)
# corr.test求相关性矩阵
occor = corr.test(wt,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
# ?corr.test

occor.r = occor$r # 取相关性矩阵R值

occor.p = occor$p # 取相关性矩阵p值
occor.r[occor.p>0.05] = 0
# occor.r[occor.p>0.05|abs(occor.r)<0.1] = 0 

#我是通过将otu和通路分别放到两个坐标轴实现图形空间的压缩
# wt3<-occor.r[c((dim(otu_more)[2]+1):dim(occor.r )[1]),c(1:dim(otu_more)[2])]
# wt3 = t(wt3)

# pdf("./a3_RDA_CCA_cor/cor_moreOTU_env.pdf")
# tiff(file="./a3_RDA_CCA_cor/cor_moreOTU_env.tif", res = 300, compression = "none", width=360,height=280,units= "mm")
library("corrplot")
wt3 = occor.r 
corrplot(wt3,cl.cex = 0.4,tl.cex = 0.6,tl.srt = 45)#,addCoef.col = "black",
library("ggcorrplot")
###方形可视化
ggcorrplot(wt3)


#  圆形可视化
p = ggcorrplot(wt3, method = "circle",outline.color	= "white",lab = TRUE)

p







