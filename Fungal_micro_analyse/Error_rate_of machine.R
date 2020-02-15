



wt  = read.csv("./Focmodel_rate.csv")
head(wt)


library(reshape2)

wt2 = melt(wt,id.vars = "X",variable.name = "group",value.name = "rate")
head(wt2)
wt2$rate = 1- wt2$rate
p = ggplot(wt2,aes(x = X,y = rate,group = group,color = group,fill = group)) + 
  # geom_bar(aes(fill = group),stat = "identity",position = "dodge",color = "black")+
  geom_line(size = 1)+geom_point(pch = 21,size = 4,color ="black" )+
  theme_classic() 
p = p +mytheme

p

p=p+theme(axis.text.x=element_text(angle=60,vjust=1, hjust=1))

FileName2 <- paste("rote_of_different_FUN_tax_model1",".pdf", sep = "")

ggsave(FileName2, p, width = 10, height =8, device = cairo_pdf, family = "Times New Roman" )

