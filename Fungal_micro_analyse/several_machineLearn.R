library("phyloseq")
 #b比较多种机器学习方法对枯萎病发病健康群落的分类预测效果
# saveRDS(ps7,"./ps_model.rds")
ps7 = readRDS("./ps_final_model.rds")
ps7 
# ps7 = filter_taxa(ps7, function(x) sum(x ) > 0 , TRUE);ps7
# ps7  = transform_sample_counts(ps7, function(x) x / sum(x) );ps7



mapping = as.data.frame(sample_data(ps7))
mapping$SampleType
table(mapping$SampleType)
table(mapping$Description)
vegan_otu <-  function(physeq){
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}
otutab = as.data.frame(t(vegan_otu(ps7)))
dim(otutab)
head(otutab)
mapping7 = as.data.frame(sample_data(ps7))
head(mapping7)




####随机森林--这是我选择的机器学习方法
library(randomForest)
library(caret)
library(pROC)
otutab_t = as.data.frame(t(otutab))
# Set classification info.
otutab_t$group = factor(mapping7$SampleType,levels= c("D_bulk","H_bulk"))
colnames(otutab_t) = paste("OTU",colnames(otutab_t),sep = "")
# set random seed for reproducible
set.seed(315)##不错···········目前最好
# RandomForest Classification
model_add= randomForest(OTUgroup ~ ., data=otutab_t, importance=TRUE, proximity=TRUE)
print(model_add)


# model = model_orig
ps6 = readRDS("./add_sample_for_predict/owe//a9_usearch_otu_table/ps_owe.rds")
ps6
map = as.data.frame(sample_data(ps6))
map
map$SampleType = paste(map$SampleTypeDH,map$zone,sep = "_")
sample_data(ps6) = map

result = predict_rand(ps7,ps6,model = model)
result[[2]]
result[[1]]




##下面对多种机器学习方法进行比对


test = otutab_t
#####随机森林###########################################
#将分组变量提到第一列
library(tidyverse)
test = select(test,OTUgroup,everything())
head(test)
train = test
folds<-createFolds(y=test[,1],k=5)
folds
max=0
num=0
fc<-as.numeric()
mod_pre<-as.numeric()
for(i in 1:5){
  fold_test<-train[folds[[i]],]
  fold_train<-train[-folds[[i]],]
  model<-randomForest(OTUgroup~.,data=fold_train, importance=TRUE, proximity=TRUE)
  model_pre<-predict(model,newdata = fold_test,type="prob")
  fc<-append(fc,as.factor(fold_test$OTUgroup))
  mod_pre<-append(mod_pre,model_pre[,1])
}
df<-cbind(fc,as.numeric(mod_pre))

# a = df









###############################SVM#######################
library(e1071)
max=0
num=0
fc<-as.numeric()
mod_pre<-as.numeric()
i = 1
for(i in 1:5){
  fold_test<-train[folds[[i]],]
  # head(fold_test)
  fold_train<-train[-folds[[i]],]
  model<-svm(OTUgroup~.,data=fold_train,probability=TRUE)
  model
  model_pre<-predict(model,newdata = fold_test,decision.values = TRUE, probability = TRUE)
  fc<-append(fc,as.numeric(fold_test$OTUgroup))
  mod_pre<-append(mod_pre,as.numeric(attr(model_pre, "probabilities")[,2]))
}
df<-cbind(df,cbind(fc,mod_pre))



# 
# model<-svm(OTUgroup~.,data=fold_train,probability=TRUE)
# model
# model_pre<-predict(model,newdata = fold_test,decision.values = TRUE, probability = TRUE)
# fc<-append(fc,as.numeric(fold_test$OTUgroup))
# mod_pre<-append(mod_pre,as.numeric(attr(model_pre, "probabilities")[,2])) 
# 
# pre_svm <- predict(model,type='response',newdata = fold_test)
# pre_svm 
# obs_p_svm = data.frame(prob=pre_svm,obs=fold_test$OTUgroup)
# ###输出混淆矩阵
# table(fold_test$OTUgroup,pre_svm,dnn=c("真实值","预测值"))
# ###绘制ROC曲线
# svm_roc <- roc(fold_test$OTUgroup,as.numeric(pre_svm))
# 
# unclass(svm_roc)
# plot(svm_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM模型ROC曲线 kernel = radial')


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          x <- 1:10        
attr(x,"dim") <- c(2, 5)

x







####################逻辑回归
max=0
num=0
fc<-as.numeric()
mod_pre<-as.numeric()
for(i in 1:5){
  fold_test<-train[folds[[i]],]
  fold_train<-train[-folds[[i]],]
  model<-glm(OTUgroup~.,family='binomial',data=fold_train)
  model
  model_pre<-predict(model,type='response',newdata=fold_test)
  model_pre
  # pred <- prediction( model_pre,fold_test$OTUgroup)
  # performance(pred,'auc')@y.values #AUC值
  # perf <- performance(pred,'tpr','fpr')
  # plot(perf)
  fc<-append(fc,fold_test$OTUgroup)
  mod_pre<-append(mod_pre,as.numeric(model_pre))
}


df<-cbind(df,cbind(fc,mod_pre))
df




# setwd("E:/D/Shared_Folder/DATA_get_wilt_analyse/model_for_online_monitor/Fwilt_bac_16s")


pdf("./ROC—ITS.pdf",height=14,width=14)
par(cex = 4)
mycol <- c("slateblue","sea green3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")
x<-plot.roc(df[,1],df[,2],
            smooth=FALSE,
            lwd=4,
            ylim=c(0,1),
            xlim=c(1,0),
            legacy.axes=TRUE,
            main="",
            col=mycol[2])

x<-plot.roc(df[,3],df[,4],
            smooth=FALSE,
            add=TRUE,
            lwd=4,
            ylim=c(0,1),
            xlim=c(1,0),
            legacy.axes=TRUE,
            main="",
            col=mycol[3])

x<-plot.roc(df[,5],df[,6],
            smooth=FALSE,
            add=TRUE,
            lwd=4,
            ylim=c(0,1),
            xlim=c(1,0),
            legacy.axes=TRUE,
            main="",
            col=mycol[4])

legend.name <- c(paste("RF","AUC",0.96,sep=" "),paste("SVM","AUC",0.87,sep=" "),paste("LR","AUC",0.76,sep=" "))
legend("bottomright", 
       legend=legend.name,
       col = mycol[2:4],
       lwd = 2,
       bty="n")


dev.off()




pdf("./ROC_RF.pdf",height=14,width=14)
par(cex = 4)
mycol <- c("slateblue","sea green3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")
x<-plot.roc(df[,1],df[,2],
            smooth=FALSE,
            lwd=4,
            ylim=c(0,1),
            xlim=c(1,0),
            legacy.axes=TRUE,
            main="",
            col=mycol[2])


legend.name <- c(paste("RF","AUC",0.9677,sep=" "))
legend("bottomright", 
       legend=legend.name,
       col = mycol[2:4],
       lwd = 2,
       bty="n")
dev.off()










##########################
ps6 = readRDS("./add_sample_for_predict/owe//a9_usearch_otu_table/ps_owe.rds")

map = as.data.frame(sample_data(ps6))
map
map$SampleType = paste(map$SampleTypeDH,map$zone,sep = "_")
sample_data(ps6) = map
result = predict_rand(ps7,ps6,model = model)

result[[2]]
test = result[[4]]
test = select(test,OTUgroup,everything())
head(test)



model<-glm(OTUgroup~.,data=train,family=binomial(link=logit))
model_pre<-predict(model,type='response',newdata=test)
pdf("/Users/baiyunfan/desktop/RF.pdf",height=6,width=6)
x<-plot.roc(test[,1],model_pre,
            smooth=F,
            lwd=2,
            ylim=c(0,1),
            xlim=c(1,0),
            legacy.axes=T,
            main="",
            col=mycol[4])


legend.name <- paste("LR","AUC",0.9557,sep=" ")
legend("bottomright", 
       legend=legend.name,
       col = mycol[4],
       lwd = 2,
       bty="n")
dev.off()



##贝叶斯
library(e1071)
max=0
num=0
fc<-as.numeric()
mod_pre<-as.numeric()
i = 1
fold_test<-train[folds[[i]],]
 # head(fold_test)
fold_train<-train[-folds[[i]],]
model<-naiveBayes(OTUgroup~.,data=fold_train)
pre_svm <- predict(model,newdata = fold_test)

library(ROCR)
obs_p_svm = data.frame(prob=pre_svm,obs=fold_test$OTUgroup)
###输出混淆矩阵
table(fold_test$OTUgroup,pre_svm,dnn=c("真实值","预测值"))
###绘制ROC曲线
svm_roc <- roc(fold_test$OTUgroup,as.numeric(pre_svm))
#求取roc值
svm_roc$auc
(auc1 = auc(svm_roc))
plot(svm_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM模型ROC曲线 kernel = radial')



##逻辑回归
i = 2
fold_test<-train[folds[[i]],]
fold_train<-train[-folds[[i]],]
model<-glm(OTUgroup~.,family='binomial',data=fold_train)
model
model_pre<-predict(model,type='response',newdata=fold_test)
model_pre

library(ROCR)
obs_p_svm = data.frame(prob=model_pre,obs=fold_test$OTUgroup)
###输出混淆矩阵
table(fold_test$OTUgroup,model_pre,dnn=c("真实值","预测值"))
###绘制ROC曲线
svm_roc <- roc(fold_test$OTUgroup,as.numeric(model_pre))
#求取roc值
svm_roc$auc#0.7675
(auc1 = auc(svm_roc))


###############svm#########

model<-svm(OTUgroup~.,data=fold_train,probability=TRUE)
model
model_pre<-predict(model,newdata = fold_test,type='response')
library(ROCR)
obs_p_svm = data.frame(prob=model_pre,obs=fold_test$OTUgroup)
###输出混淆矩阵
table(fold_test$OTUgroup,model_pre,dnn=c("真实值","预测值"))
###绘制ROC曲线
svm_roc <- roc(fold_test$OTUgroup,as.numeric(model_pre))
#求取roc值
svm_roc$auc#0.8797
(auc1 = auc(svm_roc))


###随机森林
i= 1
fold_test<-train[folds[[i]],]
fold_train<-train[-folds[[i]],]
model<-randomForest(OTUgroup~.,data=fold_train, importance=TRUE, proximity=TRUE)
model_pre<-predict(model,newdata = fold_test,type='response')
library(ROCR)
obs_p_svm = data.frame(prob=model_pre,obs=fold_test$OTUgroup)
###输出混淆矩阵
table(fold_test$OTUgroup,model_pre,dnn=c("真实值","预测值"))
###绘制ROC曲线
svm_roc <- roc(fold_test$OTUgroup,as.numeric(model_pre))
#求取roc值
svm_roc$auc#0.9677
(auc1 = auc(svm_roc))





##mboost算法很慢，所以先不运行

library(e1071)
max=0
num=0
fc<-as.numeric()
mod_pre<-as.numeric()
i = 1
for(i in 1:5){
  fold_test<-train[folds[[i]],]
  # head(fold_test)
  fold_train<-train[-folds[[i]],]
  model<-blackboost(OTUgroup~.,data=fold_train,family = Binomial())
  model
  model_pre<-predict(model,newdata = fold_test,decision.values = TRUE, probability = TRUE)
  fc<-append(fc,as.numeric(fold_test$OTUgroup))
  mod_pre<-append(mod_pre,as.numeric(model_pre))
}
df<-cbind(df,cbind(fc,mod_pre))

dim(model_pre)
df






