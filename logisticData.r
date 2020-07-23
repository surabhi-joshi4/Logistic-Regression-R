


########################Data import#######################

data(BreastCancer, package="mlbench")
bca<- BreastCancer[complete.cases(BreastCancer), ] 
str(bca)

############################## Header and Univariate Analysis#############################
names(bca)
attach(bca)
summary(bca)
head(bca,10)
table(bca$Class)


#########################if else use for class=malignant for result######################

bca$Class <- ifelse(bca$Class == "malignant", 1, 0)
bca$Class <- factor(bca$Class, levels = c(0, 1))


logit <- glm(Class ~ Cell.shape, data=bca,family='binomial')
summary(logit)
anova(logit,test='Chisq')

################################Data partition###################
set.seed(45)
library(caret)
Train <- createDataPartition(bca$Class, p=0.7, list=FALSE)
training <- bca[ Train, ]
testing <- bca[ -Train, ]


summary(training)

############################## Missing Value #############################
sapply(training,function(x) sum(is.na(x)))


###### Identify outlier of Cell.shape Variable #####
boxplot(training$Cell.shape)


bca$Class <- ifelse(bca$Class == "malignant", 1, 0)
bca$Class <- factor(bca$Class, levels = c(0, 1))

logit <- glm(Class ~ Cell.shape+Cl.thickness+Cell.size,family='binomial',
             data=training)
summary(logit)

###################### Variable Selection Method #####################  

logit2 <- step(glm(Class ~ Cell.shape+Cl.thickness+Cell.size,family='binomial', data=training),direction = "both")
summary(logit2)
anova(logit2,test='Chisq')

#######################check concordance and discordance of model#################


Acc=function(model){
  Data = cbind(model$y, model$fitted.values) 
  ones = Data[Data[,1] == 1,]
  zeros = Data[Data[,1] == 0,]
  conc=matrix(0, dim(zeros)[1], dim(ones)[1])
  disc=matrix(0, dim(zeros)[1], dim(ones)[1])
  ties=matrix(0, dim(zeros)[1], dim(ones)[1])
  for (j in 1:dim(zeros)[1])
  {
    for (i in 1:dim(ones)[1])
    {
      if (ones[i,2]>zeros[j,2])
      {conc[j,i]=1}
      else if (ones[i,2]<zeros[j,2])
      {disc[j,i]=1}
      else if (ones[i,2]==zeros[j,2])
      {ties[j,i]=1}
    }
  }
  Pairs=dim(zeros)[1]*dim(ones)[1]
  PercentConcordance=(sum(conc)/Pairs)*100
  PercentDiscordance=(sum(disc)/Pairs)*100
  PercentTied=(sum(ties)/Pairs)*100
  return(list("Percent Concordance"=PercentConcordance,"Percent Discordance"=PercentDiscordance,"Percent Tied"=PercentTied,"Pairs"=Pairs))
}

Acc(logit2)

##################### # odds Ratio ##################### 
exp(coef(logit2))

##################### ## Prediction on testing data ##################### 
testing$probs <-predict(logit2, testing, type='response')
testing$Predict<-as.factor(ifelse(testing$probs>0.70,1,0))

###################### Accuracy of testing data  #####################
library(e1071)
table(testing$Predict, testing$Class)
confusionMatrix(testing$Predict,testing$Class,positive = "1")

######################## Roc Curve  #####################
library(ROCR)
# Make predictions on training set
predictTrain = predict(logit2,testing, type="response")
# Prediction function
ROCRpred = prediction(predictTrain, testing$Class)
# Performance function
ROCRperf = performance(ROCRpred, "tpr", "fpr")
# Plot ROC curve
plot(ROCRperf)

######################### AUC ( index of accuracy) #####################

library(ROCR)
pred = prediction(testing$probs, testing$Class)
as.numeric(performance(pred, "auc")@y.values)

