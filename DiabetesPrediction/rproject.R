library('rpart')
library(caret)
library(party)
library(varhandle)
library(scales)

# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Rscript rproject.R --fold n --train data.csv --report performance.csv --predict predict.csv", call.=FALSE)
}

# parse parameters
j<-grep("-",c(args[1:length(args)]))


if(length(j)!=4){
  stop("you should type 4 argument", call.=FALSE)
}else{
  
  if((args[j[1]]=="--fold")|(args[j[2]]=="--fold")|(args[j[3]]=="--fold")|(args[j[4]]=="--fold")){
    foldcheck<-1
  }else{
    stop("missing fold command", call.=FALSE)
  }
  
  
  if((args[j[1]]=="--train")|(args[j[2]]=="--train")|(args[j[3]]=="--train")|(args[j[4]]=="--train")){
    traincheck<-1
  }else{
    stop("missing train.csv", call.=FALSE)
  }
  
  
  if((args[j[1]]=="--report")|(args[j[2]]=="--report")|(args[j[3]]=="--report")|(args[j[4]]=="--report")){
    reportcheck<-1
  }else{
    stop("missing --report performance.csv", call.=FALSE)
  }
  
  
  if((args[j[1]]=="--predict")|(args[j[2]]=="--predict")|(args[j[3]]=="--predict")|(args[j[4]]=="--predict")){
    predictcheck<-1
  }else{
    stop("missing --predict predict.csv", call.=FALSE)
  }
  
}



for(i in 1:4){
  judge<-j[i]
  
  if(args[judge]=="--fold"){
    FOLD<-args[judge+1]
  }else if(args[judge]=="--train"){
    TRdf<-args[judge+1]
  }else if(args[judge]=="--report"){
    Re_out<-args[judge+1]
  }else if(args[judge]=="--predict"){
    PRE_out<-args[judge+1]
  }else{
    stop("Unknown flag", call.=FALSE)
  }
  
}



FOLD<-as.numeric(FOLD)
gy <-read.csv(TRdf,na.strings="")

df <- subset(gy, select = c(Pregnancies, Glucose, BloodPressure, SkinThickness, BMI, Age,Outcome))



#####################################################
#boxplot(data$Age)

#Pregnancies   delete  >=14.5
#Glucose         delete  <=10
#BloodPressure  delete   <=10
#SkinThickness  delete   >=80
#Insulin  delete      >=600 
#BMI   delete    <=10 and >=60
#DiabetesPedigreeFunction  delete  
#Age  delete   >=75

#Outcome
df <- df[df$Pregnancies <= 14.5,]
df <- df[df$Glucose >= 10,]
df <- df[df$BloodPressure >= 10,]
df <- df[df$SkinThickness <= 80 ,]
#df <- df[df$Insulin <= 600 ,]
df <- df[df$BMI <= 60 & df$BMI >= 10,]
#df <- df[df$DiabetesPedigreeFunction <= 65,]
df <- df[df$Age <= 75,]



####################################################
#df$f1<- rescale(df$Pregnancies)
#df$f2<- rescale(df$Glucose)
#df$f3<- rescale(df$BloodPressure)
#df$f4<- rescale(df$SkinThickness)
#df$f5<- rescale(df$BMI)
#df$f6<- rescale(df$Age)

#df <- subset(df, select = c(df$f1,df$f2,df$f3,df$f4,df$f5,df$f6,df$Outcome))


#####################################################
#shuffle the data
set.seed(9850)
g<-runif(nrow(df))
data<-df[order(g),]

RR<-dim(data)[1]
CC<-dim(data)[2]
################################################################
#data$Pregnancies[data$Pregnancies==""] <- NA
#data$Glucose[data$Glucose==""] <- NA
#data$BloodPressure[data$BloodPressure==""] <- NA
#data$SkinThickness[data$SkinThickness==""] <- NA
#data$Insulin[data$Insulin==""] <- NA
#data$BMI[data$BMI==""] <- NA
#data$DiabetesPedigreeFunction[data$DiabetesPedigreeFunction==""] <- NA
##data$age[data$age==""] <- NA

#data$Pregnancies[data$Pregnancies=="NA"] <- NA
#data$Glucose[data$Glucose=="NA"] <- NA
#data$BloodPressure[data$BloodPressure=="NA"] <- NA
#data$SkinThickness[data$SkinThickness=="NA"] <- NA
#data$Insulin[data$Insulin=="NA"] <- NA
#data$BMI[data$BMI=="NA"] <- NA
#data$DiabetesPedigreeFunction[data$DiabetesPedigreeFunction=="NA"] <- NA
##data$age[data$age=="NA"] <- NA

#sapply(data,function(x) sum(is.na(x)))
##########################################

#####################################################################
#interval of each fold
interval<-as.integer(RR/FOLD)

#get start point of test data
start<-1

#get end point of test data
end<-interval

###################################################3
#get test data  #real test data #no train
Test<-data[start:end,]

#get input data #real input #need to be split into train data and validation part
INPUT<-data[(end+1):RR,]
Realfold<-FOLD-1
####################################################

#print(class(INPUT$Pregnancies))
#print(class(INPUT$Glucose))
#print(class(INPUT$BloodPressure))
#print(class(INPUT$SkinThickness))
#print(class(INPUT$Insulin))
#print(class(INPUT$BMI))
#print(class(INPUT$DiabetesPedigreeFunction))
#print(class(INPUT$Age))

#print(class(INPUT$Outcome))   #integer
INPUT$Outcome=as.factor(INPUT$Outcome)
Test$Outcome=as.factor(Test$Outcome)
######################################################
#model1
DTree<-function(TRData,VAData)
{
  col<-dim(TRData)[2]
  
  TRgoal<-TRData[,col]
  TRattr<-TRData[,1:(col-1)]
  
  TEgoal<-VAData[,col]
  TEattr<-VAData[,1:(col-1)]
  
  fit<- rpart(TRgoal~., data =TRattr, method = 'class',model=TRUE,control=rpart.control(maxdepth=4))
  p1<-predict(fit,TEattr,type = "class")
  
  t1<-table(TEgoal,p1)
  
  acc1<-sum(diag(t1))/sum(t1)
  acc1<-format(round(acc1,2),nsmall=2)
  
  return(acc1)
}


#save2<-0

#model2
#RF<-function(TRData2,VAData2)
#{
#  set.seed(415)
#  col2<-dim(TRData2)[2]
  
#  TRgoal2<-TRData2[,col2]
#  TRattr2<-TRData2[,1:(col2-1)]
  
#  TEgoal2<-VAData2[,col2]
#  TEattr2<-VAData2[,1:(col2-1)]
  
  
#  assign("save2", TRgoal2, envir = .GlobalEnv)
  
#  fit2 <- cforest(formula=save2~.,data = TRattr2, controls=cforest_unbiased(ntree=2000, mtry=5))
#  p2 <- predict(object = fit2, newdata = TEattr2, OOB=TRUE, type = "response")
  
#  cla<-ifelse(p2>=0.5,1,0)
  
#  t2<-table(TEgoal2,cla)
  
#  acc2<-sum(diag(t2))/sum(t2)
#  acc2<-format(round(acc2,2),nsmall=2)
  
#  return(acc2)
#}

#model3
GLM<-function(TRData3,VAData3)
{
  col3<-dim(TRData3)[2]
  
  TRgoal3<-TRData3[,col3]
  TRattr3<-TRData3[,1:(col3-1)]
  
  TEgoal3<-VAData3[,col3]
  TEattr3<-VAData3[,1:(col3-1)]
  
  fit <- glm(TRgoal3~.,data =TRattr3,family=binomial(link = "logit"))
  p3<-predict(fit,newdata = TEattr3 ,type="response")
  
  cla2<-ifelse(p3>=0.5,1,0)
  
  t3<-table(TEgoal3,cla2)
  
  acc3<-sum(diag(t3))/sum(t3)
  acc3<-format(round(acc3,2),nsmall=2)
  
  return(acc3)
  
}


##################################################################################

#savek<-0
#savek2<-0

#cross validation
KCROSS<-function(TR,foo)
{
  #same methodology as previous
  RR1<-dim(TR)[1]
  INTE<-as.integer(RR1/foo)
  
  #initialize accumulated validation accuracy
  ValidAcc1<-0
  ValidAcc2<-0
  
  #for final return
  Validcheck<-0
  traincheck<-0
  
  for(i in 0:(foo-1))
  {
    #start of validation data
    ST1<-1+i*INTE
    
    #end of validation data
    if(i!=(foo-1))
    {
      EN1<-(i+1)*INTE
    }else
    {
      EN1<-RR1
    }
    
    #get validation data
    VALID<-TR[ST1:EN1,]
    
    #get the split of training data
    if(i==0)
    {
      REALTR<-TR[(EN1+1):RR1,]  
    }else if(i==(foo-1))
    {
      REALTR<-TR[1:(ST1-1),]
    }else
    {
      DD1<-TR[1:(ST1-1),]
      DD2<-TR[(EN1+1):RR1,]
      REALTR<-rbind(DD1,DD2)
    }
    
    #use training data and validation data to get the results of each model 
    m1<-DTree(REALTR,VALID)
    #m2<-RF(REALTR,VALID)
    m2<-GLM(REALTR,VALID)
    
    #accumulate the validation accuracy
    ValidAcc1<-ValidAcc1+as.numeric(m1)
    ValidAcc2<-ValidAcc2+as.numeric(m2)
    
  }
  
  #calculate the average validation accuracy of every surogate model
  ValidAcc1<-ValidAcc1/foo
  ValidAcc2<-ValidAcc2/foo
  ValidAcc1<-format(round(ValidAcc1,2),nsmall=2)
  ValidAcc2<-format(round(ValidAcc2,2),nsmall=2)
  
  #which model
  
  choice<-0
  
  #select training and testing attribute and goal
  CC1<-dim(TR)[2]
  
  TRGG<-TR[,CC1]
  TRATT<-TR[,1:(CC1-1)]
  
  #assign("savek", TRGG, envir = .GlobalEnv)
  
  #judge which model is better
  if(ValidAcc1>ValidAcc2)
  {
    #fit the whole train data to model1
    fitT<- rpart(TRGG~., data =TRATT, method = 'class',model=TRUE,control=rpart.control(maxdepth=4))
    
    #get the predictive result of training data
    PTRAIN1<-predict(fitT,TRATT,type = "class")
   
    #table the result
    TTRAIN1<-table(TRGG,PTRAIN1)
   
    #calculate training accuracy
    ACCTRAIN1<-sum(diag(TTRAIN1))/sum(TTRAIN1)
    ACCTRAIN1<-format(round(ACCTRAIN1,2),nsmall=2)
    
    #confirm the return value
    traincheck<-ACCTRAIN1
    Validcheck<-ValidAcc1
    
    choice<-1
  }else
  {
    
    fit2T <- glm(TRGG~.,data =TRATT,family=binomial(link = "logit"))
    
    PTRAIN2<-predict(fit2T,newdata = TRATT ,type="response")
    
    claTR<-ifelse(PTRAIN2>=0.5,1,0)
    
    TTRAIN2<-table(TRGG,claTR)
    
    ACCTRAIN2<-sum(diag(TTRAIN2))/sum(TTRAIN2)
    ACCTRAIN2<-format(round(ACCTRAIN2,2),nsmall=2)
    
    traincheck<-ACCTRAIN2
    Validcheck<-ValidAcc2
    
    choice<-2
  }
  
  #list the result
  my_list<-list(traincheck,Validcheck,choice)
  
  return(my_list)
  
}


####################################################################################

#retrieve value from cross validation
ANS<-KCROSS(INPUT,Realfold)

#retrieve value
trAcc<-as.numeric(unlist(ANS[1]))
vaAcc<-as.numeric(unlist(ANS[2]))
chh<-as.numeric(unlist(ANS[3]))
#################################################################

##levels(TEST$Title) <- levels(df$Title)

coll<-dim(INPUT)[2]


fTRgoal<-INPUT[,coll]
fTRattr<-INPUT[,1:(coll-1)]

fTEattr<-Test[,1:(coll-1)]
fTEgoal<-Test[,coll]

finaltestacc<-0

prepp<-0

if(chh==1)
{
  
  Ffit<- rpart(fTRgoal~., data =fTRattr, method = 'class',model=TRUE,control=rpart.control(maxdepth=4))
  final2<-predict(Ffit,fTEattr,type = 'class')
  
  testtable2<-table(fTEgoal,final2)
  
  testacc2<-sum(diag(testtable2))/sum(testtable2)
  testacc2<-format(round(testacc2,2),nsmall=2)
  
  finaltestacc<-testacc2
  
  prepp<-final2
}else{
  
  Ffit2 <- glm(fTRgoal~.,data =fTRattr,family=binomial(link = "logit"))
  p2<-predict(Ffit2,newdata = fTEattr ,type="response")
  final<-ifelse(p2>=0.5,1,0)
  
  testtable<-table(fTEgoal,final)
  
  testacc<-sum(diag(testtable))/sum(testtable)
  testacc<-format(round(testacc,2),nsmall=2)
  
  finaltestacc<-testacc
  
  prepp<-final
}



#######################################################################

#accuracy output
newdf<-data.frame(accuracy="ave.",training=trAcc,validation=vaAcc,testing=finaltestacc)
write.csv(newdf,file=Re_out,row.names=FALSE)

#prediction result
#newdf2<-data.frame(Id=TESTID, Probability=final)
#names(newdf2)<-c("Id","Probability")
#write.csv(newdf2,file=PRE_out,row.names=FALSE)





