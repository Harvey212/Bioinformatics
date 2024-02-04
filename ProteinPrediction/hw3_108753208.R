library('rpart')
library(caret)
library(party)


# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript hw3_studentID.R --fold k --input Archaeal_tfpssm.csv --output performance.csv", call.=FALSE)
}

# parse parameters
j<-grep("-",c(args[1:length(args)]))

foldcheck<-0
inputcheck<-0
outputcheck<-0

if(length(j)!=3){
  stop("you should type 3 argument", call.=FALSE)
}else{
  
  if((args[j[1]]=="--fold")|(args[j[2]]=="--fold")|(args[j[3]]=="--fold")){
    foldcheck<-1
  }else{
    stop("missing fold command", call.=FALSE)
  }
  
  if((args[j[1]]=="--input")|(args[j[2]]=="--input")|(args[j[3]]=="--input")){
    inputcheck<-1
  }else{
    stop("missing input", call.=FALSE)
  }
  
  if((args[j[1]]=="--output")|(args[j[2]]=="--output")|(args[j[3]]=="--output")){
    outputcheck<-1
  }else{
    stop("missing output", call.=FALSE)
  }
  
}


for(i in 1:3){
  judge<-j[i]
  
  if(args[judge]=="--fold"){
      FOLD<-args[judge+1]
  }else if(args[judge]=="--output"){
    out_f<-args[judge+1]
  }else if(args[judge]=="--input"){
    df<-args[judge+1]
  }else{
    stop("Unknown flag", call.=FALSE)
  }
  
}

#model1
DTree<-function(TRData,VAData)
{
  fit<- rpart(TRData[,1]~., data =TRData, method = 'class',model=TRUE,control=rpart.control(maxdepth=2))
  p1<-predict(fit,VAData,type = "class")
  t1<-table(VAData[,1],p1)
  acc1<-sum(diag(t1))/sum(t1)
  acc1<-format(round(acc1,2),nsmall=2)
  
  return(acc1)
}

#model2
DTree2<-function(TRData2,VAData2)
{
  fit2<- rpart(TRData2[,1]~., data =TRData2, method = 'class',model=TRUE,control=rpart.control(maxdepth=1))
  p2<-predict(fit2,VAData2,type = "class")
  t2<-table(VAData2[,1],p2)
  acc2<-sum(diag(t2))/sum(t2)
  acc2<-format(round(acc2,2),nsmall=2)
  
  return(acc2)
}


#RF<-function(TRData1,VAData1)
#{
  #cf1<-cforest(TRData1[,1]~., data=TRData1, control=cforest_unbiased(mtry=NULL,ntree=5))
  #p2<- predict(cf1, newdata =VAData1, type = "response")
  #t2<-table(VAData1[,1],p2)
  #acc2<-sum(diag(t2))/sum(t2)
  #acc2<-format(round(acc2,2),nsmall=2)
  
  #return(acc2)
#}


#cross validation
KCROSS<-function(TR,TE,foo)
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
  testcheck<-0
  
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
    m2<-DTree2(REALTR,VALID)
    
    #accumulate the validation accuracy
    ValidAcc1<-ValidAcc1+as.numeric(m1)
    ValidAcc2<-ValidAcc2+as.numeric(m2)
    
  }
  
  #calculate the average validation accuracy of every surogate model
  ValidAcc1<-ValidAcc1/foo
  ValidAcc2<-ValidAcc2/foo
  ValidAcc1<-format(round(ValidAcc1,2),nsmall=2)
  ValidAcc2<-format(round(ValidAcc2,2),nsmall=2)
  
  #judge which model is better
  if(ValidAcc1>ValidAcc2)
  {
    #fit the whole train data to model1
    fit<- rpart(TR[,1]~., data =TR, method = 'class',model=TRUE,control=rpart.control(maxdepth=2))
    
    #get the predictive result of training data
    PTRAIN1<-predict(fit,TR,type = "class")
    #get the predictive result of testing data
    PTEST1<-predict(fit,TE,type = "class")
    
    #table the result
    TTRAIN1<-table(TR[,1],PTRAIN1)
    TTEST1<-table(TE[,1],PTEST1)
    
    #calculate training accuracy
    ACCTRAIN1<-sum(diag(TTRAIN1))/sum(TTRAIN1)
    ACCTRAIN1<-format(round(ACCTRAIN1,2),nsmall=2)
    
    #calculate testing accuracy
    ACCTEST1<-sum(diag(TTEST1))/sum(TTEST1)
    ACCTEST1<-format(round(ACCTEST1,2),nsmall=2)
    
    #confirm the return value
    traincheck<-ACCTRAIN1
    testcheck<-ACCTEST1
    Validcheck<-ValidAcc1
  }else
  {
    #cf1<-cforest(TR[,1]~., data=TR, control=cforest_unbiased(mtry=NULL,ntree=5))
    fit2<- rpart(TR[,1]~., data =TR, method = 'class',model=TRUE,control=rpart.control(maxdepth=1))
    
    #PTRAIN2<- predict(cf1, newdata =TR, type = "response")
    #PTEST2<- predict(cf1, newdata =TE, type = "response")
    PTRAIN2<-predict(fit2,TR,type = "class")
    PTEST2<-predict(fit2,TE,type = "class")
    
    TTRAIN2<-table(TR[,1],PTRAIN2)
    TTEST2<-table(TE[,1],PTEST2)
    
    ACCTRAIN2<-sum(diag(TTRAIN2))/sum(TTRAIN2)
    ACCTRAIN2<-format(round(ACCTRAIN2,2),nsmall=2)
    
    ACCTEST2<-sum(diag(TTEST2))/sum(TTEST2)
    ACCTEST2<-format(round(ACCTEST2,2),nsmall=2)
    
    traincheck<-ACCTRAIN2
    testcheck<-ACCTEST2
    Validcheck<-ValidAcc2
  }
  
  #list the result
  my_list<-list(traincheck,testcheck,Validcheck)
  
  return(my_list)
  
}



INPUT <-read.csv(df,header=FALSE)
FOLD<-as.numeric(FOLD)

row<-dim(INPUT)[1]
col<-dim(INPUT)[2]

colone<-INPUT[1:row,1]
colrest<-INPUT[1:row,3:col]

#colmerge<-cbind(colone,colrest)
#NEWd<-cbind(INPUT$V2,colmerge)
NEWd<-cbind(INPUT$V2,colrest)

#shuffle the data
set.seed(9850)
g<-runif(nrow(NEWd))
INPUTr<-NEWd[order(g),]

#total row
RR<-dim(INPUTr)[1]

#interval of each fold
interval<-as.integer(RR/FOLD)

#initialize lists of each desired column
SET<-c()
TTRR<-c()
VA<-c()
TTEE<-c()

#initialize average value
TRave<-0
VLave<-0
TEave<-0


for(i in 0:(FOLD-1))
{
  #get start point of test data
  start<-1+i*interval
  
  #get end point of test data
  if(i!=(FOLD-1))
  {
    end<-(i+1)*interval
  }else
  {
    end<-RR
  }
  
  #cat("(",start,",",end,")")
  #get test data
  TEST<-INPUTr[start:end,]
  #get train data
  if(i==0)
  {
    TRAIN<-INPUTr[(end+1):RR,]
  }else if(i==(FOLD-1))
  {
    TRAIN<-INPUTr[1:(start-1),]
  }else
  {
    D1<-INPUTr[1:(start-1),]
    D2<-INPUTr[(end+1):RR,]
    TRAIN<-rbind(D1,D2)
  }
  
  TRVAfold<-FOLD-1
  #retrieve value from cross validation
  ANS<-KCROSS(TRAIN,TEST,TRVAfold)
  
  INDEX<-paste0("fold",(i+1))
  
  #retrieve value
  first<-as.numeric(unlist(ANS[1]))
  second<-as.numeric(unlist(ANS[2]))
  third<-as.numeric(unlist(ANS[3]))
  
  SET<-c(SET,INDEX)
  TTRR<-c(TTRR,first)
  TTEE<-c(TTEE,second)
  VA<-c(VA,third)
  #accumulate average
  TRave<-TRave+first
  TEave<-TEave+second
  VLave<-VLave+third
  
}


#unlist the list
TTRR<-as.numeric(unlist(TTRR))
TTEE<-as.numeric(unlist(TTEE))
VA<-as.numeric(unlist(VA))

#calculate average
TRave<-TRave/FOLD
VLave<-VLave/FOLD
TEave<-TEave/FOLD

#round average
TRave<-format(round(TRave,2),nsmall=2)
VLave<-format(round(VLave,2),nsmall=2)
TEave<-format(round(TEave,2),nsmall=2)

#combine data frame
out1<-data.frame(set=SET,training=TTRR,validation=VA,test=TTEE)
out2<-data.frame("ave.",TRave,VLave,TEave)
names(out2)<-c("set","training","validation","test")
newdf<-rbind(out1,out2)

write.csv(newdf,file=out_f,row.names=FALSE)




