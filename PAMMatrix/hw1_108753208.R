# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Rscript hw1_studentID.R --input input_path_mut.txt --pam x --output output_path_pamx.txt", call.=FALSE)
}

# parse parameters
j<-grep("-",c(args[1:length(args)]))

pamcheck<-0
inputcheck<-0
outputcheck<-0

if(length(j)!=3){
  stop("you should type 3 argument", call.=FALSE)
}else{
  
  if((args[j[1]]=="--pam")|(args[j[2]]=="--pam")|(args[j[3]]=="--pam")){
    pamcheck<-1
  }else{
    stop("missing pam", call.=FALSE)
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
  
  if(args[judge]=="--pam"){
    pamnum<-args[judge+1]
    pamnum<-as.numeric(pamnum)
  }else if(args[judge]=="--output"){
    out_f<-args[judge+1]
  }else if(args[judge]=="--input"){
    in_f<-args[judge+1]
  }else{
    stop("Unknown flag", call.=FALSE)
  }
}

###############################################################
#background probability
checkfi<-function(re)
{ 
  if(re=='G')
  {
    TEM<-0.089
  }else if(re=='A')
  {
    TEM<-0.087
  }else if(re=='L')
  {
    TEM<-0.085
  }else if(re=='K')
  {
    TEM<-0.081
  }else if(re=='S')
  {
    TEM<-0.070
  }else if(re=='V')
  {
    TEM<-0.065
  }else if(re=='T')
  {
    TEM<-0.058
  }else if(re=='P')
  {
    TEM<-0.051
  }else if(re=='E')
  {
    TEM<-0.050
  }else if(re=='D')
  {
    TEM<-0.047
  }else if(re=='R')
  {
    TEM<-0.041
  }else if(re=='N')
  {
    TEM<-0.040
  }else if(re=='F')
  {
    TEM<-0.040
  }else if(re=='Q')
  {
    TEM<-0.038
  }else if(re=='I')
  {
    TEM<-0.037
  }else if(re=='H')
  {
    TEM<-0.034
  }else if(re=='C')
  {
    TEM<-0.033
  }else if(re=='Y')
  {
    TEM<-0.030
  }else if(re=='M')
  {
    TEM<-0.015
  }else if(re=='W')
  {
    TEM<-0.010
  }else
  {
    TEM<-0
  }
  
 return(TEM)
}

###############################################################
#read PAM1 from data
M1<-read.table(in_f,sep='',skip=1,header=T,stringsAsFactors = FALSE)
M1<-as.matrix(M1)

#basic info for M1
row<-dim(M1)[1]
col<-dim(M1)[2]
rowname<-rownames(M1)
colname<-colnames(M1)

#change back the value #preprocess
for (i in 1:row)
{
  for(j in 1:col)
  {
    M1[i,j]=M1[i,j]/10000
  }  
}

#calculate M^n
Mn <-diag(nrow=row)

for (i in 1:pamnum)
{
  Mn=Mn%*%M1
}
###################################################
#to calculate PAM
PAM <- matrix(rep(0,row*col),nrow = row, ncol = col, dimnames = list(rowname,colname))
for (i in 1:row)
{
  rr=rowname[i]
  fi=checkfi(rr)
  for (j in 1:col)
  {
    PAM[i,j]=as.integer(10*log10((Mn[i,j]/fi)))
  }
}

######################################
##output PAM250 as a file
write.table(PAM, out_f, sep = " ",row.names = TRUE, col.names = TRUE, quote = FALSE)


