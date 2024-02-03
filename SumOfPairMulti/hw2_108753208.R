######################################
# the reference code of program2 
######################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Biostrings")
######################################
# Rscript hw2_studentID.R --input test.fasta --score pam250.txt --gopen -10 --gextend -2
#will generate the answer 1047
########################################
# initial
######################################
# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Rscript hw2_studentID.R --input test.fasta --score pam250.txt --gopen -10 --gextend -2", call.=FALSE)
}

# parse parameters
j<-grep("--",c(args[1:length(args)]))

inputcheck<-0
scorecheck<-0
gapopencheck<-0
gapextendcheck<-0

if(length(j)!=4){
  stop("you should type 4 argument", call.=FALSE)
}else{
  
  if((args[j[1]]=="--input")|(args[j[2]]=="--input")|(args[j[3]]=="--input")|(args[j[4]]=="--input")){
    inputcheck<-1
  }else{
    stop("missing input", call.=FALSE)
  }
  
  if((args[j[1]]=="--score")|(args[j[2]]=="--score")|(args[j[3]]=="--score")|(args[j[4]]=="--score")){
    scorecheck<-1
  }else{
    stop("missing score", call.=FALSE)
  }
  
  if((args[j[1]]=="--gopen")|(args[j[2]]=="--gopen")|(args[j[3]]=="--gopen")|(args[j[4]]=="--gopen")){
    gapopencheck<-1
  }else{
    stop("missing gopen", call.=FALSE)
  }
  
  if((args[j[1]]=="--gextend")|(args[j[2]]=="--gextend")|(args[j[3]]=="--gextend")|(args[j[4]]=="--gextend")){
    gapextendcheck<-1
  }else{
    stop("missing gextend", call.=FALSE)
  }
  
  
}



for(i in 1:4){
  judge<-j[i]
  
  if(args[judge]=="--input"){
    fastafile<-args[judge+1]
  }else if(args[judge]=="--score"){
    pamscorefile<-args[judge+1]
  }else if(args[judge]=="--gopen"){
    gonum<-args[judge+1]
    gonum<-as.numeric(gonum)
  }else if(args[judge]=="--gextend"){
    genum<-args[judge+1]
    genum<-as.numeric(genum)
  }else{
    stop("Unknown flag", call.=FALSE)
  }
}

###############################################################
#
M1<-read.table(pamscorefile,sep='',skip=9,header=T,stringsAsFactors = FALSE)
M1<-as.matrix(M1)
row<-dim(M1)[1]
col<-dim(M1)[2]

library("Biostrings",verbose=F,quietly=T)

# read fasta file
ff <- readAAStringSet(fastafile)
seq_name = names(ff)
sequence = paste(ff)
seq_num=length(seq_name)
###################################################################
total<-0
library(combinat)
if(seq_num!=2)
{
  d=combn(seq_num,2)
}else{
  d<-matrix(1:2, nrow = 2, ncol = 1)
}
w=dim(d)
s=w[2] #colnum

for (i in 1:s)
{
  q=d[,i]
  fir=q[1]
  fir=as.numeric(fir)
  
  sec=q[2]
  sec=as.numeric(sec)
  
  fir_seq=sequence[fir]
  sec_seq=sequence[sec]
  
  seq_len=nchar(fir_seq)
  ############################################
  temp_fir<-''
  temp_sec<-''
  for (n in 1:seq_len)
  {
    x1<-substring(fir_seq, n, n)
    x2<-substring(sec_seq, n, n)
    
    if((x1=='-') && (x2=='-'))
    {
    }else{
      temp_fir=paste(temp_fir,x1,sep="")
      temp_sec=paste(temp_sec,x2,sep="")
    }
  }
  fir_seq<-temp_fir
  sec_seq<-temp_sec
  seq_len=nchar(fir_seq)
  #print(temp_fir)
  #print(temp_sec)
  #################################################
  
  
  #firocheck<-0
  #secocheck<-0
  fir_buff<-c()
  sec_buff<-c()
  
  #temm1<-append(temm1,fin1)
  #
  
  for (k in 1: seq_len)
  {
    com1=substring(fir_seq, k, k)
    com2=substring(sec_seq, k, k)
    
    if ((com1!='-') && (com2!='-'))
    {
      total=total+M1[com1,com2]
    }else{
      if(com1=='-')
      {
        fir_buff<-append(fir_buff,k)
      }
      if(com2=='-')
      {
        sec_buff<-append(sec_buff,k)
      }
      
    }
    #if(com1!='-')
    #{
    #  firocheck<-0
    #}
    
    #if(com2!='-')
    #{
    #  secocheck<-0
    #}
    
    
    
    #if ((com1!='-') & (com2!='-'))
    #{
    #  total=total+M1[com1,com2]
    #}else
    #{
      ##########
      #only one of which will be -
      ########
    #  if((com1=='-'))
    #  {
    #    if(firocheck==0)
    #    {
    #      total=total+gonum
    #      firocheck<-1
    #    }else{
    #      total=total+genum
    #    }  
    #  }
      
    #  if((com2=='-'))
    #  {
    #    if(secocheck==0)
    #    {
    #      total=total+gonum
    #      secocheck<-1
    #    }else{
    #      total=total+genum
    #    }  
    #  }
        
      
      
      
    #}
    
  }
  
  ############################
  
  if(length(fir_buff)!=0)
  {
    con1<-TRUE
    last_ind<-fir_buff[1]
    total=total+gonum
  }
  
  if(length(fir_buff)>1)
  {
    for(m in 2:length(fir_buff))
    {
      now_ind<-fir_buff[m]
      if(now_ind==(last_ind+1))
      {
        con1<-TRUE
      }else{
        con1<-FALSE
      }
      last_ind<-now_ind
      if(con1)
      {
        total=total+genum
      }else{
        total=total+gonum
      }
    }
    
  }
  
  ########################
  if(length(sec_buff)!=0)
  {
    con2<-TRUE
    last_ind2<-sec_buff[1]
    total=total+gonum
  }
  
  if(length(sec_buff)>1)
  {
      
    for(m in 2:length(sec_buff))
    {
      now_ind2<-sec_buff[m]
      if(now_ind2==(last_ind2+1))
      {
        con2<-TRUE
      }else{
        con2<-FALSE
      }
      last_ind2<-now_ind2
      if(con2)
      {
        total=total+genum
      }else{
        total=total+gonum
      }
    }
  }
  
}

print(total)


