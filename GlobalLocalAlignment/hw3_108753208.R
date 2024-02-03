######################################
# the reference code of program3 
######################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Biostrings")
######################################
# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript hw3_<your student ID>.R --input test.fasta --score pam250.txt --aln global --gap -10 --output test_output.fasta", call.=FALSE)
}

# parse parameters
j<-grep("--",c(args[1:length(args)]))

inputcheck<-0
scorecheck<-0
alncheck<-0
gapcheck<-0
outputcheck<-0

if(length(j)!=5){
  stop("you should type 5 argument", call.=FALSE)
}else{
  
  if((args[j[1]]=="--input")|(args[j[2]]=="--input")|(args[j[3]]=="--input")|(args[j[4]]=="--input")|(args[j[5]]=="--input")){
    inputcheck<-1
  }else{
    stop("missing input", call.=FALSE)
  }
  
  if((args[j[1]]=="--score")|(args[j[2]]=="--score")|(args[j[3]]=="--score")|(args[j[4]]=="--score")|(args[j[5]]=="--score")){
    scorecheck<-1
  }else{
    stop("missing score", call.=FALSE)
  }
  
  if((args[j[1]]=="--aln")|(args[j[2]]=="--aln")|(args[j[3]]=="--aln")|(args[j[4]]=="--aln")|(args[j[5]]=="--aln")){
    alncheck<-1
  }else{
    stop("missing aln", call.=FALSE)
  }
  
  if((args[j[1]]=="--gap")|(args[j[2]]=="--gap")|(args[j[3]]=="--gap")|(args[j[4]]=="--gap")|(args[j[5]]=="--gap")){
    gapcheck<-1
  }else{
    stop("missing gap", call.=FALSE)
  }
  
  if((args[j[1]]=="--output")|(args[j[2]]=="--output")|(args[j[3]]=="--output")|(args[j[4]]=="--output")|(args[j[5]]=="--output")){
    outputcheck<-1
  }else{
    stop("missing output", call.=FALSE)
  }
  
}


for(i in 1:5){
  judge<-j[i]
  
  if(args[judge]=="--input"){
    fastafile<-args[judge+1]
  }else if(args[judge]=="--score"){
    pamscorefile<-args[judge+1]
  }else if(args[judge]=="--aln"){
    align<-args[judge+1]
  }else if(args[judge]=="--gap"){
    gap<-args[judge+1]
    gap<-as.numeric(gap)
  }else if(args[judge]=="--output"){
    outputname<-args[judge+1]
  }else{
    stop("Unknown flag", call.=FALSE)
  }
}

###############################################################

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



seqL1<-nchar(sequence[1])
seqL2<-nchar(sequence[2])

###############################################################################3


globalalign<-function()
{
  mm <- matrix(rep(0,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1))  #store score
  dir<-matrix(rep(0,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1)) #store direction
  #dir
  #1:up
  #2:left
  #3:diag
  
  for(i in 2:(seqL2+1))
  {
    mm[i,1]=mm[(i-1),1]+gap
    dir[i,1]=1
  }
  
  for(j in 2:(seqL1+1))
  {
    mm[1,j]=mm[1,(j-1)]+gap
    dir[1,j]=2
  }
  
  
  for(i in 2:(seqL2+1))
  {
    com2<-substring(sequence[2], i-1, i-1)
    for(j in 2:(seqL1+1))
    {
      com1<-substring(sequence[1], j-1, j-1)
      fromup<-mm[(i-1),j]+gap
      fromleft<-mm[i,(j-1)]+gap
      
      dia<-mm[(i-1),(j-1)]+M1[com2,com1]
      
      temp<-0
      tempdir<-0
      if(fromup>fromleft)
      {
        temp<-fromup
        tempdir<-1
      }else{
        temp<-fromleft
        tempdir<-2
      }
      
      if(temp>dia)
      {
        mm[i,j]<-temp
        dir[i,j]<-tempdir
      }else
      {
        mm[i,j]<-dia
        dir[i,j]<-3
      }
      
    }  
  }
  
  
  #paste(a,b,c)
  indx<-seqL1+1
  indy<-seqL2+1
  fin1<-''
  fin2<-''
  
  
  while (!((indx==1)&(indy==1)))
  {
    
    check<-dir[indy,indx]
    if(check==1)
    {
      aseq2<-substring(sequence[2], indy-1, indy-1)
      aseq1<-'-'  
      indy=indy-1
    }
    
    if(check==2)
    {
      aseq2<-'-'
      aseq1<-substring(sequence[1], indx-1, indx-1)
      indx=indx-1
    }
    
    if(check==3)
    {
      aseq2<-substring(sequence[2], indy-1, indy-1)
      aseq1<-substring(sequence[1], indx-1, indx-1)
      indy=indy-1
      indx=indx-1
    }
    
    fin1<-paste(aseq1,fin1,sep="")
    fin2<-paste(aseq2,fin2,sep="")
    
  }
  
  
  aa<-AAStringSet(fin1)
  names(aa)<-paste(seq_name[1],sep="")
  
  bb<-AAStringSet(fin2)
  names(bb)<-paste(seq_name[2],sep="")
  
  
  writeXStringSet(c(aa,bb),file=outputname,format="fasta")
  #print('1')
  
}

#######################################################################################

localalign<-function()
{
  bestx<-c()
  besty<-c()
  bestscore<-0
  
  mm <- matrix(rep(0,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1))
  dir<-matrix(rep(4,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1))
  #dir
  #1:up
  #2:left
  #3:diag
  #4:non
  
  for(i in 2:(seqL2+1))
  {
    mm[i,1]=0
    dir[i,1]=4
  }
  
  for(j in 2:(seqL1+1))
  {
    mm[1,j]=0
    dir[1,j]=4
  }
  
  
  for(i in 2:(seqL2+1))
  {
    
    com2<-substring(sequence[2], i-1, i-1)
    for(j in 2:(seqL1+1))
    {
      com1<-substring(sequence[1], j-1, j-1)
      
      fromup<-mm[(i-1),j]+gap
      fromleft<-mm[i,(j-1)]+gap
      
      dia<-mm[(i-1),(j-1)]+M1[com2,com1]
      
      temp<-0
      tempdir<-0
      if(fromup>fromleft)
      {
        temp<-fromup
        tempdir<-1
      }else{
        temp<-fromleft
        tempdir<-2
      }
      
      if(temp>dia)
      {
        if(temp>=0)
        {
          ######################3
          #when temp(from previous to now)==0
          #case1: above +, then minus10=>0   #record 4 or 1?1
          #case2: diagonal 0 +0   #record 4 or 3?3
          #case3: diagonal+ - #record 4 or 3?3
          ########################3
          mm[i,j]<-temp
          dir[i,j]<-tempdir
        }else
        {
            mm[i,j]<-0
            dir[i,j]<-4
        }
        
      }else
      {
        if(dia>=0)
        {
          mm[i,j]<-dia
          dir[i,j]<-3
        }else
        {
          mm[i,j]<-0
          dir[i,j]<-4
        }
        
      }
      
      
      if(mm[i,j]>bestscore)
      {
        
        bestx <- c(i)
        besty <- c(j)
        
        bestscore<-mm[i,j]
      }else{
        if(mm[i,j]==bestscore)
        {
          bestx <- append(bestx, i)
          besty <- append(besty, j)
        }
      }
      
      
      
    }  
  }
  
  
  temm1<-c()
  temm2<-c()
  
  
  for(i in 1:length(bestx))
  {
    nowx<-bestx[i]
    nowy<-besty[i]
    #paste(a,b,c)
    indx<-nowy #col
    indy<-nowx #row
    fin1<-''
    fin2<-''
    stopcheck<-FALSE
    
    while (!stopcheck)
    {
    
      check<-dir[indy,indx]
      if(check==1)
      {
        aseq2<-substring(sequence[2], indy-1, indy-1)
        aseq1<-'-'  
        indy=indy-1
      }
    
      if(check==2)
      {
        aseq2<-'-'
        aseq1<-substring(sequence[1], indx-1, indx-1)
        indx=indx-1
      }
    
      if(check==3)
      {
        aseq2<-substring(sequence[2], indy-1, indy-1)
        aseq1<-substring(sequence[1], indx-1, indx-1)
        indy=indy-1
        indx=indx-1
      }
    
      fin1<-paste(aseq1,fin1,sep='')
      fin2<-paste(aseq2,fin2,sep='')
      
      if(dir[indy,indx]==4)
      {
        stopcheck<-TRUE
      }
    
    }
    
    temm1<-append(temm1,fin1)
    temm2<-append(temm2,fin2)
    
  }
  
  
  bestL<-0
  bestind<-c()
  
  for(i in 1:length(temm1))
  {
  L1=nchar(temm1[i])
  
  if(L1>bestL)
  {
    bestind<-c(i)
    bestL<-L1
  }else
  {
    if(L1==bestL)
    {
      bestind<-append(bestind,i)
    }
  }
  }
  
  
  ans<-c()
  
  for(i in 1:length(bestind))
  {
    myind<-bestind[i]
    s1<-temm1[myind]
    
    aa<-AAStringSet(s1)
    myname<-paste('protein',toString(i),sep='')
    names(aa)<-paste(myname)
    
    ans<-append(ans,aa)
  }
  
  for(i in 1:length(bestind))
  {
    myind<-bestind[i]
    s2<-temm2[myind]
    
    bb<-AAStringSet(s2)
    myname<-paste('protein',toString(i),sep='')
    names(bb)<-paste(myname)
    
    ans<-append(ans,bb)
  }
  
  writeXStringSet(ans,file=outputname,format="fasta")
  #print('2')

}




########################################################
if(align=='global')
{
  globalalign()
}else{
  localalign()
}
##########################################################









