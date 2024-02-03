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
  stop("USAGE: Rscript hw4_studentID.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}

# parse parameters
j<-grep("--",c(args[1:length(args)]))

inputcheck<-0
scorecheck<-0
alncheck<-0

opencheck<-0
openpenal<-0

extendcheck<-0
extendpenal<-0

outputcheck<-0
outputname<-''


if(length(j)!=6){
  stop("you should type 6 argument", call.=FALSE)
}else{
  
  if((args[j[1]]=="--input")|(args[j[2]]=="--input")|(args[j[3]]=="--input")|(args[j[4]]=="--input")|(args[j[5]]=="--input")|(args[j[6]]=="--input")){
    inputcheck<-1
  }else{
    stop("missing input", call.=FALSE)
  }
  
  if((args[j[1]]=="--score")|(args[j[2]]=="--score")|(args[j[3]]=="--score")|(args[j[4]]=="--score")|(args[j[5]]=="--score")|(args[j[6]]=="--score")){
    scorecheck<-1
  }else{
    stop("missing score", call.=FALSE)
  }
  
  if((args[j[1]]=="--aln")|(args[j[2]]=="--aln")|(args[j[3]]=="--aln")|(args[j[4]]=="--aln")|(args[j[5]]=="--aln")|(args[j[6]]=="--aln")){
    alncheck<-1
  }else{
    stop("missing aln", call.=FALSE)
  }
  
  if((args[j[1]]=="--gap_open")|(args[j[2]]=="--gap_open")|(args[j[3]]=="--gap_open")|(args[j[4]]=="--gap_open")|(args[j[5]]=="--gap_open")|(args[j[6]]=="--gap_open")){
    opencheck<-1
  }else{
    stop("missing gap_open", call.=FALSE)
  }
  
  if((args[j[1]]=="--gap_extend")|(args[j[2]]=="--gap_extend")|(args[j[3]]=="--gap_extend")|(args[j[4]]=="--gap_extend")|(args[j[5]]=="--gap_extend")|(args[j[6]]=="--gap_extend")){
    extendcheck<-1
  }else{
    stop("missing gap_extend", call.=FALSE)
  }
  
  if((args[j[1]]=="--output")|(args[j[2]]=="--output")|(args[j[3]]=="--output")|(args[j[4]]=="--output")|(args[j[5]]=="--output")|(args[j[6]]=="--output")){
    outputcheck<-1
  }else{
    stop("missing output", call.=FALSE)
  }
  
}


for(i in 1:6){
  judge<-j[i]
  
  if(args[judge]=="--input"){
    fastafile<-args[judge+1]
  }else if(args[judge]=="--score"){
    pamscorefile<-args[judge+1]
  }else if(args[judge]=="--aln"){
    align<-args[judge+1]
  }else if(args[judge]=="--gap_open"){
    gap<-args[judge+1]
    openpenal<-as.numeric(gap)
  }
  else if(args[judge]=="--gap_extend"){
    gap<-args[judge+1]
    extendpenal<-as.numeric(gap)
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
  #dir<-matrix(rep(0,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1)) #store direction
  IX <- matrix(rep(0,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1))  #for penalty of seq1
  IY <- matrix(rep(0,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1))  #for penalty of seq2
##########################################  

  for(j in 2:(seqL1+1))
  {
    IX[1,j]=-Inf
  }
  
  for(i in 2:(seqL2+1))
  {
    IX[i,1]=openpenal+(i-1)*extendpenal
    #dir[i,1]=1
  }
  ################################################
  for(j in 2:(seqL1+1))
  {
    IY[1,j]=openpenal+(j-1)*extendpenal
    #dir[1,j]=2
  } 
  
  for(i in 2:(seqL2+1))
  {
    IY[i,1]=-Inf
  }
  ####################################################
  for(j in 2:(seqL1+1))
  {
    mm[1,j]=-Inf
  }
  for(i in 2:(seqL2+1))
  {
    mm[i,1]=-Inf 
  }
  

  ############################################
  #dir
  #1:up
  #2:left
  #3:diag
  
  
  for(i in 2:(seqL2+1))
  {
    com2<-substring(sequence[2], i-1, i-1)
    for(j in 2:(seqL1+1))
    {
      com1<-substring(sequence[1], j-1, j-1)
      ##########################################
      dia1<-mm[(i-1),(j-1)]+M1[com2,com1]
      dia2<-IX[(i-1),(j-1)]+M1[com2,com1]
      dia3<-IY[(i-1),(j-1)]+M1[com2,com1]
      
      temp<-0
      tempdir<-0
      if(dia2>dia3)
      {
        temp<-dia2
        #tempdir<-1  #record which matrix you should go next, not this time
      }else{
        temp<-dia3
        #tempdir<-2
      }
      
      if(temp>=dia1)
      {
        mm[i,j]<-temp
        #dir[i,j]<-tempdir
      }else
      {
        mm[i,j]<-dia1
        #dir[i,j]<-3
      }
      
      ####################################
      fromup1<-mm[(i-1),j]+openpenal+extendpenal
      fromup2<-IX[(i-1),j]+extendpenal
      fromup3<-IY[(i-1),j]+openpenal+extendpenal
      
      te<-0
      if(fromup2>fromup3)
      {
        te<-fromup2  
      }else
      {
        te<-fromup3
      }
      
      if(te>=fromup1)
      {
        IX[i,j]<-te  
      }else{
        IX[i,j]<-fromup1
      }
      
      ######################################
      
      fromleft1<-mm[i,(j-1)]+openpenal+extendpenal
      fromleft2<-IY[i,(j-1)]+extendpenal
      fromleft3<-IX[i,(j-1)]+openpenal+extendpenal
      
      te2<-0
      
      if(fromleft2>fromleft3)
      {
        te2<-fromleft2  
      }else{
        te2<-fromleft3
      }
      
      if(te2>=fromleft1)
      {
        IY[i,j]<-te2
      }else{
        IY[i,j]<-fromleft1
      }
      ###########################################
      
    }  
  }
  
  
#############################################################  
  indx<-seqL1+1
  indy<-seqL2+1
  fin1<-''
  fin2<-''
  
  while (!((indx==1)&(indy==1)))
  {
    
    seex<-IX[indy,indx]
    seey<-IY[indy,indx]
    seem<-mm[indy,indx]
    
    newdir<-0
    temdir<-0
    temm<-0
    
    if(seex>seey)
    {
      temdir<-1
      temm<-seex
    }else{
      temdir<-2
      temm<-seey
    }
    
    if(temm>=seem)
    {
      newdir<-temdir
    }else{
      newdir<-3
    }
    
    if(newdir==1)
    {
      aseq2<-substring(sequence[2], indy-1, indy-1)
      aseq1<-'-'  
      indy=indy-1
    }
    
    if(newdir==2)
    {
      aseq2<-'-'
      aseq1<-substring(sequence[1], indx-1, indx-1)
      indx=indx-1
    }
    
    if(newdir==3)
    {
      aseq2<-substring(sequence[2], indy-1, indy-1)
      aseq1<-substring(sequence[1], indx-1, indx-1)
      indy=indy-1
      indx=indx-1
    }
    
    fin1<-paste(aseq1,fin1,sep="")
    fin2<-paste(aseq2,fin2,sep="")
    
  }
  
  #print(fin1)
  #print(fin2)
  
  
  
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
  
  mm <- matrix(rep(0,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1))  #store score
  dir<-matrix(rep(0,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1)) #store direction
  IX <- matrix(rep(0,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1))  #for penalty of seq1
  IY <- matrix(rep(0,(seqL1+1)*(seqL2+1)),nrow=(seqL2+1), ncol=(seqL1+1))  #for penalty of seq2
  #dir
  #1:up
  #2:left
  #3:diag
  #4:non
  
  dir[1,1]<-4
  ##################################################
  for(j in 2:(seqL1+1))
  {
    IX[1,j]=-Inf
  }
  
  for(i in 2:(seqL2+1))
  {
    IX[i,1]=openpenal+(i-1)*extendpenal
    #dir[i,1]=1
  }
  ################################################
  for(j in 2:(seqL1+1))
  {
    IY[1,j]=openpenal+(j-1)*extendpenal
    #dir[1,j]=2
  } 
  
  for(i in 2:(seqL2+1))
  {
    IY[i,1]=-Inf
  }
  ####################################################
  for(j in 2:(seqL1+1))
  {
    mm[1,j]=0 #-Inf
  }
  for(i in 2:(seqL2+1))
  {
    mm[i,1]=0 #-Inf 
  }
  
  #########################################################################
  
  tes1<-0
  tes2<-0
  
  
  for(i in 2:(seqL2+1))
  {
    
    com2<-substring(sequence[2], i-1, i-1)
    for(j in 2:(seqL1+1))
    {
      com1<-substring(sequence[1], j-1, j-1)
      ################################################
      
      dia1<-mm[(i-1),(j-1)]+M1[com2,com1]
      dia2<-IX[(i-1),(j-1)]+M1[com2,com1]
      dia3<-IY[(i-1),(j-1)]+M1[com2,com1]
      
      temp<-0
      #tempdir<-0
      if(dia2>dia3)
      {
        temp<-dia2
        #tempdir<-1  #record which matrix you should go next, not this time
      }else{
        temp<-dia3
        #tempdir<-2
      }
      
      if(temp>=dia1)
      { 
        if(temp>=0)
        {
          mm[i,j]<-temp
        }else{
          mm[i,j]<-0
          dir[i,j]<-4
        }
        #dir[i,j]<-tempdir
      }else
      {
        if(dia1>=0)
        {
          mm[i,j]<-dia1
        }else{
          mm[i,j]<-0
          dir[i,j]<-4
        }
        #dir[i,j]<-3
      }
      
      ############################################
      fromup1<-mm[(i-1),j]+openpenal+extendpenal
      fromup2<-IX[(i-1),j]+extendpenal
      fromup3<-IY[(i-1),j]+openpenal+extendpenal
      
      te<-0
      if(fromup2>fromup3)
      {
        te<-fromup2  
      }else
      {
        te<-fromup3
      }
      
      if(te>=fromup1)
      {
        IX[i,j]<-te  
      }else{
        IX[i,j]<-fromup1
      }
      
      
      ######################################
      
      fromleft1<-mm[i,(j-1)]+openpenal+extendpenal
      fromleft2<-IY[i,(j-1)]+extendpenal
      fromleft3<-IX[i,(j-1)]+openpenal+extendpenal
      
      te2<-0
      
      if(fromleft2>fromleft3)
      {
        te2<-fromleft2  
      }else{
        te2<-fromleft3
      }
      
      if(te2>=fromleft1)
      {
        IY[i,j]<-te2
      }else{
        IY[i,j]<-fromleft1
      }
      ###########################################
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
      
      #########################
      if(IX[i,j]>tes1)
      {
        tes1<-IX[i,j]
      }
      
      if(IY[i,j]>tes2){
        tes2<-IY[i,j]
      }
      ############################
      
    }  
  }
  
  #print(bestscore)
  #print(tes1)
  #print(tes2)
  
  
  temm1<-c()
  temm2<-c()
  
  for(i in 1:length(bestx))
  {
    nowx<-bestx[i]
    nowy<-besty[i]
    
    indx<-nowy #col
    indy<-nowx #row
    fin1<-''
    fin2<-''
    stopcheck<-FALSE
    #print(nowx)
    #print(nowy)
    #bestt<-0
    
    while (!stopcheck)
    {
      
      
      seex<-IX[indy,indx]
      seey<-IY[indy,indx]
      seem<-mm[indy,indx]
      
      newdir<-0
      temdir<-0
      temm<-0
      bestt<-0
      
      if(seex>seey)
      {
        temdir<-1
        temm<-seex
      }else{
        temdir<-2
        temm<-seey
      }
      
      if(temm>=seem)
      {
        newdir<-temdir
        bestt<-temm
      }else{
        newdir<-3
        bestt<-seem
      }
      ##################################
      ##make sense?
      ##if bestt=0 and newdir=1 or 2? stopcheck=False
      ##if bestt=0 and newdir=3? see if dir=4 then stopcheck=TRUE else FALSE
      if(bestt<=0)
      {
        stopcheck<-TRUE
        ###################
        #if(bestt==0)
        #{
        #  if((newdir==1)|(newdir==2))
        #  {
        #    stopcheck=FALSE
        #  }else{
        #    if(dir[i,j]!=4)
        #    {
        #      stopcheck=FALSE
        #    }
        #  }  
        #}
        #######################3
      }
      
      ######################################
      if(!stopcheck)
      {
        if(newdir==1)
        {
          aseq2<-substring(sequence[2], indy-1, indy-1)
          aseq1<-'-'  
          indy=indy-1
        }
        
        if(newdir==2)
        {
          aseq2<-'-'
          aseq1<-substring(sequence[1], indx-1, indx-1)
          indx=indx-1
        }
        
        if(newdir==3)
        {
          aseq2<-substring(sequence[2], indy-1, indy-1)
          aseq1<-substring(sequence[1], indx-1, indx-1)
          indy=indy-1
          indx=indx-1
        }
        
        fin1<-paste(aseq1,fin1,sep='')
        fin2<-paste(aseq2,fin2,sep='')
        
      }
      
      
      #######################################################
      
    }
    
    temm1<-append(temm1,fin1)
    temm2<-append(temm2,fin2)
  }
  
  
  #########################################
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
  #############################################
  
  
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
#seq1<-'LFVALYDFVASGDNTLSITKGEKLRVLGYNHNGE--WCEAQTKNGQGWVPSNYIT'
#seq2<-'VIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGYVPRNLLG'
#mk<-0
#lock<-0
#for(i in 1:nchar(seq1))
#{
#  com1<-substring(seq1, i, i)
#  com2<-substring(seq2, i, i)
#  if((com1!='-')&(com2!='-'))
#  {
#    mk<-mk+M1[com1,com2]
#    if(lock==1)
#    {
#      mk<-mk-(extendpenal)+openpenal
#      lock=0  
#    }
#  }else{
#    lock<-1
#    mk<-mk+(extendpenal)
#  }
  
#}
#print(mk)

#################################################################




