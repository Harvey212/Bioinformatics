import argparse
import pandas as pd
import numpy as np
import math
from itertools import chain

#######################################################
parser = argparse.ArgumentParser(prog='hw5.py')
parser.add_argument("--input",required=True)#
parser.add_argument("--method",required=True)#
parser.add_argument("--output",required=True)#

args = parser.parse_args()
data = pd.read_csv(args.input,delim_whitespace=True,header=None,skiprows=1)
dd=np.array(data)

mymethod=args.method
outputfile=args.output
##########################################
def findMin(table):
	minn=math.inf
	minRow=-1
	minCol=-1
	
	############################################
	rowstart=0
	rowend=table.shape[0]-1

	for i in range(rowstart,rowend):
		colstart=i+1
		colend=table.shape[1]
		for j in range(colstart,colend):
			v=table[i][j]
			if v<minn:
				minn=v
				minRow=i
				minCol=j
	################################################


	return minRow,minCol




def upgma(table):
	row=table.shape[0]
	col=table.shape[1]

	labels=[]
	sets=[]

	temptable=table
	origintable=table

	for m in range(row):    #[0,1,2,3,4,5]
		labels.append(m)

	for k in range(row): #[[0],[1],[2],[3]]
		temp=[]
		temp.append(k)
		sets.append(temp)

	#######################################3
	buff=[]
	for n in range(row):
		buff.append(0)

	########################################
	#print(chr(65))
	base=65
	labelname=[]
	for mm in range(row):
		labelname.append(str(chr(base)))
		base+=1

	while len(sets)>1:
		#################
		rr,cc=findMin(temptable) ##rr: indexrow, cc:indexcol

		restindex=[]
		for k in range(len(sets)):
			if (k!=rr) and (k!=cc):
				restindex.append(k)
		####################################
		mydist=(temptable[rr][cc])/2
		tarA=mydist-buff[rr]
		tarA=f'{round(tarA,3):.3f}'#math.ceil(tarA * 1000)/1000.0#
		tarB=mydist-buff[cc]
		tarB=f'{round(tarB,3):.3f}'#math.ceil(tarB * 1000)/1000.0#

		tarRest=[]
		for mm in restindex:
			tarRest.append(buff[mm])

		###################################
		tempset=[]
		setA=sets[rr]
		setB=sets[cc]
		for a in setA:
			tempset.append(a)

		for b in setB:
			tempset.append(b)

		finset=[]
		finset.append(tempset)
		for c in restindex:
			finset.append(sets[c])

		sets=finset
		###########################################
		newtable=np.zeros((len(sets),len(sets)))

		Rst=0
		Rend=newtable.shape[0]-1

		for i in range(Rst,Rend):
			Cst=i+1
			Cend=newtable.shape[1]
			for j in range(Cst,Cend):
				sa=sets[i]
				sb=sets[j]
				totalcost=0

				for aa in sa:
					for bb in sb:
						cost=origintable[aa][bb]
						totalcost+=cost

				totalcost=totalcost/(len(sa)*len(sb))

				newtable[i][j]=totalcost
				newtable[j][i]=totalcost

		temptable=newtable

		############################
		fir=labels[rr]
		sec=labels[cc]

		com=[]
		com.append(fir)
		com.append(sec)

		fir1=labelname[rr]
		sec1=labelname[cc]
		comname='('+sec1+':'+str(tarB)+','+fir1+':'+str(tarA)+')'
		###############################
		buff=[]
		buff.append(mydist)#buff.append(math.ceil(mydist*1000)/1000.0)
		for kk in tarRest:
			buff.append(kk)#buff.append(math.ceil(kk*1000)/1000.0)#

		###################################3
		rest=[]
		restname=[]
		for j in restindex:
			rest.append(labels[j])
			restname.append(labelname[j])
		##############################################
		labels=[]
		labels.append(com)

		labelname=[]
		labelname.append(comname)

		for s in rest:
			labels.append(s)

		for ss in restname:
			labelname.append(ss)

		##################################3

	finans=labelname[0]+';'
	
	return finans

def createQ(table):
	row=table.shape[0]
	col=table.shape[1]
	Qtable=np.zeros((row,col))

	n=row

	for i in range(row):
		for j in range(col):
			if i!=j:
				summ1=0
				for k in range(n):
					summ1+=table[i][k]

				summ2=0
				for k in range(n):
					summ2+=table[j][k]

				Qtable[i][j]=(n-2)*table[i][j]-summ1-summ2

	return Qtable



def NJ(table):
	row=table.shape[0]
	col=table.shape[1]

	labels=[]
	sets=[]

	temptable=table

	for m in range(row):    #[0,1,2,3,4,5]
		labels.append(m)

	for k in range(row): #[[0],[1],[2],[3]]
		temp=[]
		temp.append(k)
		sets.append(temp)

	base=65
	labelname=[]
	for mm in range(row):
		labelname.append(str(chr(base)))
		base+=1

	#######################################
	#print(temptable)
	########################################

	while len(sets)>2:
		#################
		QQ=createQ(temptable)##rr: indexrow, cc:indexcol
		rr,cc=findMin(QQ)
		##########
		#if len(sets)==4:
		#	rr=1
		#	cc=2
		#print(QQ)
		#lala=[]
		#lala.append(rr)
		#lala.append(cc)
		#print(lala)
		#print(temptable)
		############
		
		restindex=[]
		for k in range(len(sets)):
			if (k!=rr) and (k!=cc):
				restindex.append(k)

		orisetnum=len(sets)
		###################################
		tempset=[]
		setA=sets[rr]
		setB=sets[cc]
		for a in setA:
			tempset.append(a)

		for b in setB:
			tempset.append(b)

		finset=[]
		finset.append(tempset)
		for c in restindex:
			finset.append(sets[c])

		sets=finset
		###########################################
		restdist=[]  #dist that is not related to the distance concerned
		relateddistA=[] #dist(a,k)
		relateddistB=[] #dist(b,k) 
		stardist=temptable[rr][cc] #dist(a,b)

		for i in range(0,temptable.shape[0]):
			tem=[]
			for j in range(0,temptable.shape[1]):
				if (i!=j):
					if (i!=rr) and (i!=cc) and (j!=rr) and (j!=cc):
						tem.append(temptable[i][j])

					if (i==rr) and (j!=cc) and (j!=rr):
						relateddistA.append(temptable[i][j])

					if (i==cc) and (j!=cc) and (j!=rr):
						relateddistB.append(temptable[i][j])

			if (i!=rr) and (i!=cc):
				restdist.append(tem)
		###########################################
		newtable=np.zeros((len(sets),len(sets)))

		for j in range(1,newtable.shape[1]):
			specialcost=0.5*(relateddistA[(j-1)]+relateddistB[(j-1)]-stardist)
			newtable[0][j]=specialcost
			newtable[j][0]=specialcost

		countrow=0
		for i in range(1,newtable.shape[0]):
			mycount=restdist[countrow]
			countcol=0
			for j in range(1,newtable.shape[1]):
				if i!=j:
					newtable[i][j]=mycount[countcol]
					countcol+=1

			countrow+=1

		temptable=newtable
		###################################
		buffk1=0
		for g in range(len(relateddistA)):
			buffk1+=relateddistA[g]

		buffk1+=stardist

		buffk2=0
		for g in range(len(relateddistB)):
			buffk2+=relateddistB[g]

		buffk2+=stardist

		buffA=0
		buffB=0


		buffA=0.5*stardist+(1/(2*(orisetnum-2)))*((buffk1)-(buffk2))
		buffB=stardist-buffA
		buffA=f'{round(buffA,3):.3f}'
		buffB=f'{round(buffB,3):.3f}'

		############################
		fir=labels[rr]
		sec=labels[cc]

		com=[]
		com.append(fir)
		com.append(sec)

		fir1=labelname[rr]
		sec1=labelname[cc]
		
		###################################3
		rest=[]
		restname=[]
		for j in restindex:
			rest.append(labels[j])
			restname.append(labelname[j])
		##############################################
		labels=[]
		labels.append(com)

		for s in rest:
			labels.append(s)
		#########################################
		if len(sets)!=2:
			comname='('+fir1+':'+str(buffA)+','+sec1+':'+str(buffB)+')'
			labelname=[]
			labelname.append(comname)

			for ss in restname:
				labelname.append(ss)
		else:
			lastname=restname[0]
			lastdist=temptable[0][1]
			lastdist=f'{round(lastdist,3):.3f}'
			ff='('+fir1+':'+str(buffA)+','+sec1+':'+str(buffB)+','+lastname+':'+str(lastdist)+')'
			labelname=[]
			labelname.append(ff)
	
	#print('NJ')
	#print(labelname[0])
	finans=labelname[0]+';'
	return finans



if mymethod=='UPGMA':
	myans=upgma(dd)
else:
	#NJ(dd)
	myans=NJ(dd)

with open(outputfile, 'w') as f:
	f.write(myans)
