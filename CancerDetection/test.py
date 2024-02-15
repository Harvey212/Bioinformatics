import pandas as pd
import numpy as np
import random
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn import svm
#from pyearth import Earth #for mar
from sklearn.naive_bayes import GaussianNB
from sklego.linear_model import LADRegression
#functional tree?
from pgmpy.models import BayesianNetwork
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import ElasticNet
from sklearn import linear_model
from sklearn.multiclass import OneVsRestClassifier

col1=['Age','hsa-miR-29a-3p','hsa-miR-30d-5p','hsa-miR-200a-3p','hsa-miR-200c-3p','hsa-miR-320d','hsa-miR-320c','hsa-miR-450b-5p','hsa-miR-203a','hsa-miR-486-3p','hsa-miR-1246','hsa-miR-1307-5p']
col2=['Age','hsa-miR-16-2-3p','hsa-miR-200a-3p','hsa-miR-200c-3p','hsa-miR-320b','hsa-miR-320d']
col3=['Age','hsa-miR-23b-3p','hsa-miR-29a-3p','hsa-miR-32-5p','hsa-miR-92a-3p','hsa-miR-150-5p','hsa-miR-200a-3p','hsa-miR-200c-3p','hsa-miR-203a','hsa-miR-320c','hsa-miR-320d','hsa-miR-335-5p','hsa-miR-450b-5p','hsa-miR-1246','hsa-miR-1307-5p']

def showResult(col):
	df = pd.read_csv('final_dset_combine.csv')
	#att=list(df.columns)
	##########################################
	X=df[col]
	y=df['outcome_cat']
	rows=X.shape[0]
	co=X.shape[1]
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.245, random_state=42)
	########################################
	##Linear discriminant analysis
	clf = LDA()
	y_pred =clf.fit(X_train, y_train).predict_proba(X_test)
	cc=roc_auc_score(y_test,y_pred,multi_class='ovo')
	cc=round(cc,2)
	print(f"AUC for Linear discriminant analysis is {cc}.")
	######################################
	##Logistic regression
	clf = LogisticRegression(random_state=1)
	y_pred =clf.fit(X_train, y_train).predict_proba(X_test)
	cc=roc_auc_score(y_test,y_pred,multi_class='ovo')
	cc=round(cc,2)
	print(f"AUC for Logistic regression is {cc}.")
	######################################
	##Neural network
	clf = MLPClassifier(solver='lbfgs', alpha=1e-5, random_state=1)
	y_pred =clf.fit(X_train, y_train).predict_proba(X_test)
	cc=roc_auc_score(y_test,y_pred,multi_class='ovo')
	cc=round(cc,2)
	print(f"AUC for Neural network is {cc}.")
	##################################
	##Support vector machine
	clf = svm.SVC(decision_function_shape='ovo',probability=True)
	y_pred =clf.fit(X_train, y_train).predict_proba(X_test)
	cc=roc_auc_score(y_test,y_pred,multi_class='ovo')
	cc=round(cc,2)
	print(f"AUC for Support vector machine is {cc}.")
	################################
	##Naive Bayes classifier
	gnb = GaussianNB()
	y_pred = gnb.fit(X_train, y_train).predict_proba(X_test)
	cc=roc_auc_score(y_test,y_pred,multi_class='ovo')
	cc=round(cc,2)
	print(f"AUC for Naive Bayes classifier is {cc}.")
	########################################
	##Random forest
	clf = RandomForestClassifier()
	y_pred=clf.fit(X_train, y_train).predict_proba(X_test)
	cc=roc_auc_score(y_test,y_pred,multi_class='ovo')
	cc=round(cc,2)
	print(f"AUC for Random forest is {cc}.")

print('Result for Significance-based selection:')
showResult(col1)
print('\n')
print('Result for Correlation-based feature subset selection:')
showResult(col2)
print('\n')
print('Result for Expression fold change selection:')
showResult(col3)





##failed model
##################################
#clf=LADRegression().fit(X_train, y_train)
#clf.fit(X_train, y_train)
#####################################
#clf = linear_model.MultiTaskElasticNet(alpha=0.1)
#y_pred = clf.fit(X_train, y_train).predict(X_test)
#print(y_pred)