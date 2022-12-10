import sys,os,glob
from termcolor import cprint
from ast import literal_eval
import numpy as np
from io import StringIO
import pandas as pd
import re
def child_analyze(subdir,fi):
	param=re.split("[(*)]",subdir)[:-1]
	param_name,param_value=param[::2],param[1::2]
	assert len(param_name)==len(param_value)
	for col,val in zip(param_name,param_value):
		if not col in df.columns:
			df[col] = pd.Series(np.zeros(len(df)), index=df.index)
		df.loc[fi,col]=val
	with open('shellinfo.txt','r') as f:
		for li in range(6):
			line=f.readline()
			tp,value=line.rstrip('\n').split(':')
			df.loc[fi,tp]=value

df=pd.DataFrame(columns=['R','gamma','Rp','shelltype','radius','pen_v','hex_v','tot_v','subunits'])
fi=0
if len(sys.argv)==1: childdir=os.listdir() 
else: childdir=sys.argv[1:]
for subdir in childdir:
	if os.path.isdir(subdir):
		os.chdir(subdir)
		if glob.glob('final.gsd'):
			cprint(subdir,'yellow')
			child_analyze(subdir,fi)
			fi+=1
		os.chdir('..')
print(df)
df.to_csv('Allshellinfo.csv')