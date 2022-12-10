import sys,os,glob
from termcolor import cprint
from ast import literal_eval
import numpy as np
import pandas as pd
from io import StringIO
def child_analyze(subdir,file):
	i=0
	f=open(file,'r')
	line=f.readline()
	while line:
		if "frame:" in line:
			print(line,'time:',i)
			frame=line[:-5].split(':')[1]
		if 'shellE:' in line:
			i+=1
			if i==1:
				col=("frame subunit"+f.readline()).split()
				pdf=pd.DataFrame(columns=col)
			else: f.readline()
			Edata=f.readline()[2:]
			tnum=f.readline().split(':')[1]
			data=frame+' '+tnum+Edata
			df=",".join(col)+"\n"+",".join(data.split())
			pdf=pdf.append(pd.read_csv(StringIO(df)),ignore_index=True)
		line=f.readline()
	f.close()
	print(pdf)
	pdf.to_csv(subdir+'.csv')

if len(sys.argv)==1: childdir=os.listdir() 
else: childdir=sys.argv[1:]
for subdir in childdir:
	if os.path.isdir(subdir):
		os.chdir(subdir)
		# if glob.glob('final.gsd'):
		for resultfile in glob.glob('out'):
			cprint(subdir,'yellow')
			child_analyze(subdir,resultfile)
		os.chdir('..')