# @File: IOframe.py
# @Author: Siyu Li (email:lisiyu211@gmail.com)
# @Date:   2018-03-20 16:57:48
# @Last Modified time: 2018-07-03 12:50:39
import numpy as np
import pandas as pd 
from ast import literal_eval


def frame_out(frame,sys):
	shell=sys.shell
	genome=sys.sys.particles[0]
	with open('shell.csv','a') as f:
		f.write("#frame:%d\n"%frame)
		f.write("#genome:{},{},{}\n".format(genome.position[0],genome.position[1],genome.position[2]))
		f.write("#vertex:%d\n"%len(shell.vertex))
		shell.vertex.to_csv(f)
		f.write("#line:%d\n"%len(shell.line))
		shell.line.to_csv(f)
		f.write("#triangle:%d\n"%len(shell.triangle))
		shell.triangle.to_csv(f)
def frame_in(param,shell,genome):
	with open('../{}.csv'.format(param['infile']),'r') as f:
		frame_row=0
		for row_num, line in enumerate(f, 1):
			if 'frame' in line:
				frame_row = row_num
				print('frame_row:',frame_row)
			if row_num>frame_row and line[0]=='#':
				data_row=row_num
				data_name,data_nrow=line.split(':')
				if data_name[1:]=='genome':
					genome.particles=np.array([[float(x) for x in data_nrow.split(',')]])
				else:
					print(data_name[1:],int(data_nrow))
					df=pd.read_csv('../{}.csv'.format(param['infile']), skiprows = data_row, index_col=0,header=0, nrows=int(data_nrow))
					if data_name[1:]=='vertex':
						for vi,row in df.iterrows():
							[x,y,z,l,t,edge]=row
							if l=='inf':
								shell.vertex.loc[vi]=[x,y,z,np.inf,np.inf,np.inf]
							else:
								shell.vertex.loc[vi]=[x,y,z,literal_eval(l),literal_eval(t),edge]
					elif data_name[1:]=='line':shell.line=df
					elif data_name[1:]=='triangle':shell.triangle=df
					else: raise NameError('Unknown input data')