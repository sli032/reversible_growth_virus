# @File: genome.py
# @Author: Siyu Li (email:lisiyu211@gmail.com)
# @Date:   2018-03-20 17:13:12
# @Last Modified time: 2018-07-03 12:20:17
import hoomd
from hoomd import md
import numpy as np
import pandas as pd

class poly():
	"""docstring for polymer"""
	def __init__(self,param,cor,bonds):
		# super(polymer, self).__init__()
		self.n=param['len']
		self.a=2*param['Rg']
		self.rp=param['Rg']
		self.particles=cor
		self.typeid=0
		# self.pids=self.typeid*self.n
		self.bodyid=-1
		self.bonds=bonds
		self.bondid=0
	def gen_print(self):
		print("n:{},rp:{},particles:{},typeID:{},bodyID:{},bondID:{}".format(self.n,self.rp,self.particles,self.typeid,self.bodyid,self.bondid))
class sphere():
	"""docstring for polymer"""
	def __init__(self,param,center):
		# super(polymer, self).__init__()
		self.n=1
		self.rp=param['Rg']
		self.particles=center
		self.typeid=0
		# self.pids=self.typeid*self.n
		self.bodyid=-1
		self.bondid=0
	def copy_hoomd(self,sys):
		#this function is not tested yet
		self.particles=sys.particles[0].position
	def gen_print(self):
		print("n:{},rp:{},particles:{},typeID:{},bodyID:{},bondID:{}".format(self.n,self.rp,self.particles,self.typeid,self.bodyid,self.bondid))
#creat polymer
# n:polymer length(monomer number)
# a:kuhn length(bond length)
# r:bead radius
# typeid:polymer particle typeid, define for later potential pair
# bodyid:polymer particle bodyid, define for hoomd use
def create_polymer(param):
	# ,typeid_uni=True,bodyid_uni=True):
	"""docstring for Polymer"""
	n=param['len']
	row=param['row_n']
	cor=np.empty((n,3))
	start=np.array([-(row-1)*param['Rg'],-(param['turn_n']-1)*param['Rg'],param['Height']])
	step=np.array([2*param['Rg'],0.,0.])
	turnstep=np.array([0.,2*param['Rg'],0.])
	upstep=np.array([0.,0.,2*param['Rg']])
	cor[0]=start
	#define coordinate for polymer particles
	for i in range(1,n):
		if i % row==0:
			step=-step
			if (i/row) % param['turn_n']==0:
				cor[i]=cor[i-1]+upstep
				turnstep=-turnstep
			else:
				cor[i]=cor[i-1]+turnstep
		else: cor[i]=cor[i-1]+step
	#define bonds for polymer particles
	bonds=np.int_([[i,i+1] for i in range(n-1)])
	return poly(param,cor,bonds)
def create_sphere(param):
	center=np.array([[0.,0.,param['Rg']+param['Rs']]])
	return sphere(param,center)