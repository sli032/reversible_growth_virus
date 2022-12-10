# from virus.capsid import get_all_grows,shell_grow,shell_detach,get_all_close
# import virus
import hoomd
import virus.capsid as cp
import numpy as np
from termcolor import colored,cprint
# import virus.param as param

ONEDGE=32767

def mc_grow(self):
	shell=self.shell
	shellE=self.shellE

	growlist=shell.get_all_grows(grow=1)
	if len(growlist)>0:
		[growi]=growlist.sample(weights=None,random_state=np.random.RandomState()).index

		cprint("growlist:\n{}\ngrowsite:{}".format(growlist,growi),"yellow")

		self.shell=shell.shell_grow(growlist,growi)
		self.relax()
		L,Ltwin=len(growlist),len(self.shell.detach_tlist())

		delE=self.shellE.loc[0,'Eshell']-shellE.loc[0,'Eshell']
		Pa=(L/Ltwin)/self.param['lambda']*np.exp((self.param['mu']-delE)/self.param['kT']) if Ltwin>0 else 1.
		cprint("oldT:{},newT:{},oldL:{},newL:{}\noldE:\n{}\nnewE:\n{}\ndelE:{},Pa{}".format(len(shell.triangle),len(self.shell.triangle),L,Ltwin,shellE,self.shellE,delE,Pa),'yellow')
		if np.random.rand()<min(1, Pa): 
			cprint("Grow accepted",'green')
			self.dump()
			return 3
		else:
			cprint("Grow rejected",'red')
			self.shell=shell
			self.shellE=shellE
			return 0
	else: cprint("no available grow",'red')
	return 0
		


def mc_detach(self):
	shell=self.shell
	shellE=self.shellE

	detach_tlist=shell.detach_tlist()
	cprint("detach_tlist:{}".format(detach_tlist),"yellow")
	if detach_tlist and len(shell.triangle)>1:
		ti=np.random.choice(detach_tlist)

		self.shell=shell.shell_detach(ti)
		self.relax()

		growlist=self.shell.get_all_grows(grow=1)
		L,Ltwin=len(detach_tlist),len(growlist)

		delE=self.shellE.loc[0,'Eshell']-shellE.loc[0,'Eshell']
		Pr=(L/Ltwin)*self.param['lambda']*np.exp((-self.param['mu']-delE)/self.param['kT']) if Ltwin>0 else 1.
		cprint("oldT:{},newT:{},oldL:{},newL:{}\noldE:\n{}\nnewE:\n{}\ndelE:{},Pr:{}".format(len(shell.triangle),len(self.shell.triangle),L,Ltwin,shellE,self.shellE,delE,Pr),'yellow')

		if np.random.rand()<min(1, Pr): 
			cprint("Detach accepted",'green')
			self.dump()
			return -1
		else:
			cprint("Detach rejected",'red')
			self.shell=shell
			self.shellE=shellE
			return 0
	else:
		cprint("no trimer could be detached!!",'red')
		return 0.

def mc_move(self):
	N_vmove,N_tmove,growvlist,diffusetlist=self.mc_info()
	cprint("MC moving.. choosing vmove or tmove:\nN_vmove:{}, N_tmove:{}".format(N_vmove,N_tmove),'yellow')
	convernum=0
	while (N_vmove+N_tmove)>0:
		shellE=self.shellE
		if np.random.rand()<N_vmove/(N_vmove+N_tmove):
			cprint("\nVertex Moving...",'yellow')
			movestatus=self.mc_vmove(growvlist)
		else:
			cprint("\nTrimer Moving...",'yellow')
			movestatus=self.mc_tmove(diffusetlist)
		self.printenergy(movestatus)
		if movestatus:
			N_vmove,N_tmove,growvlist,diffusetlist=self.mc_info()
			print("test E:",self.shellE.loc[0,'Eshell']," ",shellE.loc[0,'Eshell'],' ',self.shellE.loc[0,'Eshell']-shellE.loc[0,'Eshell'])
			if np.abs(self.shellE.loc[0,'Eshell']-shellE.loc[0,'Eshell'])<self.param['kT']:
					convernum+=1
			else: convernum=0
		else: convernum+=1
		self.time+=1
		print("converge steps:",convernum)

		if convernum>=self.param['convernum']: return 1
	cprint("no more availalbe mc move!",'red')
	return 0
	

def mc_vmove(self,growlist):
	shell=self.shell
	shellE=self.shellE

	if len(growlist)>0:
		[growi]=growlist.sample(weights=None,random_state=np.random.RandomState()).index

		cprint("growlist:\n{}\nmovesite:{}".format(growlist,growi),"yellow")

		self.shell=shell.shell_grow(growlist,growi)
		self.relax()

		delE=self.shellE.loc[0,'Eshell']-shellE.loc[0,'Eshell']
		Pm=np.exp(-delE/self.param['kT'])
		cprint("oldE:\n{}\nnewE:\n{}\ndelE:{},Pm:{}".format(shellE,self.shellE,delE,Pm),'yellow')
		if np.random.rand()<min(1, Pm): 
			cprint("Vertex move accepted",'green')
			self.dump()
			return 1
		else:
			cprint("Vertex move rejected",'red')
			self.shell=shell
			self.shellE=shellE
			return 0
	else: cprint("no available vertex move",'red')
	return 0
def mc_tmove(self,diffuse_tlist):
	shell=self.shell
	shellE=self.shellE

	if diffuse_tlist:
		ti=np.random.choice(diffuse_tlist)

		shellmiddle=shell.shell_detach(ti)

		growlist=shellmiddle.get_all_grows(grow=1)
		oldi=find_old_position(growlist,shell,ti)

		cprint("growlist:\n{}".format(growlist),"yellow")

		if len(growlist[growlist.index!=oldi])>0:
			[growi]=growlist[growlist.index!=oldi].sample(weights=None,random_state=np.random.RandomState()).index

			cprint("oldi:{},growi:{}".format(oldi,growi),'yellow')

			self.shell=shellmiddle.shell_grow(growlist,growi)
			self.relax()

			delE=self.shellE.loc[0,'Eshell']-shellE.loc[0,'Eshell']
			Pm=np.exp(-delE/self.param['kT'])
			cprint("oldE:\n{}\nnewE:\n{}\ndelE:{},Pm:{}".format(shellE,self.shellE,delE,Pm),'yellow')
			if np.random.rand()<min(1, Pm): 
				cprint("Trimer diffuse accepted",'green')
				self.dump()
				return 1
			else:
				cprint("Trimer diffuse rejected",'red')
				self.shell=shell
				self.shellE=shellE
				return 0
		else:
			cprint("no growsite available when moving triangle {}".format(ti),'red')
	else: cprint("no trimer could be moved",'red')
	return 0

def find_old_position(growlist,shell,ti):
	v0,v1,v2,l0,l1,l2,onedge=shell.triangle.loc[ti]
	l0p,l1p,l2p=np.abs(shell.line.loc[[l0,l1,l2],'edge'])
	if l0p+l1p+l2p<2*ONEDGE:
		for i,row in growlist[growlist['growtype']=='v'].iterrows():
			if not shell.vertex.loc[row['growindex'],'edge']==ONEDGE:
			# if and (gv in [v0,v1,v2]):
				return i
	else:
		for i,row in growlist[growlist['growtype']=='l'].iterrows():
			# if row['growindex'] in [l0p,l1p,l2p]:
			if not shell.line.loc[row['growindex'],'edge']==ONEDGE:
				return i
	return 32767
