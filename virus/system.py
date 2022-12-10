# @File: system.py
# @Author: Siyu Li (email:lisiyu211@gmail.com)
# @Date:   2018-03-20 16:57:51
# @Last Modified time: 2019-06-25 15:13:43
import hoomd
from hoomd import md
import virus
from . import _system
import numpy as np
import pandas as pd
from termcolor import colored,cprint
import sys
# import virus.param as param

ONEDGE=32767


def initial_system(snap,shell,param):
	sys=hoomd.init.read_snapshot(snap)
	setup_hoomd(param)

	logrec=np.array(['N','potential_energy','kinetic_energy','pair_lj_energy','pair_morse_energy','pair_yukawa_energy','bond_harmonic_energy'])
	logchoice=np.array([True,True,True,param['LJ'],param['MORSE'],param['YUKAWA'],param['BONDHARM']])
	log=hoomd.analyze.log(filename=None, quantities=logrec[logchoice], period=None)
	return virusystem(sys,log,shell,param)

class virusystem(object):
	"""docstring for system"""
	def __init__(self, sys,log,shell,param):
		self.sys = sys
		self.log = log
		self.shell=shell
		self.shellE=[]
		self.param=param
		self.time=0
		self.frame=0
	def dump(self,file=None):
		cprint("dump gsd, frame:{}".format(self.frame),'green',attrs=['reverse'])
		if file is None: file=self.param['gsdfile']
		hoomd.dump.gsd(filename=file,group=hoomd.group.all(),period=None,static=[])
		self.frame+=1
		# dump gsd each time after relaxation by setting perios=param['total'], otherwise set period=None and bring back another dump.gsd in the loop
		# hoomd.dump.gsd(filename="sphere.gsd",group=hoomd.group.all(),period=param['total'],static=[])
	def relax(self,printshell=False,printhoomd=True):
		# shell=self.shell
		cprint('Shell Hoomd:\n','yellow')

		snap=self.sys.take_snapshot(all=True)
		self.shell.update_shell_info()
		shstart,shend=virus.snap.assign_obj_to_snap(snap,self.shell,cover=True,dihedral=True)

		print('shstart,shend:',shstart,shend)
		if self.param['constituents']:
			cps=virus.constituents.create_constituents(self.shell,self.param['Rc'])
			cpstart,cpend=virus.snap.assign_obj_to_snap(snap,cps,cover=False,particle=True,bond=True,dihedral=False)
		if printhoomd: 
			cprint("snap before relax:\n",'yellow')
			virus.snap.snapprint(snap,particle=True,bond=True,dihedral=True)
		self.sys.restore_snapshot(snap)
		relaxstatus=virus.run(self.param)
		self.shell.copy_hoomd(self.sys,shstart)
		if printhoomd:
			snap=self.sys.take_snapshot(all=True)
			cprint("snap after relax:\n",'yellow')
			virus.snap.snapprint(snap,particle=True,bond=True,dihedral=True)
		if printshell: self.shell.shellprint()
		if not relaxstatus:
			self.dump("RelaxFail.gsd")
			self.shell.shellprint()
			virus.IOframe.frame_out(-1,self)
			cprint("Relax Fail!!! pay attention!",'red')
			# retry
			self.shell.update_shell_info()
			genome=virus.genome.create_sphere(self.param)
			snap=virus.snap.initial_snap(self.param['L'],genome,self.shell)
			self.sys=hoomd.init.read_snapshot(snap)
			setup_hoomd(self.param)
			logrec=np.array(['N','potential_energy','kinetic_energy','pair_lj_energy','pair_morse_energy','pair_yukawa_energy','bond_harmonic_energy'])
			logchoice=np.array([True,True,True,self.param['LJ'],self.param['MORSE'],self.param['YUKAWA'],self.param['BONDHARM']])
			self.log=hoomd.analyze.log(filename=None, quantities=logrec[logchoice], period=None)
			# snap=self.sys.take_snapshot(all=True)
			# shstart,shend=virus.snap.assign_obj_to_snap(snap,self.shell,cover=True,dihedral=True)
			# self.sys.restore_snapshot(snap)
			status=virus.run(self.param)
			self.shell.copy_hoomd(self.sys,shstart)
			virus.snap.snapprint(snap,particle=True,bond=True,dihedral=True)
			if not status:
				sys.exit("relax fail again!!!")
			cprint("new relax E:\n{}".format(self.energy()),'green')
			self.dump("RelaxFail.gsd")
		# self.shell=shell
		self.shellE=self.energy()

	def energy(self):
		log=self.log
		Ekey=self.param['growrec']
		shell=self.shell
		dic={
		'tot':'potential_energy',
		'k': 'kinetic_energy', 
		'lj': 'pair_lj_energy', 
		'morse': 'pair_morse_energy',
		'yukawa': 'pair_yukawa_energy', 
		'stretch': 'bond_harmonic_energy'}
		# Elist=[0.]*len(Ekey)
		Elist=pd.DataFrame(np.zeros((1,len(Ekey))),columns=Ekey)
		Estretch=log.query(dic['stretch'])
		Ebend=log.query(dic['tot'])-log.query(dic['lj'])-log.query(dic['morse'])-log.query(dic['yukawa'])-Estretch
		Eltension=shell.line_tension()
		Epent=shell.pentenergy()
		Ehp=shell.shell_Ehp()
		for i,key in enumerate(Ekey):
			if key=='tot':Elist[key]=log.query(dic['tot'])+Eltension+Epent+Ehp
			elif key=='Eshell':Elist[key]=Estretch+Ebend+Eltension+Epent+Ehp
			elif key=='Elasticity':Elist[key]=Estretch+Ebend
			elif key=='stretch':Elist[key]=Estretch
			elif key=='bend':Elist[key]=Ebend
			elif key=='hp':Elist[key]=Ehp
			elif key=='ltension':Elist[key]=Eltension
			elif key=='Epent':Elist[key]=Epent
			else: Elist[key]=log.query(dic[key])
		return Elist
	def printenergy(self,status):
		cprint('shellE:\n{}\nSubunits:{}\nEshell:{}\nEshell/Subunits:{}\nElasticity:{}\nElasticity/Subunits:{}'.format(self.shellE,len(self.shell.triangle),self.shellE.loc[0,'Eshell'],self.shellE.loc[0,'Eshell']/len(self.shell.triangle),self.shellE.loc[0,'Elasticity'],self.shellE.loc[0,'Elasticity']/len(self.shell.triangle)),'green',attrs=['reverse'] if status else [])
		if status: self.shell.shellprint()

	def mc_grow(self):
		return _system.mc_grow(self)
	def mc_detach(self):
		return _system.mc_detach(self)
	def mc_info(self):
		shell=self.shell
		growlist=shell.get_all_grows(grow=2)
		diffusetlist=shell.diffuse_tlist()
		cprint("\n\n\nMC system growlist:\n{}\n\diffusetlist:{}".format(growlist,diffusetlist),'green')
		growtlist=growlist[growlist['choice']=='G']
		growvlist=growlist[growlist['choice']=='S']
		N_vmove,N_tmove=len(growvlist),len(diffusetlist)*len(growtlist)
		# N_tmove=len(diffusetlist)*(len(growtlist)-2) 
		#estimate growtlist without filter
		return N_vmove,N_tmove,growvlist,diffusetlist
	def mc_move(self):
		return _system.mc_move(self)
	def mc_vmove(self,growlist):
		return _system.mc_vmove(self,growlist)
	def mc_tmove(self,detach_tlist):
		return _system.mc_tmove(self,detach_tlist)
	
def setup_hoomd(param):
	nl=hoomd.md.nlist.cell()
	if param['LJ']:
		lj = hoomd.md.pair.lj(0,nlist=nl)
		lj.set_params(mode='shift')
		#Repulsive LJ part, 'G' with any types other than itself
		if param['Rg']==0.:lj.pair_coeff.set('G', 'S', epsilon=0., sigma=0., alpha=0.)
		else: lj.pair_coeff.set('G', 'S', epsilon=param['lj_eps'], sigma=param['Rg']+param['Rs'], alpha=param['lj_alpha'], r_cut=param['lj_rcut'])
		#Standard LJ, 'G' with itself
		lj.pair_coeff.set('G', ['G','C'], epsilon=0., sigma=0., alpha=0.)
		lj.pair_coeff.set('C', 'C', epsilon=param['lj_eps'], sigma=2*param['Rc'], alpha=param['lj_alpha'], r_cut=2*param['Rc'])
		lj.pair_coeff.set('C', 'S', epsilon=param['lj_eps'], sigma=param['Rc']+param['Rs'], alpha=param['lj_alpha'], r_cut=param['Rc']+param['Rs'])
		# No Interaction
		lj.pair_coeff.set('S', 'S', epsilon=0., sigma=0., alpha=0.)
	if param['MORSE']:
		morse=hoomd.md.pair.morse(r_cut=param['morse_rcut'],nlist=nl)
		morse.pair_coeff.set('S', 'S', D0=param['morse_D0'], alpha=param['morse_alpha'], r0=param['morse_R0'],r_on=0.8)
		morse.pair_coeff.set('G', ['S','C'], D0=0., alpha=0., r0=0.)
		morse.pair_coeff.set('G', 'G', D0=0., alpha=0., r0=0.)
		morse.pair_coeff.set('C', ['S','G','C'], D0=0.,alpha=0.,r0=0.)
	if param['YUKAWA']:
		yukawa=hoomd.md.pair.yukawa(0,nlist=nl)
		yukawa.set_params(mode='xplor')
		# Attractive
		yukawa.pair_coeff.set('G','S', epsilon=-2.*param['yukawa_eps'], kappa=1., r_cut=param['yukawa_rcut'], r_on=param['yukawa_ron'])
		# Repulsive
		yukawa.pair_coeff.set('G','G', epsilon=param['yukawa_eps'], kappa=1., r_cut=param['yukawa_rcut'], r_on=param['yukawa_ron'])
		# No Interaction
		yukawa.pair_coeff.set('S','S', epsilon=0., kappa=0., r_cut=0., r_on=0.)
		yukawa.pair_coeff.set('C', ['C','G','S'], llepsilon=0., kappa=0., r_cut=0.,r_on=0.)
	print('k_s:',param['sh_bondk'],',k_genome:',param['gen_bondk'],'k_shcp:',param['shcp_bondk'],',k_b:',param['sh_dihedk'],'theta0:',param['theta0'])
	if param['BONDHARM']:
		harmonic = hoomd.md.bond.harmonic()
		harmonic.bond_coeff.set('shell', k=param['sh_bondk'], r0=param['l0'])
		harmonic.bond_coeff.set('genome', k=param['gen_bondk'], r0=2*param['Rg'])
		harmonic.bond_coeff.set('shcp',k=param['shcp_bondk'],r0=param['l0']*np.sqrt(3)/3)
		# btable = bond.table(width=1000)
		# btable.bond_coeff.set('gensh', func=bond_harmonic, rmin=0., rmax=Rp+Rx, coeff=dict(kappa=330, r0=Rp+Rx))
	if param['DIHEDHARM']:
		dtable = hoomd.md.dihedral.table(width=10000)
		dtable.dihedral_coeff.set('capsomer', func=dihedral_harmonic, coeff=dict(kappa=param['sh_dihedk'], theta0=np.pi-param['theta0']))
	print("forces:",hoomd.context.current.forces)
	print(hoomd.context.current.forces[0].bond_coeff.values)
	print(hoomd.context.current.forces[1].dihedral_coeff.values)
	print("sys:",hoomd.context.current.system)

def bond_harmonic(r, rmin, rmax, kappa, r0):
   V = 0.5 * kappa * (r-r0)**2;
   F = -kappa*(r-r0);
   return (V, F)
def dihedral_harmonic(theta, kappa, theta0):
	# print('theta:',theta)
	V = kappa * (1 - np.cos(theta-theta0));
	F = - kappa * np.sin(theta-theta0)
	return (V, F)