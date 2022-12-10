# @File: run.py
# @Author: Siyu Li (email:lisiyu211@gmail.com)
# @Date:   2018-03-20 18:08:51
# @Last Modified time: 2019-07-15 17:25:54
import hoomd
from hoomd import md
import hoomd.deprecated as hoomdold
import numpy as np
import pandas as pd
import json
import virus
import sys,os,glob
from termcolor import colored,cprint

ONEDGE=32767
Succ=1
Fail=0

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open(param['outfile'], "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass 

def prep_run():
	# initial shell
	shell=virus.capsid.create_shell(param)
	if param['infile']=='NA':
		genome=virus.genome.create_sphere(param)
		shell.initial()
		# virus.capsid.shell_grow(shell)
		# shell.attach(1)
		# shell.attach(2)
		# shell.attach(3)
	else: 
		print('Import %s.csv'%param['infile'])
		genome=virus.genome.create_sphere(param)
		virus.IOframe.frame_in(param,shell,genome)
	shell.update_shell_info()
	shell.shellprint()

	if param['constituents']:
		cps=virus.constituents.create_constituents(shell,param['Rc'])
		snap=virus.snap.initial_snap(param['L'],genome,shell,cps)
	else:
		snap=virus.snap.initial_snap(param['L'],genome,shell)
	# initial hoomdsnap
	virus.snap.snapprint(snap)
	
	return virus.system.initial_system(snap,shell,param)

def run():
	shell_status=0,0
	action=False
	virusys=prep_run()
	virusys.dump()
	virusys.relax(printhoomd=True)

	print("initial shell:")
	virusys.printenergy(1)

	while (not virusys.shell.test_shell_close()) & (virusys.time<param['maxtime']):
		
		# Monte Carlo Move
		virusys.mc_move()
		mergestatus=test_merge(virusys)

		# Grand Canonical MC move
		if mergestatus!=-1:
			if np.random.rand()<0.5:
				cprint("\nTime:{}\nGrowing...".format(virusys.time),'green',attrs=['reverse'])
				action=virusys.mc_grow()
				virusys.printenergy(action)
			else:
				cprint("\nTime:{}\nDetaching...".format(virusys.time),'red',attrs=['reverse'])
				action=virusys.mc_detach()
				virusys.printenergy(action)
		mergestatus=test_merge(virusys)
		# if action: update_sys(virusys)
		virus.IOframe.frame_out(virusys.time,virusys)
		

	virusys.dump("final2.gsd")
	virus.IOframe.frame_out(virusys.time,virusys)
	cprint("Assemble Done!\nFinal Shell:",'white',attrs=['reverse'])
	virusys.printenergy(1)

def update_sys(sys):
	# if len(sys.shell.triangle) % 10==0: 
	# 	sys.param['mu']+=param['delmu']
	# 	cprint("new mu:{}".format(sys.param['mu']),'yellow')
	if len(sys.shell.triangle) % int(np.round(0.8*param['idealT']))==0:
		print('paramkT:',param['kT'])
		sys.param['kT']=0.1*param['kT']
		cprint("Annealing trimer num:{} new temp: {}".format(int(np.round(0.8*param['idealT'])),sys.param['kT']),'yellow')
		# sys.param['mu']+=param['delmu']
		# cprint("new mu:{}".format(sys.param['mu']),'yellow')
# def mcrun(sys):
	# status: 1: tmove, 2:vmove, 3:grow, -1:detach, 0:nothing
	# trialn,moven,movestatus=0,0,0
	# action=False

	# # Monte Carlo Move
	# sys.mc_move()
	# mergestatus=test_merge(sys)

	# if mergestatus!=-1:
	# 	if np.random.rand()<0.5:
	# 		cprint("\nTime:{}\nGrowing...".format(sys.time),'green',attrs=['reverse'])
	# 		action=sys.mc_grow()
	# 		sys.printenergy(action)
	# 	else:
	# 		cprint("\nTime:{}\nDetaching...".format(sys.time),'red',attrs=['reverse'])
	# 		action=sys.mc_detach()
	# 		sys.printenergy(action)
	# mergestatus=test_merge(sys)

	# Grand Canonical MC move
	# if mergestatus!=-1:
	# 	cprint("GCMC moving")
	# 	action,initialtime=False,sys.time
	# 	while (not action) and (sys.time-initialtime<20): 
	# 		if np.random.rand()<0.5:
	# 			cprint("\nTime:{}\nGrowing...".format(sys.time),'green',attrs=['reverse'])
	# 			if sys.mc_grow(): action=True
	# 			printenergy(sys.shell,sys.shellE,action)
	# 		else:
	# 			cprint("\nTime:{}\nDetaching...".format(sys.time),'red',attrs=['reverse'])
	# 			if sys.mc_detach(): action=True
	# 			printenergy(sys.shell,sys.shellE,action)
	# 		sys.time+=1
		
		
def test_merge(sys):
	cprint("test merging...",'yellow')
	mergestatus=sys.shell.shell_merge()
	while mergestatus==1:
		cprint('Force Merge Happen','yellow')
		sys.relax()
		sys.printenergy(1)
		sys.dump()
		mergestatus=sys.shell.shell_merge()
	if mergestatus==-1:
		cprint("Assemble done from merge")
		sys.relax()
		sys.dump("final1.gsd")
	else: cprint("no more merge!",'yellow')
	return mergestatus
		
param=sys.argv[1]
cprint(param,'green')
with open(param) as f: param = json.load(f)
if param['outfile']!='NA':
	sys.stdout=Logger()
run()