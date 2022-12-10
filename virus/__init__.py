# @File: __init__.py
# @Author: Siyu Li (email:lisiyu211@gmail.com)
# @Date:   2018-03-20 16:58:29
# @Last Modified time: 2019-02-26 10:36:04
from virus import capsid
from virus import genome
from virus import system
from virus import IOframe
from virus import snap
from virus import constituents
import numpy as np
import hoomd
from hoomd import md
import hoomd.deprecated as hoomdold
from termcolor import colored,cprint

# class virus(object):
# 	"""docstring for virus"""
# 	def __init__(self, param):
# 		self.param = param	
def run(param):
	cprint('Hoomd Run:\n','yellow')
	if param['METHOD']=='minimize':
		fire=hoomd.md.integrate.mode_minimize_fire(dt=param['dt'], ftol=param['ftol'], Etol=param['Etol'])
		nve=hoomd.md.integrate.nve(group=hoomd.group.nonrigid())
		count=0
		while not(fire.has_converged()) and count<param['maxrun']:
			print('couverge run:',count)
			hoomd.run(param['total'])
			count+=1
		nve.disable()
		hoomd.context.current.thermos[-1].disable()
		if count>=param['maxrun']:
			cprint('Converge Fail!','red')
			# param['ftol']=param['ftol']*100.
			# param['Etol']=param['Etol']*100.
			return 0
		return 1
	elif param['METHOD']=='lagenvin':
		hoomd.run(param['total'])
		return 2
	# hoomdold.dump.xml(filename="final%s.xml"%frame, group=hoomd.group.all(),all=True)
