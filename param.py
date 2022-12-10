#!/usr/bin/python2
# @File: param.py
# @Author: Siyu Li (email:lisiyu211@gmail.com)
# @Date:   2018-03-20 15:05:37
# @Last Modified time: 2019-03-14 14:06:43

# import itertools
import json
# import math
# import pathlib
# import string
import os.path
import numpy as np
import itertools
import yaml
import shutil


def give_child_name(dic,namekey):
	name=''
	for key in namekey:
		name+='%s(%s)'%(key,dic[key])
	return name
def generate_child_param(name,dic,path):
	with open(path+'/%s.json'%name, 'w') as file:
		json.dump(dic,file,indent=4, sort_keys=True, separators=(',', ': '))
def default_param(dic,key,value):
	if not dic.has_key(key):
		dic[key]=value
def train_child_param(dic):	
	#Genome section:
	default_param(dic,'row_n',     round(dic['len']**(1./3)))
	default_param(dic,'turn_n',    round(dic['len']**(1./3)))
	default_param(dic,'Height',    dic['R']-np.sqrt(abs((dic['R']-dic['Rs'])**2 - dic['Rg']**2*(dic['turn_n']**2 + dic['row_n']**2))))
	default_param(dic,'Rg',        0.)
	

	#Shell section:
	default_param(dic,'gamma', 1.)
	default_param(dic,'sh_dihedk', 1.)
	default_param(dic,'sh_bondk',  dic['sh_dihedk']*dic['gamma'])
	default_param(dic,'hp_eps',    0.001*dic['sh_bondk'])
	default_param(dic,'kT',        dic['hp_eps'])
	default_param(dic,'mu',        -8.*dic['hp_eps'])


	default_param(dic,'gen_bondk', 0.)
	default_param(dic,'shcp_bondk', 0.)
	default_param(dic,'theta0',    2.*np.arcsin(dic['l0']/np.sqrt(12*dic['R']**2 - 3*dic['l0']**2)))
	default_param(dic,'idealT',    16.*np.pi*(dic['Rg'] if dic['LJ'] else dic['R'])**2/np.sqrt(3))
	default_param(dic,'thetagrow', dic['theta0']/np.pi*180.)
	dic['thetagrow']=dic['thetagrow']/180.*np.pi
	dic['theta_lower']=dic['theta_lower']/180.*np.pi
	default_param(dic,'sigmaE',    dic['kT'])
	# dic['theta_upper']=dic['theta_upper']/180.*np.pi
	# dic['theta_tri']=dic['theta_tri']/180.*np.pi


	#Dynamic section:
	default_param(dic,'L',          20 * dic['R'])
	# default_param(dic,'lj_sigma',   max(dic['Rg']+dic['Rs'],np.sqrt((dic['l0']**2)/3. + dic['Rg']**2)))
	# default_param(dic,'lj_sigma',   dic['Rg']+dic['Rs'])
	default_param(dic,'lj_alpha',   2.)
	default_param(dic,'lj_rcut',    5.)
	
	default_param(dic,'yukawa_eps', dic['bjerrum']*dic['debye']*np.exp(2*dic['Rg']/dic['debye'])/(2*dic['Rg']+dic['debye']))
	default_param(dic,'yukawa_kappa',1./dic['debye'])
	default_param(dic,'yukawa_rcut', 3.* dic['debye'])
	default_param(dic,'yukawa_ron',  dic['debye'])
	
	# Monte Carlo
	# default_param(dic,'ZZ', np.exp(dic['mu']/dic['kT'])/dic['lambda'])
	
	return dic

prompt = '> '
path=os.path.abspath(os.getcwd())
#read parent parameter file
with open(path+'/parent_param.yaml', 'r') as f:
	    doc = yaml.load(f)
value=[]
subkey=[]
overwriteflag=False
for key in doc.keys():
	if (key!='Name') & (key!='Path'):
		subkey+=doc[key].keys()
		value+=doc[key].values()
#give birth to child parameter:
for i,sub_value in enumerate(itertools.product(*value)):
	#generate child dictionary:
	child_dict=dict(zip(subkey,sub_value))
	child_name=give_child_name(child_dict,doc['Name']['key'])
	print(child_name)
	child_path=os.path.join(path,child_name)
	print(path,child_path)
	if not os.path.exists(child_path):
		os.makedirs(child_path)
	elif overwriteflag:
		shutil.rmtree(child_path)
		os.makedirs(child_path)
	else: 
		print('Child Directory Exist\nDo you want to cover?')
		cover=raw_input(prompt)
		if cover=='yes':shutil.rmtree(child_path);os.makedirs(child_path)
		elif cover=='allyes': shutil.rmtree(child_path);os.makedirs(child_path);overwriteflag=True
		elif cover=='no': print('overwrite dir reject')
		elif cover=='allno':print('overwrite reject for all'); break
		else: print('unknown command')
	child_dict=train_child_param(child_dict)
	generate_child_param(child_name,child_dict,child_path)
