# @File: ovitoanalyze.py
# @Author: Siyu Li (email:lisiyu211@gmail.com)
# @Date:   2018-03-19 18:20:01
# @Last Modified time: 2018-07-24 18:11:04
from ovito.io import *
from ovito.data import *
from ovito.modifiers import *
from ovito.vis import *
import numpy as np
import sys,os,glob
from io import StringIO
def shellE():
	with open('out','r') as f:
		for rownum,line in enumerate(f,start=1):
			if 'Final ShellE' in line:
				print(line,'blue')
				finalE=line.split('[')[-1]
				finalE=np.loadtxt(StringIO(finalE))
				Elasticity=sum(finalE[3:5])
			if 'Final Subunit' in line:
				Svalue=int(line.split(':')[1])
				print(line,'blue')
			if 'Elasticity' in line:
				if 'Final Elasticity' in line:
					Evalue=float(line.split(':')[1])
				print(line,'blue')
	print('E/units:%.5f/%d=%.5f'%(Elasticity,Svalue,Elasticity/Svalue))
	print('E/units:%.5f/%d=%.5f'%(Evalue,Svalue,Evalue/Svalue))
	return Elasticity,Svalue

def child_analyze():
	#compute shape info
	node = import_file('final.gsd')
	#find vertex neighbors
	parN=node.source.number_of_particles
	bond_property=node.source.bonds.array
	table=np.zeros((parN,parN),dtype='b')
	for i,j in bond_property:
		table[i,j]=1
	vneigh=np.sum(table,axis=0)
	tnum=int(round(len(bond_property)/6))

	#vertex connectivity determinant
	vertex=np.nonzero(vneigh!=0)[0]
	vnum=len(vertex)
	matrix=table[vertex[:,None],vertex]
	det=int(round(np.linalg.det(matrix)))
	
	#find pentamer & hexamer number, & triangle number
	vpen=np.nonzero(vneigh==5)[0]
	vhex=np.nonzero(vneigh==6)[0]
	vgen=0
	vtrash=np.nonzero(vneigh==0)[0][1:]
	npen=len(vpen)
	nhex=len(vhex)
	print("pen_v:{}, hex_v:{}, tot_v:{}, subunits:{}, det:{}".format(npen,nhex,vnum,tnum,det))
	
	pos=node.source.particle_properties['Position'].marray
	center=sum(pos[vertex])/vnum
	radius=np.mean(np.linalg.norm(pos[vertex]-center,axis=1))
	print('radius:',radius)

	#define shell type
	charv=np.int_([12,14,15,16,20,20,22,23,27,32,42,42])
	chart=np.int_([20,24,26,28,36,36,40,42,50,60,80,80])
	chardet=np.array([625,0,0,0,-14400,-24167,-152100,235224,84700,0,2621440000,3340840000])
	char=np.array(['A','B','C','D','E','E*','E**','3fold','F','G','H','H*'])
	findtype=np.nonzero(((charv==vnum)&(chart==tnum)&(chardet==det)))[0]
	if len(findtype)==0: shelltype='IRG'
	else:shelltype=char[findtype][0]
	print('shelltype:',shelltype)
	with open('shellinfo.txt', 'w') as f:
		f.write('shelltype:{}\nradius:{}\npen_v:{}\nhex_v:{}\ntot_v:{}\nsubunits:{}\ndet:{}\n'.format(shelltype,radius,npen,nhex,vnum,tnum,det))
		#read shell energy
		E,N=shellE()
		f.write('Elasticity:{},Subunits:{}\nE/S:{}\n'.format(E,N,E/N))
		f.write('Matrix:\n{}'.format(matrix))
		



if len(sys.argv)==1: childdir=os.listdir() 
else: childdir=sys.argv[1:]
for subdir in childdir:
	if os.path.isdir(subdir):
		os.chdir(subdir)
		for resultfile in glob.glob('final1.gsd'):
			print(subdir)
			child_analyze()
		os.chdir('..')

# cell = node.source.cell
# cell.display.enabled = False
# node.add_to_scene()
# vp = Viewport()
# # vp.type = Viewport.Type.PERSPECTIVE
# # vp.camera_pos = (-10, -15, 15)
# # vp.camera_dir = (2, 3, -3)
# # vp.fov = math.radians(60.0)
# settings = RenderSettings(
#     filename = "myimage.png",
#     size = (2000, 1500)
# )
# vp.render(settings)

# type_property = node.source.particle_properties.particle_type
# type_property.marray[vpen]=3
# print(type_property.array)
# print(node.source.particle_properties.particle_type.array)
# node.source.particle_properties.particle_type.changed()
# for particle in vhex:
#     # ptype_id = type_property.array[atom_index]
# for ptype in type_property.type_list:
#         if ptype.id == ptype_id:
        # print(ptype.color)
#             break