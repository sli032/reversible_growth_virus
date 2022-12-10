import hoomd
from hoomd import md
import virus
import numpy as np

# class snap(object):
# 	"""docstring for snap"""
# 	def __init__(self, arg):
# 		super(snap, self).__init__()
# 		self.arg = arg
		

def initial_snap(Boxsize,genome,shell,constituent=None):
	hoomd.context.initialize()
	snap=hoomd.data.make_snapshot(N=0,
		box=hoomd.data.boxdim(L=Boxsize),
		particle_types=['G','S','C'], #P:polymer; X:shell; C:constituent
		bond_types=['genome','shell','shcp'],
		dihedral_types=['capsomer'],
		angle_types=['ds'],
		dtype='double')
	genstart,genend=assign_obj_to_snap(snap,genome,cover=False,particle=True,bond=False,dihedral=False)
	shstart,shend=assign_obj_to_snap(snap,shell,cover=False,dihedral=True)
	print('genstart,genend:',genstart,genend)
	print('shstart,shend:',shstart,shend)
	if constituent!=None:
		cpstart,cpend=assign_obj_to_snap(snap,constituent,cover=False,particle=True,bond=True,dihedral=False)
		print('cpstart,cpend:',cpstart,cpend)
	return snap

def assign_snap_particles(snap,obj,cover):
	# particles assign in range:[pstart,pend)
	if cover: 
		for i,pid in enumerate(snap.particles.typeid):
			if pid==obj.typeid: 
				pstart=i
				break
	else: pstart=snap.particles.N
	
	pend=pstart+obj.n
	snap.particles.resize(pend)

	for i, [cor] in enumerate(zip(obj.particles),start=pstart):
		# print('cor:{},pid:{}'.format(cor,pid))
		snap.particles.diameter[i] = 2*obj.rp
		snap.particles.position[i]=cor
		snap.particles.typeid[i] = obj.typeid
		snap.particles.body[i] = obj.bodyid
	return pstart,pend
def assign_snap_bonds(snap,obj,pstart,cover):
	if cover: bstart=int(np.nonzero(snap.bonds.typeid==obj.bondid)[0][0])
	else: bstart=snap.bonds.N

	snap.bonds.resize(bstart+len(obj.bonds))

	for i, bond in enumerate(obj.bonds,start=bstart):
		snap.bonds.group[i]=pstart+bond
		snap.bonds.typeid[i]=obj.bondid
def assign_snap_dihedrals(snap,obj,pstart,cover):
	if snap.dihedrals.N>0:
		if cover: 
			dstart=int(np.nonzero(snap.dihedrals.typeid==obj.dihedralid)[0][0])
		else: dstart=snap.dihedrals.N
	else: dstart=0

	snap.dihedrals.resize(dstart+len(obj.dihedrals))
	for i, dihedral in enumerate(obj.dihedrals,start=dstart):
		snap.dihedrals.group[i]=pstart+dihedral
		snap.dihedrals.typeid[i]=obj.dihedralid

def assign_obj_to_snap(snap,obj,cover=True,particle=True,bond=True,dihedral=False):
	#!!assign bond and dihedral before particle
	#particle number needed for re-index bond and dihedral
	if particle: pstart,pend=assign_snap_particles(snap,obj,cover)
	if bond: assign_snap_bonds(snap,obj,pstart,cover)
	if dihedral: assign_snap_dihedrals(snap,obj,pstart,cover)
	return pstart,pend
def snapprint(snap,particle=True,bond=True,dihedral=True):
	if particle: 
		print('particle position:\n{} \
			   \nparticle typeid:\n{} \
			   \npartical body:\n{} \
			   \nparticale diameter:\n{}' \
			   .format(snap.particles.position[:],snap.particles.typeid[:],snap.particles.body[:],snap.particles.diameter[:]))
	if bond: 
		print('bonds group:\n{} \
			   \nbonds typeid:\n{}' \
			   .format(snap.bonds.group[:],snap.bonds.typeid[:]))
	if dihedral: 
		print('dihedrals group:\n{} \
			   \ndihedrals typeid:\n{}' \
			   .format(snap.dihedrals.group[:],snap.dihedrals.typeid[:]))