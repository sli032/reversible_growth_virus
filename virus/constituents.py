import numpy as np

class constituents():
	"""docstring for Shell"""
	def __init__(self,r,cors,bonds):
		# super(polymer, self).__init__()
		self.rp=r
		self.particles=cors
		self.n=len(self.particles)
		self.typeid=2
		# self.pids=self.typeid*self.n #typeid
		self.bodyid=-1 #bodyid
		self.bonds=bonds
		self.bondid=2
	def update(self,shell):
		cors,bonds=[],[]
		i=0
		for ti,row in shell.triangle.iterrows():
			[v0,v1,v2]=row[:3]
			[cor0,cor1,cor2]=np.array(shell.vertex.loc[[v0,v1,v2],list('xyz')])
			cors=cors+[(cor0+cor1+cor2)/3.]
			bonds=bonds+list(map(lambda x: [i,shell.pdic[x]-shell.n], [v0,v1,v2]))
			i+=1
		self.particles=cors
		self.n=len(self.particles)
		# self.pids=self.typeid*self.n #typeid
		self.bonds=np.int_(bonds)
def create_constituents(shell,Rc):
	cors,bonds=[],[]
	i=0
	print(shell.n)
	for ti,row in shell.triangle.iterrows():
		[v0,v1,v2]=row[:3]
		[cor0,cor1,cor2]=np.array(shell.vertex.loc[[v0,v1,v2],list('xyz')])
		cors=cors+[(cor0+cor1+cor2)/3.]
		bonds=bonds+list(map(lambda x: [i,shell.pdic[x]-shell.n], [v0,v1,v2]))
		i+=1
	return constituents(Rc,cors,np.int_(bonds))