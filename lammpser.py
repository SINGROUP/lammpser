import os
import math
import random
from collections import OrderedDict


class Vector3(object):

	def __init__(self,x,y,z):
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)

	def Clone(self):
		return Vector3(self.x, self.y, self.z)

	def Normalise(self):
		norm = self.x*self.x + self.y*self.y + self.z*self.z
		norm = math.sqrt(norm)
		self.x /= norm
		self.y /= norm
		self.z /= norm

	def Norm(self):
		return math.sqrt(self.SqrNorm)
	def SqrNorm(self):
		return self.x*self.x + self.y*self.y + self.z*self.z

	def __add__(self, other):
		return Vector3(self.x+other.x, self.y+other.y, self.z+other.z)
	def __sub__(self, other):
		return Vector3(self.x-other.x, self.y-other.y, self.z-other.z)


	def __str__(self):
		return str(self.x)+' '+str(self.y)+' '+str(self.z)

	@staticmethod
	def Rotation(v,q):

		p = Vector3(0,0,0)
		p.x = (q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z)*v.x
		p.x += 2*(q.x*q.y - q.w*q.z)*v.y
		p.x += 2*(q.x*q.z + q.w*q.y)*v.z

		p.y = 2*(q.x*q.y + q.w*q.z)*v.x
		p.y += (q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z)*v.y
		p.y += 2*(q.y*q.z - q.w*q.x)*v.z

		p.z = 2*(q.x*q.z - q.w*q.y)*v.x
		p.z += 2*(q.y*q.y + q.w*q.x)*v.y
		p.z += (q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z)*v.z

		return p

	def Rotate(self,q):
		p = Vector3.Rotation(self,q)
		self.x=p.x
		self.y=p.y
		self.z=p.z

	@staticmethod
	def PBC(v,box):
		p = v.Clone()
		p.x -= math.floor(p.x/box.x)*box.x
		p.y -= math.floor(p.y/box.y)*box.y
		p.z -= math.floor(p.z/box.z)*box.z
		return p

	@staticmethod
	def PBCDistance(u,v,box):

		p = u-v;
		p.x -= math.round(p.x/box.x)*box.x
		p.y -= math.round(p.y/box.y)*box.y
		p.z -= math.round(p.z/box.z)*box.z
		return p.Norm


	@staticmethod
	def Dot(u,v):
		return u.x*v.x + u.y*v.y + u.z*v.z



	@staticmethod
	def RandomDirection():

		theta = random.uniform(-math.pi, math.pi)*0.5
		phi = random.uniform(0,math.pi)*2
		v = Vector3(math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta))
		v.Normalise()
		return v

	@staticmethod
	def RandomInBox(box):
		return Vector3(random.uniform(0,box.x), random.uniform(0,box.y), random.uniform(0,box.z))


class Quaternion(object):

	@staticmethod
	def Identity():
		return Quaternion(1,0,0,0)

	@staticmethod
	def FromAngleAxis(angle, axis):
		q = Quaternion.Identity()
		axis.Normalise()
		q.w = math.cos(0.5*angle)
		q.x = axis.x * math.sin(0.5*angle)
		q.y = axis.y * math.sin(0.5*angle)
		q.z = axis.z * math.sin(0.5*angle)
		return q


	def __init__(self,w,x,y,z):
		self.w = w
		self.x = x
		self.y = y
		self.z = z

	@staticmethod
	def Random():
		q = Quaternion.FromAngleAxis(random.uniform(0,2*math.pi), Vector3.RandomDirection())
		return q





class AtomType:

	def __init__(self, name, mass, charge, epsilon=0, sigma=1):
		self.name = name
		self.mass = mass
		self.charge = charge
		self.epsilon = epsilon
		self.sigma = sigma
		
class MoleculeType:

	def __init__(self, ID):

		self.ID = ID
		self.atoms = OrderedDict()

		self.bonds = []
		self.angles = []
		self.dihedrals = []
		self.impropers = []


	def AddAtom(self, ID, typeID, position):

		if ID in self.atoms.keys():
			raise NameError("ERROR: atom["+ID+"] is already in the molecule.")

		##TODO: check that the typeID exists? naaaa
		atom = Atom(ID,typeID,position)
		self.atoms[ID] = atom
		return atom


	def Print(self):

		print 'MoleculeType ['+self.ID+']:'
		print 'Atoms:'
		for aID in self.atoms.keys():
			a = self.atoms[aID]
			a.Print()
		print 'Bonds:'
		print str(self.bonds)




class Atom:
	pass

class Molecule:
	pass

class System:
	pass


## class for an atom in a molecule template
class Atom:


	def __init__(self, ID, typeID, position=Vector3(0,0,0)):

		self.ID = ID
		self.typeID = typeID
		self.position = position
		self.index = -1

	def Print(self):
		print 'Atom ['+self.ID+'] type '+self.typeID+'\t-> '+str(self.position)
		
	def Clone(self):
		atom = Atom(self.ID, self.typeID, self.position)
		return atom

## class for a molecule in the topology
class Molecule:

	def __init__(self, ID, typeID, position=Vector3(0,0,0), rotation=Quaternion.Identity()):

		self.index = -1
		self.ID = ID
		self.typeID = typeID
		self.position = position
		self.rotation = rotation
		self.atoms = OrderedDict()


	def GetAtomPosition(self, atomID, system):

		atom = system.moleculeTypes[self.typeID].atoms[atomID]
		p = atom.position.Clone()

		p = p.Rotation(self.rotation)

		p = p+self.position

		return p

	def Initialise(self, molType, startID):

		print 'mol init',molType.ID,molType.atoms.keys()
		self.atoms = OrderedDict()
		idcnt = startID
		for aID in molType.atoms.keys(): #take all atoms from the moltype and replicate them
			a = molType.atoms[aID]

			atm = a.Clone(); #clone and set the position
			atm.position = Vector3.Rotation(a.position,self.rotation) + self.position
			atm.index = idcnt

			self.atoms[aID] = atm

			idcnt+=1
			


class System(object):

	def __init__(self, ID):

		self.ID = ID

		self.__atomtypes = OrderedDict()
		self.topology = OrderedDict()
		self.moleculeTypes = OrderedDict()

		self.box = Vector3(1,1,1)

		self.bondTypes = OrderedDict()
		self.angleTypes = OrderedDict()
		self.dihedralTypes = OrderedDict()
		self.improperTypes = OrderedDict()


		print 'System '+ID+' created.'


	def AddAtomType(self, **kwargs):

		if not ('name' in kwargs.keys()):
			raise NameError("ERROR: atomtype must have a name.")

		name = kwargs['name']

		if name in self.__atomtypes.keys():
			raise NameError("ERROR: "+name+" is already in the atom types list.")

		if not ('mass' in kwargs.keys()):
			print "WARNING: atom type mass is not given: using 1"
			mass = 1
		else:
			mass = float(kwargs['mass'])

		if not ('charge' in kwargs.keys()):
			print "WARNING: atom type charge is not given: using 0"
			charge = 0
		else:
			charge = float(kwargs['charge'])

		if not ('epsilon' in kwargs.keys()):
			print "WARNING: atom type LJ epsilon is not given: using 0"
			epsilon = 0
		else:
			epsilon = float(kwargs['epsilon'])

		if not ('sigma' in kwargs.keys()):
			print "WARNING: atom type LJ sigma is not given: using 1"
			sigma = 0
		else:
			sigma = float(kwargs['sigma'])

		self.__atomtypes[name] = AtomType(name,mass,charge,epsilon,sigma)
		print 'Atom type '+name+' created.'


	def RemoveAtomType(self, name):

		if not(name in self.__atomtypes.keys()):
			print "WARNING: "+name+" is not in the atom types list."
			return

		del self.__atomtypes[name] 
		print 'Atom type '+name+' removed.'


	def AddMoleculeType(self,name):

		if name in self.moleculeTypes.keys():
			raise NameError("ERROR: molcule type["+name+"] was already defined.")

		mol = MoleculeType(name)
		self.moleculeTypes[name] = mol

		print 'Molecule template '+name+' created.'
		return mol


	def MakeMolecule(self, ID, typeID, position=Vector3(0,0,0), rotation=Quaternion.Identity()):

		if not(typeID in self.moleculeTypes.keys()):
			raise NameError('ERROR: molecule type '+typeID+' not defined.')

		if ID in self.topology.keys():
			raise NameError('ERROR: molecule '+ID+' is already in the topology.')

		mol = Molecule(ID,typeID,position,rotation)
		self.topology[ID] = mol

		print 'molecule created and added to topology.'
		return mol

	def ExportLAMMPS(self, filename):

		f = open(filename,'w')

		f.write('system '+self.ID+' created by LAMMPSer\n\n')

		#compute the number of atoms
		counter = 0
		for molID in self.topology.keys():
			mol = self.topology[molID]
			counter += len(self.moleculeTypes[mol.typeID].atoms.keys())
		print 'total # atoms: ',counter
		f.write(str(counter)+' atoms\n')

		#compute the number of bonds
		counter = 0
		for molID in self.topology.keys():
			mol = self.topology[molID]
			counter += len(self.moleculeTypes[mol.typeID].bonds)
		print 'total # bonds: ',counter
		f.write(str(counter)+' bonds\n')

		#compute the number of angles
		counter = 0
		for molID in self.topology.keys():
			mol = self.topology[molID]
			counter += len(self.moleculeTypes[mol.typeID].angles)
		print 'total # angles: ',counter
		f.write(str(counter)+' angles\n')

		#compute the number of diheds
		counter = 0
		for molID in self.topology.keys():
			mol = self.topology[molID]
			counter += len(self.moleculeTypes[mol.typeID].dihedrals)
		print 'total # dihedrals: ',counter
		f.write(str(counter)+' dihedrals\n')

		#compute the number of improps
		counter = 0
		for molID in self.topology.keys():
			mol = self.topology[molID]
			counter += len(self.moleculeTypes[mol.typeID].impropers)
		print 'total # impropers: ',counter
		f.write(str(counter)+' impropers\n')


		#compute the number of  Types
		f.write(str(len(self.__atomtypes))+' atom types\n')
		f.write(str(len(self.bondTypes))+' bond types\n')
		f.write(str(len(self.angleTypes))+' angle types\n')
		f.write(str(len(self.dihedralTypes))+' dihedral types\n')
		f.write(str(len(self.improperTypes))+' improper types\n\n')

		f.write('0 '+str(self.box.x)+' xlo xhi\n')
		f.write('0 '+str(self.box.y)+' ylo yhi\n')
		f.write('0 '+str(self.box.z)+' zlo zhi\n\nMasses\n\n')

		#determine the type list
		atomtypelist = self.__atomtypes.keys()
		counter = 1
		for atID in atomtypelist:
			f.write(str(counter)+' '+str(self.__atomtypes[atID].mass)+'\n')
			counter += 1
		f.write('\n')
		f.write('Pair Coeffs\n\n')
		counter = 1
		for atID in atomtypelist:
			f.write(str(counter)+'\t'+str(self.__atomtypes[atID].epsilon)+'\t'+str(self.__atomtypes[atID].sigma)+'\n')
			counter += 1
		f.write('\n')

		#determine the type list - bonds
		bondtypelist = self.bondTypes.keys()
		if len(bondtypelist) > 0:
			counter = 1
			f.write('Bond Coeffs\n\n')
			for bt in bondtypelist:
				b = self.bondTypes[bt]
				f.write(str(counter)+'\t'+str(b[0])+'\t'+str(b[1])+'\n')
				counter += 1
			f.write('\n')

		angletypelist = self.angleTypes.keys()
		if len(angletypelist) > 0:
			counter = 1
			f.write('Angle Coeffs\n\n')
			for bt in angletypelist:
				b = self.angleTypes[bt]
				f.write(str(counter)+'\t')
				for p in b:
					f.write(str(p)+'\t')
				f.write('\n')
				#f.write(str(counter)+'\t'+str(b[0])+'\t'+str(b[1])+'\n')
				counter += 1
			f.write('\n')		

		dihedtypelist = self.dihedralTypes.keys()
		if len(dihedtypelist) > 0:
			counter = 1
			f.write('Dihedral Coeffs\n\n')
			for bt in dihedtypelist:
				b = self.dihedralTypes[bt]
				f.write(str(counter)+'\t')
				for p in b:
					f.write(str(p)+'\t')
				f.write('\n')
				counter += 1
			f.write('\n')		

		impstypelist = self.improperTypes.keys()
		if len(impstypelist) > 0:
			counter = 1
			f.write('Improper Coeffs\n\n')
			for bt in impstypelist:
				b = self.improperTypes[bt]
				f.write(str(counter)+'\t'+str(b[0])+'\t'+str(b[1])+'\t'+str(b[2])+'\n')
				counter += 1
			f.write('\n')		


		f.write('Atoms\n\n')
		counter = 1; molidx = 1
		for molID in self.topology.keys():
			mol = self.topology[molID]
			mol.index = molidx
			
			mol.Initialise(self.moleculeTypes[mol.typeID],counter)

			#for aID in self.moleculeTypes[mol.typeID].atoms.keys():
			for aID in mol.atoms.keys():
			#for atom in mol.atoms:
				atom = mol.atoms[aID]
				#atom = self.moleculeTypes[mol.typeID].atoms[aID]
				atomTemplate = self.__atomtypes[atom.typeID]
				f.write(str(atom.index)+'\t') #atom index
				f.write(str(mol.index)+'\t') #molecule index
				f.write(str(atomtypelist.index(atom.typeID)+1)+'\t') #atom type
				f.write(str(atomTemplate.charge)+'\t') #molecule index
				p = Vector3.PBC(atom.position, self.box)

				f.write(str(p)+'\n')
				counter += 1
			molidx += 1

		

		#write the bonds
		hasBonds = False
		for molID in self.topology.keys():
			mol = self.topology[molID]
			if len(self.moleculeTypes[mol.typeID].bonds) > 0:
				hasBonds = True
				break
		if hasBonds: #write the bonds if there are any
			f.write('\nBonds\n\n')
			counter = 1
			for molID in self.topology.keys():
				mol = self.topology[molID]
				molType = self.moleculeTypes[mol.typeID]

				for bnd in molType.bonds:
					
					bndIndex = bondtypelist.index(bnd[0])+1
					atom1ID = mol.atoms[bnd[1]].index
					atom2ID = mol.atoms[bnd[2]].index

					f.write(str(counter)+' '+str(bndIndex)+' ')
					f.write(str(atom1ID)+" "+str(atom2ID)+"\n")

					counter+=1

		#write the angles
		hasAngles = False
		for molID in self.topology.keys():
			mol = self.topology[molID]
			if len(self.moleculeTypes[mol.typeID].angles) > 0:
				hasAngles = True
				break
		if hasAngles: #write the bonds if there are any
			f.write('\nAngles\n\n')
			counter = 1
			for molID in self.topology.keys():
				mol = self.topology[molID]
				molType = self.moleculeTypes[mol.typeID]

				for bnd in molType.angles:
					
					bndIndex = angletypelist.index(bnd[0])+1
					atom1ID = mol.atoms[bnd[1]].index
					atom2ID = mol.atoms[bnd[2]].index
					atom3ID = mol.atoms[bnd[3]].index

					f.write(str(counter)+' '+str(bndIndex)+' ')
					f.write(str(atom1ID)+" "+str(atom2ID)+" "+str(atom3ID)+"\n")

					counter+=1

		#write the dihedrals
		hasDiheds = False
		for molID in self.topology.keys():
			mol = self.topology[molID]
			if len(self.moleculeTypes[mol.typeID].dihedrals) > 0:
				hasDiheds = True
				break
		if hasDiheds: #write the bonds if there are any
			f.write('\nDihedrals\n\n')
			counter = 1
			for molID in self.topology.keys():
				mol = self.topology[molID]
				molType = self.moleculeTypes[mol.typeID]

				for bnd in molType.dihedrals:
					
					bndIndex = dihedtypelist.index(bnd[0])+1
					atom1ID = mol.atoms[bnd[1]].index
					atom2ID = mol.atoms[bnd[2]].index
					atom3ID = mol.atoms[bnd[3]].index
					atom4ID = mol.atoms[bnd[4]].index

					f.write(str(counter)+' '+str(bndIndex)+' ')
					f.write(str(atom1ID)+" "+str(atom2ID)+" "+str(atom3ID)+" "+str(atom4ID)+"\n")

					counter+=1

		#write the dihedrals
		hasImps = False
		for molID in self.topology.keys():
			mol = self.topology[molID]
			if len(self.moleculeTypes[mol.typeID].impropers) > 0:
				hasImps = True
				break
		if hasImps: #write the bonds if there are any
			f.write('\nImpropers\n\n')
			counter = 1
			for molID in self.topology.keys():
				mol = self.topology[molID]
				molType = self.moleculeTypes[mol.typeID]

				for bnd in molType.dihedrals:
					
					bndIndex = dihedtypelist.index(bnd[0])+1
					atom1ID = mol.atoms[bnd[1]].index
					atom2ID = mol.atoms[bnd[2]].index
					atom3ID = mol.atoms[bnd[3]].index
					atom4ID = mol.atoms[bnd[4]].index

					f.write(str(counter)+' '+str(bndIndex)+' ')
					f.write(str(atom1ID)+" "+str(atom2ID)+" "+str(atom3ID)+" "+str(atom4ID)+"\n")

					counter+=1



		f.close()


	def ExportXYZ(self, filename):


		f = open(filename,'w')
		
		#compute the number of atoms
		counter = 0
		for molID in self.topology.keys():
			mol = self.topology[molID]
			counter += len(self.moleculeTypes[mol.typeID].atoms.keys())
		print 'total # atoms: ',counter
		f.write(str(counter)+'\n\n')

		counter = 1; molidx = 1
		for molID in self.topology.keys():
			mol = self.topology[molID]
			mol.index = molidx
			
			mol.Initialise(self.moleculeTypes[mol.typeID],counter)

			#for aID in self.moleculeTypes[mol.typeID].atoms.keys():
			for aID in mol.atoms.keys():
			#for atom in mol.atoms:
				atom = mol.atoms[aID]
				#atom = self.moleculeTypes[mol.typeID].atoms[aID]
				atomTemplate = self.__atomtypes[atom.typeID]
				f.write(str(atomTemplate.name)+'\t') #atom index
				p = Vector3.PBC(atom.position, self.box)
				f.write(str(p)+'\n')
				counter += 1
			molidx += 1




		f.close()












