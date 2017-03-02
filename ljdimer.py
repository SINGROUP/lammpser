import lammpser
from lammpser import Vector3, Quaternion
import math

# create an empty system
s = lammpser.System('asd')

# set the size of the simulation box
s.box = Vector3(25,25,25)

# declare atomic types
s.AddAtomType(name='A',mass=39.948, charge=0, epsilon=0.238, sigma=3.40)

# declare bond types - C-Cr is just a name
s.bondTypes['A-A'] = [200, 3.0] # 1.122*3.4 same parameters LAMMPS would require, in the same order

# declare a molecular type
dimer = s.AddMoleculeType('dimer')

# put atoms in the molecule - give unique ID, atom type ID, and position in the molecule
dimer.AddAtom('A1', 'A',  Vector3(0,0,0))
dimer.AddAtom('A2', 'A',  Vector3(0.3,0,0))

# create a bond - specify bond type name, atom1 ID and atom2 ID
dimer.bonds.append( ['A-A','A1','A2'] )

c = 1
for x in xrange(6):
	for y in xrange(6):
		for z in xrange(6):
			m = s.MakeMolecule('molecule_'+str(c), 'dimer')
			m.position = Vector3(x*6,y*6,z*6)
			m.rotation = Quaternion.Random()
			c += 1


# write the lammps topology
s.ExportLAMMPS('dimer_666.topo')


