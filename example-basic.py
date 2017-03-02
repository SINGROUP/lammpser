## Simple example that shows how to create a system and basic dimer molecules.

import lammpser
from lammpser import Vector3, Quaternion
import math

# create an empty system
s = lammpser.System('something')

# set the size of the simulation box
s.box = Vector3(2,2,2)

# declare atomic types
s.AddAtomType(name='C',mass=12, charge=0, epsilon=1, sigma=1)
s.AddAtomType(name='Cr',mass=12, charge=0, epsilon=1, sigma=1)

# declare bond types - C-Cr is just a name
s.bondTypes['C-Cr'] = [330, 1.2] #same parameter LAMMPS would require, in the same order

# declare angle/dihedral/improper types
# ...

# declare a molecular type
cmol = s.AddMoleculeType('myDimer')

# put atoms in the molecule - give unique ID, atom type ID, and position in the molecule
cmol.AddAtom('C1', 'C',  Vector3(0,0,0))
s.moleculeTypes['myDimer'].AddAtom('C2', 'Cr', Vector3(1,0,0))

# create a bond - specify bond type name, atom1 ID and atom2 ID
cmol.bonds.append( ['C-Cr','C1','C2'] )


# create a molecule in the system, using a declared type - give unique ID and a type ID
s.MakeMolecule('molecule1', 'myDimer')
# a copy of myDimer named molecule1 is now added to the topology of the system

# get a reference to the molecule by ID
mol1 = s.topology['molecule1']

# create another molecule
mol2 = s.MakeMolecule('molecule2', 'myDimer')

# set the position of the molecule
mol2.position = Vector3(3,0,0)
# set the orientation of the molecule, using angle-axis notation
mol2.rotation = Quaternion.FromAngleAxis(math.pi/2, Vector3(0,1,0))

# make a random rotation
mol2.rotation = Quaternion.Random()

# write the lammps topology
s.ExportLAMMPS('test.topo')


