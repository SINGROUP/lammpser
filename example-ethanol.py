## Make a box with some ethanol molecules.
# All interaction parameters are taken from CHARMM.
#
import lammpser
from lammpser import Vector3, Quaternion
import math, random

# create an empty system
s = lammpser.System('booze')

# set the size of the simulation box
s.box = Vector3(21.3234,21.3234,21.3234)

# declare atomic types
s.AddAtomType(name='CG321',mass=12, charge=0.05,	epsilon=0.0560, sigma=2*2.0100 /1.122)
s.AddAtomType(name='CG331',mass=12, charge=-0.27,	epsilon=0.0780, sigma=2*2.0500 /1.122)
s.AddAtomType(name='OG311',mass=16, charge=-0.65,	epsilon=0.1921, sigma=2*1.7650 /1.122)
s.AddAtomType(name='HGP1',mass=1.008, charge=0.42,	epsilon=0.0460, sigma=2*0.2245 /1.122)
s.AddAtomType(name='HGA2',mass=1.008, charge=0.09,	epsilon=0.0350, sigma=2*1.3400 /1.122)
s.AddAtomType(name='HGA3',mass=1.008, charge=0.09,	epsilon=0.0240, sigma=2*1.3400 /1.122)


# declare bond types - C-Cr is just a name
s.bondTypes['C-C'] = [222.50, 1.5280] # CG321  CG331   222.50     1.5280 ! PROT alkane update, adm jr., 3/2/92
s.bondTypes['C-O'] = [428.00, 1.4200] # CG321  OG311   428.00     1.4200 ! PROT methanol vib fit EMB 11/21/89
s.bondTypes['O-H'] = [545.00, 0.9600] # OG311  HGP1    545.00     0.9600 ! PROT EMB 11/21/89 methanol vib fit; og tested on MeOH EtOH,...
s.bondTypes['CH2-H'] = [309.00, 1.1110] # CG321  HGA2    309.00     1.1110 ! PROT alkane update, adm jr., 3/2/92
s.bondTypes['CH3-H'] = [322.00, 1.1110] # CG331  HGA3    322.00     1.1110 ! PROT alkane update, adm jr., 3/2/92

s.angleTypes['CH2-OH-H'] = [50.00, 106, 0,0] 	# CG321  OG311  HGP1     50.00    106.00 
s.angleTypes['OH-CH2-H'] = [45.90, 108.89, 0,0]	# OG311  CG321  HGA2     45.90    108.89 
s.angleTypes['H-CH2-H'] = [35.50, 109.00, 5.40, 1.802] # HGA2   CG321  HGA2     35.50    109.00     5.40   1.802
s.angleTypes['OH-CH2-CH3'] = [75.70, 110.10, 0,0] # CG331  CG321  OG311    75.70    110.10
s.angleTypes['CH3-CH2-H'] = [34.60, 110.10, 22.53, 2.17900] # CG331  CG321  HGA2     34.60    110.10   22.53   2.17900
s.angleTypes['CH2-CH3-H'] = [34.60, 110.10, 22.53, 2.17900] # CG321  CG331  HGA3     34.60    110.10   22.53   2.17900 ! PROT alkane update, adm jr., 3/2/92
s.angleTypes['H-CH3-H'] = [35.50, 108.40, 5.40, 1.802] # HGA3   CG331  HGA3     35.50    108.40    5.40   1.80200 ! PROT alkane update, adm jr., 3/2/92

s.dihedralTypes['type1'] = [0.18, 3, 0,0.0]
s.dihedralTypes['type2'] = [0.16, 3, 0,0.0]
s.dihedralTypes['type3'] = [0.16, 3, 0,0.0]
s.dihedralTypes['type41'] = [1.1300,1,0,0.0] # CG331  CG321  OG311  HGP1       1.1300  1     0.00 ! og ethanol
s.dihedralTypes['type42'] = [0.1400,2,0,0.0] # CG331  CG321  OG311  HGP1       0.1400  2     0.00 ! og ethanol
s.dihedralTypes['type43'] = [0.2400,3,0,0.0] # CG331  CG321  OG311  HGP1       0.2400  3     0.00 ! og ethanol


# declare a molecular type
ethanol = s.AddMoleculeType('ethanol')

# put atoms in the molecule - give unique ID, atom type ID, and position in the molecule
ethanol.AddAtom('C1', 'CG321',  Vector3(10.15,11,9.6))
ethanol.AddAtom('O1', 'OG311',  Vector3(9,10.6,10.4))
ethanol.AddAtom('HO1', 'HGP1',  Vector3(9,11.27,11.1))
ethanol.AddAtom('H11', 'HGA2',  Vector3(10.95,11.5,10.2))
ethanol.AddAtom('H12', 'HGA2',  Vector3(9.87,12,9))

ethanol.AddAtom('C2', 'CG331',  Vector3(10.69,10,8.7))
ethanol.AddAtom('H21', 'HGA3',  Vector3(11.45,10.5,8))
ethanol.AddAtom('H22', 'HGA3',  Vector3(11.17,9.13,9.21))
ethanol.AddAtom('H23', 'HGA3',  Vector3(9.9,9.5,8.18))

#dihedrals
# type 1 dihed   HGA2   CG321  OG311  HGP1       0.1800  3     0.00 ! og methanol
ethanol.dihedrals.append( ['type1', 'HO1','O1','C1','H11'] )
ethanol.dihedrals.append( ['type1', 'HO1','O1','C1','H12'] )

# type 2 dihed   OG311  CG321  CG331  HGA3       0.1600  3     0.00 ! PROT rotation barrier in Ethane (SF)
ethanol.dihedrals.append( ['type2', 'O1','C1','C2','H21'] )
ethanol.dihedrals.append( ['type2', 'O1','C1','C2','H22'] )
ethanol.dihedrals.append( ['type2', 'O1','C1','C2','H23'] )

#type 3 dihed    HGA2   CG321  CG331  HGA3       0.1600  3     0.00 ! PROT rotation barrier in Ethane (SF)
ethanol.dihedrals.append( ['type3', 'H11','C1','C2','H21'] )
ethanol.dihedrals.append( ['type3', 'H11','C1','C2','H22'] )
ethanol.dihedrals.append( ['type3', 'H11','C1','C2','H23'] )
ethanol.dihedrals.append( ['type3', 'H12','C1','C2','H21'] )
ethanol.dihedrals.append( ['type3', 'H12','C1','C2','H22'] )
ethanol.dihedrals.append( ['type3', 'H12','C1','C2','H23'] )

ethanol.dihedrals.append( ['type41', 'C2','C1','O1','HO1'] )
ethanol.dihedrals.append( ['type42', 'C2','C1','O1','HO1'] )
ethanol.dihedrals.append( ['type43', 'C2','C1','O1','HO1'] )


# angles
ethanol.angles.append( ['CH2-OH-H', 'C1', 'O1', 'HO1'] )
ethanol.angles.append( ['OH-CH2-H', 'O1', 'C1', 'H11'] )
ethanol.angles.append( ['OH-CH2-H', 'O1', 'C1', 'H12'] )
ethanol.angles.append( ['H-CH2-H', 'H11', 'C1', 'H12'] )
ethanol.angles.append( ['OH-CH2-CH3', 'O1', 'C1', 'C2'] )
ethanol.angles.append( ['CH3-CH2-H', 'C2', 'C1', 'H11'] )
ethanol.angles.append( ['CH3-CH2-H', 'C2', 'C1', 'H12'] )

ethanol.angles.append( ['CH2-CH3-H', 'C1', 'C2', 'H21'] )
ethanol.angles.append( ['CH2-CH3-H', 'C1', 'C2', 'H22'] )
ethanol.angles.append( ['CH2-CH3-H', 'C1', 'C2', 'H23'] )

ethanol.angles.append( ['H-CH3-H', 'H21', 'C2', 'H22'] )
ethanol.angles.append( ['H-CH3-H', 'H21', 'C2', 'H23'] )
ethanol.angles.append( ['H-CH3-H', 'H22', 'C2', 'H23'] )

# create a bond - specify bond type name, atom1 ID and atom2 ID
ethanol.bonds.append( ['C-C','C1','C2'] )
ethanol.bonds.append( ['C-O','C1','O1'] )
ethanol.bonds.append( ['O-H','O1','HO1'] )
ethanol.bonds.append( ['CH2-H','C1','H11'] )
ethanol.bonds.append( ['CH2-H','C1','H12'] )
ethanol.bonds.append( ['CH3-H','C2','H21'] )
ethanol.bonds.append( ['CH3-H','C2','H22'] )
ethanol.bonds.append( ['CH3-H','C2','H23'] )

# place 27 molecules in a cubic arrangement into the box
c = 1
nc = 3
for x in xrange(nc):
        for y in xrange(nc):
		for z in xrange(nc):
			m = s.MakeMolecule('molecule_'+str(c), 'ethanol')
			m.position = Vector3(x*s.box.x/nc + 10, y*s.box.y/nc+ 10, z*s.box.z/nc + 10)
			m.rotation = Quaternion.Random()
			c += 1 #counter



# write the lammps topology
s.ExportLAMMPS('ethanol1.topo')
s.ExportXYZ('ethanol1.xyz')


