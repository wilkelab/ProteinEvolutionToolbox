#!/usr/bin/python

from Bio.PDB import *

def getResNumber(residue):
	id = residue.get_id()
	assert(id[0]==' ') # just make sure we're working with a properly numbered PDB file
	assert(id[2]==' ')
	return id[1]


pdbfile='Renumbered_Structure.pdb'

parser=PDBParser()

structure=parser.get_structure('', pdbfile)
model=structure[0]

for chain in model:
    for r1 in chain:
        i = 0
        for r2 in chain:
            if getResNumber(r2) != getResNumber(r1):
                i += 1/(float(r1['CA']-r2['CA'])**2)
        print(i)
