###################################################################################
## This script was written by Austin Meyer in the lab of Dr. Claus Wilke at the  ##
## University of Texas at Austin.  For any questions or concerns you can email   ##
## Austin at austin.g.meyer@gmail.com.                                           ##
###################################################################################

######################################################################
## Some of these includes are probably unecessary.                  ##
######################################################################
import sys, subprocess

from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

######################################################################
## This function accepts the full unfettered sequence file full pdb ##
## file, and calls the above translate_align function to do its     ##
## job.  Then, it calls the rsa determining function above and the  ##
## header function.  It combines the information to output a single ##
## file with the codon, residue, and rsa.  There are a bunch of     ##
## other columns just to satisfy the format of the script Mike has  ##
## already written.                                                 ##
######################################################################

def main():

  args =  sys.argv
  pdb =  args[1]

  structure = parsePDBStructure( pdb )
  (new_pdb, chain_renumbered_pdb) = renumberResidues( structure )
    
  return 0

######################################################################
## This function returns the structure attribute of a PDB file      ##
## from the PDB parser module.                                      ##
######################################################################

def parsePDBStructure( pdb_id ):
    parser = PDBParser()
    structure = parser.get_structure('test_rsa', pdb_id)
    return structure

######################################################################
## This function allows me to renumber the residues in a chain      ##
## to fix a particular odd numbering problem in Neuraminidase       ##
######################################################################

def renumberResidues( structure ):
  model = structure[0]

  i = 1
  for chain in model:
    for residue in chain:
      residue.id = (' ', i, ' ')
      i+=1
    
  new_pdb = 'Renumbered_Structure.pdb'
  w = PDBIO()
  w.set_structure(structure)
  w.save(new_pdb)

  return (new_pdb, structure)

## Execute the main function
if __name__ == "__main__":
    main()

