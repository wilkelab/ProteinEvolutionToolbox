#!/usr/bin/python

import sys

from prody import *
from numpy import *

def loadFiles(pdb_name, dcd_name):
  dcd = Trajectory(dcd_name)
  structure = parsePDB(pdb_name)

  dcd.setCoords(structure)

  dcd.link(structure)

  return (dcd, structure)

def calculateRMSDs(dcd, structure, skip):

  rmsd_all_sites = zeros(len(structure.calpha))

  residue_number_in_pdb = (structure.calpha).getResnums()

  for site in range(len(residue_number_in_pdb)):
    print(residue_number_in_pdb[site])
    
    site_positions = structure.select('resnum ' + str(residue_number_in_pdb[site]))
    
    locations = zeros(3 * (len(dcd) - skip)).reshape(len(dcd) - skip, 3)
    
    dcd.reset()
    for i, frame in enumerate(dcd):
      if i < skip:
        continue
      else:
        locations[i - skip] = calcCenter(site_positions)

    center = mean(locations, axis = 0)
    rmsd_all_sites[site] = sqrt(mean(sum((locations - center)**2, axis=1)))
   
  return rmsd_all_sites

def calculateRMSDsCA(dcd, structure, skip):

  rmsd_all_sites = zeros(len(structure.calpha))

  residue_number_in_pdb = (structure.calpha).getResnums()

  for site in range(len(residue_number_in_pdb)):
    print(residue_number_in_pdb[site])
    
    site_positions = structure.select('ca resnum ' + str(residue_number_in_pdb[site]))
    
    locations = zeros(3 * (len(dcd) - skip)).reshape(len(dcd) - skip, 3)
    
    dcd.reset()
    for i, frame in enumerate(dcd):
      if i < skip:
        continue
      else:
        locations[i - skip] = site_positions.getCoords()

    center = mean(locations, axis = 0)
    rmsd_all_sites[site] = sqrt(mean(sum((locations - center)**2, axis=1)))
   
  return rmsd_all_sites
  
def generateRMSDFile(pdb_file, dcd_file, out_file, skip, ref_type):
  (dcd, structure) = loadFiles(pdb_file, dcd_file)

  #Declare variable
  all_rmsds = []
  
  if ref_type == 'COM':
    all_rmsds = calculateRMSDs(dcd, structure, skip)
  elif ref_type == 'CA':
    all_rmsds = calculateRMSDsCA(dcd, structure, skip)
  else:
    print "You did not pick either CA for C-alphas or COM for center of mass in the RMSF calculation"

  output_handle = open(out_file, 'w')

  output_handle.write('rmsd\n')

  for i in all_rmsds:
    output_handle.write(str(i) + '\n')

  output_handle.close()

  return 0

def main():
  if len( sys.argv ) != 6:
    print '''

    You screwed up choosing input files.

    '''
    print "     ", sys.argv[0], "<pdb file> <dcd file> <skip> <CA or COM> <output file>"
    
  else:
    pdb_file = sys.argv[1]
    dcd_file = sys.argv[2]
    skip     = int(sys.argv[3])
    ref_type = sys.argv[4]
    out_file = sys.argv[5]
    

    generateRMSDFile( pdb_file, dcd_file, out_file, skip, ref_type )

if __name__ == "__main__":
  main()

