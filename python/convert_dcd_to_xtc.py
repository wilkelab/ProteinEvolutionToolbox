#!/usr/bin/env python

from mdtraj import trajectory
import sys

def main():
  if len( sys.argv ) != 4:
    print '''

    You don't have the right number of arguments.

    '''
  else:
    pdb_name = sys.argv[1]
    dcd_name = sys.argv[2]
    xtc_name = sys.argv[3]

    convert_dcd_xtc(pdb_name, dcd_name, xtc_name)

def convert_dcd_xtc(pdb_name, dcd_name, xtc_name):
  t = trajectory.load(dcd_name, top=pdb_name)
  t.save(xtc_name)
  
## Execute the main function
if __name__ == "__main__":
  main()
