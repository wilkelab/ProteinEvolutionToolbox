#!/usr/bin/env python

from prody import *
from prody.utilities import which
from multiprocessing import Pool

import numpy as np
import os, sys, shutil, glob, subprocess


def main(): 
  if len( sys.argv ) > 4:
    print '''

    You don't have the right number of arguments.

    '''
  elif len( sys.argv ) ==  3:
    location = str(sys.argv[1])
    distance = str(sys.argv[2])

  files = get_minimal_filenames('.', '.pdb')

  print(files)

  location = '.'
  distance = '5'

  for pdb in files:
    parameterize and combine files
    parameterize_structures(pdb)

    combined_name = pdb + "_c_b.pdb"

  pdb_files = get_filenames(location, 'c_b.pdb')

  #Count contacts
  num_contacts = get_contacts_files(location, pdb_files, distance)
  print(num_contacts)

  write_to_file('contacts.txt', num_contacts, ',')

def write_model(filename):
  model_setup = 'package require autopsf\nmol load pdb ' + filename + '\nautopsf -protein ' + filename 

  output = open('model.tcl', 'w')
  output.write(model_setup)
  output.close()

def write_combination(receptor_name_pdb, receptor_name_psf, ligand_name_pdb, ligand_name_psf):
  inplace_change(receptor_name_pdb, ' P1', ' R1')
  inplace_change(receptor_name_pdb, ' P2', ' R2')
  inplace_change(receptor_name_pdb, ' P3', ' R3')
  inplace_change(receptor_name_pdb, ' P4', ' R4')
  inplace_change(receptor_name_psf, ' P1', ' R1')
  inplace_change(receptor_name_psf, ' P2', ' R2')
  inplace_change(receptor_name_psf, ' P3', ' R3')
  inplace_change(receptor_name_psf, ' P4', ' R4')

  inplace_change(ligand_name_pdb, ' P1', ' L1')
  inplace_change(ligand_name_pdb, ' P2', ' L2')
  inplace_change(ligand_name_pdb, ' P3', ' L3')
  inplace_change(ligand_name_pdb, ' P4', ' L4')
  inplace_change(ligand_name_psf, ' P1', ' L1')
  inplace_change(ligand_name_psf, ' P2', ' L2')
  inplace_change(ligand_name_psf, ' P3', ' L3')
  inplace_change(ligand_name_psf, ' P4', ' L4')

  model_setup = 'package require psfgen\nresetpsf\nreadpsf ' + receptor_name_psf + '\nreadpsf ' + ligand_name_psf + '\n'
  model_setup += 'coordpdb ' + receptor_name_pdb + '\ncoordpdb ' + ligand_name_pdb +'\n'
  model_setup += 'writepsf ' + receptor_name_pdb[:4] + '_c_b.psf\nwritepdb ' + receptor_name_pdb[:4] + '_c_b.pdb'

  output = open('combine.tcl', 'w')
  output.write(model_setup)
  output.close()

def parameterize_structures(file_prefix):
  receptor_name = file_prefix + '_r_b.pdb'
  ligand_name = file_prefix + '_l_b.pdb'

  write_model(receptor_name)
  execute_line = "vmd -dispdev text -eofexit < model.tcl"
  subprocess.call(execute_line, shell=True)
  
  write_model(ligand_name)
  execute_line = "vmd -dispdev text -eofexit < model.tcl"
  subprocess.call(execute_line, shell=True)

  receptor_name_pdb = file_prefix + '_r_b_autopsf.pdb'
  receptor_name_psf = file_prefix + '_r_b_autopsf.psf'
  ligand_name_pdb = file_prefix + '_l_b_autopsf.pdb'
  ligand_name_psf = file_prefix + '_l_b_autopsf.psf'
  
  write_combination(receptor_name_pdb, receptor_name_psf, ligand_name_pdb, ligand_name_psf)

  write_model(ligand_name)
  execute_line = "vmd -dispdev text -eofexit < combine.tcl"
  subprocess.call(execute_line, shell=True)

  subprocess.call('rm *autopsf*', shell=True)

def anm_conf_sampling(name, pdb):
  #Set up the ANM
  structure_anm = ANM('structure')
  structure_anm.buildHessian(pdb.ca)
  structure_anm.calcModes(n_modes=3)
  structure_anm_ext, structure_all = extendModel(structure_anm, pdb.ca, pdb, norm=True)

  #Conformational sampling
  ens = sampleModes(structure_anm_ext, atoms=pdb.protein, n_confs=4, rmsd=1.0)  
  pdb.addCoordset(ens.getCoordsets())
  pdb.all.setBetas(0)
  pdb.ca.setBetas(1)

  dir_name = name[0:4] + '_ens'

  if os.path.exists(dir_name):
    shutil.rmtree(dir_name)
  
  os.makedirs(dir_name)
  
  for i in range(1, pdb.numCoordsets()):
    fn = os.path.join(dir_name, name[0:4] + '_' + str(i) + '.pdb') 
    writePDB(fn, pdb, csets=i)

  return(0)

def conf_opt_setup(name):
  namd2 = which('namd2')
  par = os.path.join('/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmpar1.2', 'par_all27_prot_lipid_na.inp')

  dir_name = name[0:4] + '_opt'

  if os.path.exists(dir_name):
    shutil.rmtree(dir_name)

  os.makedirs(dir_name)

  conf = open('min.conf').read()

  for pdb in glob.glob(os.path.join(name[0:4] + '_ens', '*.pdb')):
    fn = os.path.splitext(os.path.split(pdb)[1])[0]
    pdb = os.path.join('..', pdb)
    out = open(os.path.join(dir_name, fn + '.conf'), 'w')
    out.write(conf.format(
      out=fn, pdb=pdb, par=par))
    out.close()

def conf_opt(name):
  os.chdir(name + '_opt')

  cmds=[]

  for conf in glob.glob('*.conf'):
    inplace_change(conf, 'structure       ../structure', 'structure       ../' + name + '_c_b.psf')
    fn = os.path.splitext(conf)[0]
    cmds.append('namd2 ' + conf + ' > ' + fn + '.log')

  pool=Pool(2)

  signals = pool.map(os.system, cmds)
  
  os.chdir('../')

  return(0)
 
def rename_files(location, files, new_extension):
  os.chdir(location)
  for i in files:
    shutil.copy(i, str(i[:-5]) + new_extension)
 
  os.chdir('..')
  return(0)

def get_minimal_filenames(location, extension):
  raw_files = os.listdir(location)
 
  files = []
  for i in raw_files:
    if i.endswith(extension):
      files.append(i[:4])

  return(sorted(list(set(files))))

def get_filenames(location, extension):
  raw_files = os.listdir(location)
 
  files = []
  for i in raw_files:
    if i.endswith(extension):
      files.append(i)

  return(sorted(files))

def get_contacts_files(location, files, distance):

  num_contacts=[]

  print(files)

  for i in range(0,len(files)):
    structure=parsePDB(files[i])
    names = list(set(structure.getSegnames()))
    
    receptor_segnames = 'segment'
    ligand_segnames = 'segment'

    for i in names:
      if i.startswith( 'R' ):
        receptor_segnames += ' ' + i
      elif i.startswith( 'L' ):
        ligand_segnames += ' ' + i
    
    receptor = structure.select(receptor_segnames + " and noh").copy()
    ligand = structure.select(ligand_segnames + " and noh").copy()
    contacts=receptor.select('calpha and (same residue as within ' + distance + ' of ligand)', ligand=ligand)
    num_contacts.append(len(contacts))

  return(np.array(num_contacts))

def write_to_file(filename, data, delimiter):
  output = open(filename, 'w')

  for i in range(0, len(data)):
    if i != len(data)-1:
      output.write(str(data[i]) + delimiter)
    else:
      output.write(str(data[i]))

  print(data)

  output.close()

def inplace_change(filename, old_string, new_string):
        s=open(filename).read()
        if old_string in s:
                print 'Changing "{old_string}" to "{new_string}"'.format(**locals())
                s=s.replace(old_string, new_string)
                f=open(filename, 'w')
                f.write(s)
                f.flush()
                f.close()
        else:
                print 'No occurances of "{old_string}" found.'.format(**locals())

## Execute the main function
if __name__ == "__main__":
  main()

