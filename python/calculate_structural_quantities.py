#!/usr/bin/python

import sys

from prody import *
from numpy import *
from Bio.Cluster import kcluster, clustercentroids
from numpy import *
from scipy.cluster import vq
from pylab import *

def loadFiles(pdb_name, dcd_name, start, stop, step):
  dcd = parseDCD(dcd_name, start = start, stop = stop, step = step)
  structure = parsePDB(pdb_name)

  dcd.setCoords(structure)
  dcd.setAtoms(structure.calpha)

  dcd.iterpose()

  return (dcd, structure)

def calculateClusters(dcd, structure):

  com = array(calcCenter(dcd.getCoordsets()))

  centers, dist = vq.kmeans(com, 3)

  savetxt('centers.txt', centers, delimiter=',')
  savetxt('centroids.txt', com, delimiter=',')
  
  return 0

def calculateEDA(dcd, structure):
  eda_ensemble = EDA('temp')
  eda_ensemble.buildCovariance(dcd)
  eda_ensemble.calcModes()

  conformation = dcd[0:]
  
  projection = calcProjection(conformation, eda_ensemble[1:6])
  sqflucts = calcSqFlucts(eda_ensemble[1:3])

  savetxt('projection.txt', projection, delimiter=',')
  savetxt('sqflucts.txt', sqflucts, delimiter=',')

  return 0

def calculateRMSs(dcd, structure):
  rmsfs = dcd.getRMSFs()
  rmsds = dcd.getRMSDs()

  savetxt('rmsf.txt', rmsfs, delimiter=',')
  savetxt('rmsd.txt', rmsds, delimiter=',')

  show(plot(rmsds))

  return 0

def generateAnalysis(pdb_file, dcd_file, out_file, start, stop, step):
  (dcd, structure) = loadFiles(pdb_file, dcd_file, start, stop, step)

  calculateClusters(dcd, structure)

  calculateEDA(dcd, structure)

  calculateRMSs(dcd, structure)

  return 0

def main():
  if len( sys.argv ) != 7:
    print '''

    You screwed up choosing input files.

    '''
    print "     ", sys.argv[0], "<pdb file> <dcd file> <start> <stop> <step> <output file>"
    
  else:
    pdb_file = sys.argv[1]
    dcd_file = sys.argv[2]
    start    = int(sys.argv[3])
    stop     = int(sys.argv[4])
    step     = int(sys.argv[5])
    out_file = sys.argv[6]

    generateAnalysis( pdb_file, dcd_file, out_file, start, stop, step )

if __name__ == "__main__":
  main()

