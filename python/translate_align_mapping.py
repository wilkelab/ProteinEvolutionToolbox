###################################################################################
## This script was written by Austin Meyer in the lab of Dr. Claus Wilke at the  ##
## University of Texas at Austin.  For any questions or concerns you can email   ##
## Austin at austin.g.meyer@gmail.com.                                           ##
###################################################################################

###################################################################################
## This script takes command line arguments in the following form                ##
## python translate_align_mapping.py fasta_file pdb_structure_file output_type   ##
##                                                                               ##
## Here is an example run to output the aligned and reverse translated sequences ##
## python translate_align_mapping.py fasta.fasta structure.pdb seqs              ##
##                                                                               ##
## Here is an example run to return the sequence-structure map                   ##
## python translate_align_mapping.py fasta.fasta structure.pdb map               ##
###################################################################################

######################################################################
## Some of these includes are probably unecessary.                  ##
######################################################################
import sys, subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

######################################################################
## This is the dict for switching back and forth between three      ##
## letter and one letter amino acid definitions.  This probably     ##
## isn't necessary since it is probably implemented in biopython    ##
## somewhere.                                                       ##
######################################################################

resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }
            
######################################################################
## Here are the maximum RSA values per one letter amino acid        ##
## abbreviation.                                                    ##
######################################################################

residue_max_acc = {'A': 113.0, 'R': 241.0, 'N': 158.0, 'D': 151.0, \
                   'C': 140.0, 'Q': 189.0, 'E': 183.0, 'G': 85.0,  \
                   'H': 194.0, 'I': 182.0, 'L': 180.0, 'K': 211.0, \
                   'M': 204.0, 'F': 218.0, 'P': 143.0, 'S': 122.0, \
                   'T': 146.0, 'W': 259.0, 'Y': 229.0, 'V': 160.0}

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
  filename =  args[1]
  pdb      =  args[2]
  out_form =  args[3]

  structure = parsePDBStructure( pdb )
  (new_pdb, chain_renumbered_pdb) = renumberResidues( structure )

  (ref_seq, aligned_records, MatchingDict) = translate_align_revtranslate(filename, new_pdb)
  
  if(out_form == "seqs"):
    output_nucleotide_sequences(ref_seq, aligned_records, MatchingDict)
  elif(out_form == "map"):
    output_numerical_map(ref_seq, aligned_records, MatchingDict)
    
  return 0
  
######################################################################
## Removes all occurences of an element from a list without having  ##
## to sort the list.                                                ##
######################################################################

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]

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

######################################################################
## This function will return the amino acid sequence from the PDB   ##
## file.                                                            ##
######################################################################

def get_aa_fromPDB( pdb_id ):
  structure = parsePDBStructure( pdb_id )
  polypeptides = ""
  ppb=PPBuilder()
  for pp in ppb.build_peptides(structure):
    polypeptides += pp.get_sequence()
  return polypeptides

######################################################################
## This is sort of the crux function.  It is designed to take in a  ##
## file with a list of sequences and a pdb filename of the          ##
## structure corresponding to the sequences.  Then, the script will ##
## first translate the sequence and build a list of amino acids to  ##
## correspond to the codon list allowing it to remember the codon   ##
## for each site... I tried to build a dictionary instead, but a    ##
## list seemed to work just as well.  Then, the script will align   ##
## the amino acid sequence via a system call to muscle, and will    ##
## output the corresponding codons to the aligned sequences.        ##
##                                                                  ##
## More documentation here is probably good.                        ##
######################################################################

def translate_align_revtranslate( filename, pdb_id ):
    print("\n\nStarting Translate and Alignment Function\n\n")
    
    input_handle1 = open(filename, 'rU')
    records1 = list(SeqIO.parse(input_handle1, "fasta"))
    input_handle1.close()
   
    proteinSequences = []
    MatchingIDs = []
    MatchingDict = {}
    l = 0
    
    chain = get_aa_fromPDB( pdb_id )
    
    for i in records1:  
        j = 0
        All_codons = []
        include_seq=True
        
        while j < len(i.seq)-2:
            codons = "%c%c%c" % (i.seq[j], i.seq[j+1], i.seq[j+2])
            element = [codons, translate(codons)]
            j+=3
            if element[1] == "*" and j > len(i.seq)-2:
               del element
               continue
            elif element[1] == "*" or element[1] not in residue_max_acc:
               include_seq=False
               break
            All_codons.append(element)
            del element

        if include_seq:
          MatchingDict[i.id] = All_codons
          MatchingIDs.append(i.id)
          l+=1
          print l
        else:
          records1.remove(i)
          continue

    first_line = "%c%s\n%s\n" % (">", "ref_seq", chain)
    
    outfile2 = 'tmp_preAlign.txt'
    
    output_handle2 = open(outfile2, 'w')
    output_handle2.write(first_line)

    j = 0
    
    for an_id in MatchingIDs:
        seq = ""
        line_aa = ""
        for letter in range (0, len(MatchingDict[an_id]) ):
         seq = seq + (MatchingDict[an_id])[letter][1]
        line_aa = "%c%s\n%s\n" % (">", an_id, seq)
        output_handle2.write(line_aa)
        j+=1
    output_handle2.close()
    
    outfile3 = 'tmp_Aligned.txt'
    command = "mafft --clustalout " + outfile2 + " > " + outfile3
    subprocess.call(command, shell=True)
	
    input_handle2 = open(outfile3, 'rU')
    aligned_records = list(SeqIO.parse(input_handle2, "clustal"))
	
    ref_seq = ""
    for i in aligned_records:
        if i.id == "ref_seq":
            ref_seq = i.seq
            aligned_records.remove(i)
        else:
            continue
    
    ### Make a loop that looks for sequences with in/dels and removes the insertion relative to the reference sequence
    for a_record in aligned_records:
        
      ### Make sure we're not working with the reference sequence
      this_sequence = (a_record.seq).tomutable()
      
      for a_codon in range( 0, len(ref_seq) ):
        if a_record.seq[a_codon] == '-':
          MatchingDict[a_record.id].insert( a_codon, ['---', '-'] )
          
      seq_counter=0
      
      for a_codon in range( 0, len(ref_seq) ):
        if ref_seq[a_codon] == "-":
          dict_codon = MatchingDict[a_record.id].pop( a_codon - seq_counter )
          seq_counter+=1
          
      assert len(remove_values_from_list(ref_seq, "-")) == len(MatchingDict[a_record.id]), "\nLength mod seq: %s Length of Ref Seq: %s\n" \
        % (len(remove_values_from_list(ref_seq, "-")), len(MatchingDict[a_record.id]))
    				
    ref_seq = remove_values_from_list(ref_seq, "-")
    
    return(ref_seq, aligned_records, MatchingDict)

######################################################################
## This function is for outputting the final aligned nucleotide     ##
## sequences.  It also outputs the corresponding nucleotide file    ##
## temporarily for error checking.  If you want to preserve the tmp ##
## files you must comment out the second to last line of the        ##
## function.                                                        ##
######################################################################

def output_nucleotide_sequences(ref_seq, aligned_records, MatchingDict):

  outfile = 'Final_Aligned.txt'
  aa_outfile = 'tmp_aa_out.txt'

  output_handle = open(outfile, 'w')
  aa_out_handle = open(aa_outfile, 'w')

  for a_record in aligned_records:
 
    seq = ""
    nuc_seq = ""
    aa_seq = ""
      
    for j in range(0, len(ref_seq)):
      assert ref_seq[j]!="-"
      nuc_seq = nuc_seq + (MatchingDict[a_record.id])[j][0]
      aa_seq  = aa_seq  + (MatchingDict[a_record.id])[j][1]
    
    line = "%c%s\n%s\n" % (">", a_record.id, nuc_seq)
    output_handle.write(line)

    line = "%c%s\n%s\n" % (">", a_record.id, aa_seq)
    aa_out_handle.write(line)

  print("\n\nThe output file is name: " + outfile + "\n")
  subprocess.call('rm tmp_*', shell=True)
  return 0

######################################################################
## This function outputs a numerical map for each input nucleotide  ##
## sequence.  The sequence of comma delimited numbers provides the  ##
## corresponding site in the pdb's amino acid sequence.             ##
######################################################################

def output_numerical_map(ref_seq, aligned_records, MatchingDict):

  outfile = 'Final_Aligned.txt'
  output_handle = open(outfile, 'w')
  
  for a_record in aligned_records:
 
    seq = ""
    aa_seq = ""
    nuc_seq = ""
      
    for j in range(0, len(ref_seq)):
      assert ref_seq[j]!="-"
      codon = (MatchingDict[a_record.id])[j][0]
      
      if(len(nuc_seq) == 0):
        nuc_seq = str(j)
      else:
        nuc_seq = nuc_seq + "," + str(j)
    
    line = "%c%s\n%s\n" % (">", a_record.id, nuc_seq)
    output_handle.write(line)

  print("\n\nThe output file is name: " + outfile + "\n")
  # subprocess.call('rm tmp_*', shell=True)
  return 0

## Execute the main function
if __name__ == "__main__":
    main()

