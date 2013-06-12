## This script takes command line arguments in the following form                ##
## python translate_align_mapping.py fasta_file pdb_structure_file output_type   ##

## Here is an example run to output the aligned and reverse translated sequences ##
## python translate_align_mapping.py fasta.fasta structure.pdb seqs              ##

## Here is an example run to return the sequence-structure map                   ##
## python translate_align_mapping.py fasta.fasta structure.pdb map               ##

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

  filename = remove_duplicate_sequences(filename)
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
## This function is designed to scan for and remove duplicate       ##
## sequences that may occur in the data, leaving only the first     ##
## one that is come to in the file.  Although the two sequence      ##
## may come from different sources, it is not useful                ##
## to have duplicates.                                              ##
######################################################################

def remove_duplicate_sequences(filename):
    print("\n\nStarting Remove Duplicate Sequence Function")
    input_handle = open(filename, 'rU')
    records = list(SeqIO.parse(input_handle, "fasta"))
    input_handle.close()

    outfile = 'tmp_remove_duplicates.txt'
    output_handle = open(outfile, 'w')

    count = 0
    a = len(records)

    for i in records:
        for j in records:
            if i.id != j.id:
                if str(i.seq) == str(j.seq):
                    records.remove(j)

        count = count + 1
        line = "%i%s" % (a-count, " Remaining")
        print line
	
    count = 0
    short_sequences = []
    b = len(records)
	
    for k in records:
        short_sequences.append(k)
        count  = count + 1			
        line = "%i" % (b-count)
        #print line
			
    SeqIO.write(short_sequences, output_handle, "fasta")
    output_handle.close()
    return outfile

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
        #print i.id
        
        while j < len(i.seq)-2:
            codons = "%c%c%c" % (i.seq[j], i.seq[j+1], i.seq[j+2])
            element = [codons, translate(codons)]
            j+=3
            if element[1] == "*" and j > len(i.seq)-2:
            	del element
            	continue
            All_codons.append(element)
            del element
        if include_seq:
        	MatchingDict[i.id] = All_codons
        	MatchingIDs.append(i.id)
        	l+=1

    first_line = "%c%s\n%s\n" % (">", "ref_seq", chain)

    outfile2 = 'tmp_preAlign.txt'
    output_handle2 = open(outfile2, 'w')
    output_handle2.write(first_line)
    
    for an_id in MatchingIDs:
        seq = ""
        line_aa = ""
        for letter in range (0, len(MatchingDict[an_id]) ):
        	seq = seq + (MatchingDict[an_id])[letter][1]
        line_aa = "%c%s\n%s\n" % (">", an_id, seq)
        output_handle2.write(line_aa)
        
    output_handle2.close()

    outfile3 = 'tmp_Aligned.txt'
    command = "muscle -in " + outfile2 + " -clwstrict -out " + outfile3 + " -diags -maxiters 2 -sv -distance1 kbit20_3"
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

      seq_counter=0
      
      for a_codon in range( 0, len(ref_seq) ):
        if a_record.seq[a_codon] == "-" and ref_seq[a_codon] != "-":
          MatchingDict[a_record.id].insert( a_codon-seq_counter+1, ['---', '-'] )
          continue
        elif a_record.seq[a_codon] == "-" and ref_seq[a_codon] == "-":
          this_sequence.pop(a_codon-seq_counter)
          seq_counter+=1
          continue
        elif ref_seq[a_codon]=="-" and a_record.seq[a_codon]!="-":
          seq_codon = str(this_sequence.pop( a_codon-seq_counter ))
          dict_codon = MatchingDict[a_record.id].pop( a_codon-seq_counter )
          assert seq_codon == dict_codon[1], \
            "\nThis is the dictionary codon: %s This is the sequence codon: %s\nHere is the seq id: %s Here is the location in the sequence: %s" \
            % (dict_codon, seq_codon, a_record.id, a_codon-seq_counter)
          seq_counter+=1
          continue
        else:
          continue

      assert len(this_sequence) == len(MatchingDict[a_record.id]), "\nLength mod seq: %s Length of Ref Seq: %s\n" \
        % (len(this_sequence), len(MatchingDict[a_record.id]))
            
      a_record.seq=this_sequence.toseq()
    				
    ref_seq = remove_values_from_list(ref_seq, "-")
    
    return(ref_seq, aligned_records, MatchingDict)
    
def output_nucleotide_sequences(ref_seq, aligned_records, MatchingDict):

  outfile = 'Final_Aligned.txt'
  output_handle = open(outfile, 'w')
  
  for a_record in aligned_records:
 
    seq = ""
    nuc_seq = ""
    aa_seq = ""
      
    for j in range(0, len(ref_seq)):
      assert ref_seq[j]!="-"
      nuc_seq = nuc_seq + (MatchingDict[a_record.id])[j][0]
    
    line = "%c%s\n%s\n" % (">", a_record.id, nuc_seq)
    output_handle.write(line)

  print("\n\nThe output file is name: " + outfile + "\n")
  subprocess.call('rm tmp_*', shell=True)
  return 0

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
      
      if(codon == '---'):
        continue
      elif(len(nuc_seq) == 0):
        nuc_seq = str(j)
      else:
        nuc_seq = nuc_seq + "," + str(j)
    
    line = "%c%s\n%s\n" % (">", a_record.id, nuc_seq)
    output_handle.write(line)

  print("\n\nThe output file is name: " + outfile + "\n")
  subprocess.call('rm tmp_*', shell=True)
  return 0

## Execute the main function
if __name__ == "__main__":
    main()

