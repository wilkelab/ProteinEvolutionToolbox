""" 
    SJS 9/24/14.
    Takes in a file of nucleotide sequences, translates, aligns, and returns both a nucleotide and protein alignment.
    Default aligner is mafft --auto. Output files will be in FASTA format. If you want something else, you're on your own.
    Usage: python align.py <infile> [optional arguments]. To see optional arguments and full instructions, enter python align.py -h/--help .
"""

import os
import sys
import argparse
import subprocess
from Bio import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna


      
def parse_args():
    parser = argparse.ArgumentParser(description = 'Takes a file of nucleotide sequences and returns a protein alignment made w/ mafft, as well as nucleotide alignment back-translated from this protein alignment.')
    parser.add_argument("infile", metavar = "infile",  type = str, help = "A file containing nucleotide sequences to align. Must be in fasta, phylip, or nexus format.")
    parser.add_argument("-prot_outfile", dest = "prot_outfile", default = "protein.fasta", type = str, help = "Output file for protein alignment.")
    parser.add_argument("-nuc_outfile", dest = "nuc_outfile", default = "nucleotide.fasta", type = str, help = "Output file for nucleotide alignment.")
    parser.add_argument("-aligner", dest = "aln_exec", default = "mafft --auto", type = str, help="Any mafft algorithm (mafft, linsi, ginsi, einsi...). Default is 'mafft --auto'.")
    parser.add_argument("-options", dest="aln_opt", default = " --quiet ", type = str, help = "Optional mafft parameters.")
    
    return parser.parse_args()



def parse_translate_infile(infile):
    ''' Read in input sequence file, which can be in either fasta, phylip, or nexus.
        Translate nucleotide sequences to protein sequences.
        Returns two dictionaries - one with nucleotide seqs, one with protein seqs. All unaligned. 
    '''

    record_dict = {}
    try:
        records = list(SeqIO.parse(infile, 'fasta'))
    except:
        try:
            records = list(SeqIO.parse(infile, 'phylip-relaxed'))
        except:
            try:
                records = list(SeqIO.parse(infile, 'nexus'))
            except:
                raise AssertionError("\n\nCannot read in input sequence file. Verify that it is either in FASTA, phylip, or nexus format.")
    for rec in records:
        record_dict[str(rec.id)] = str(rec.seq)
    return translate_records(record_dict)
    
    

def translate_records(nuc_dict):
    '''  If stop codon at last codon, remove it from the nucleotide sequence.'''
    
    prot_dict = {}
    for entry in nuc_dict:
        nuc_dict[entry] = nuc_dict[entry].replace('-','') # Remove any possible gaps in nucleotide sequence
        nucseq = nuc_dict[entry]
        assert(len(nucseq)%3 == 0), "Nucleotide sequence length is not a multiple of three. I cannot translate, so I'm quitting."
        
        prot_seq = ''
        for i in range(0,len(nucseq),3):
            codon = nucseq[i:i+3]
            try:
                amino = str( Seq.Seq(codon, generic_dna).translate() )
            except:
                raise AssertionError("\n\nCould not translate input nucleotide codon, quitting.")
            if amino == '*':
                # If stop codon is the last codon, just remove it
                if i == len(nucseq)-3:
                    nuc_dict[entry] = nuc_dict[entry][:-3]
                else:
                    raise AssertionError("\n\n You have internal stop codons, quitting.")
            else:
                prot_seq += amino
        prot_dict[entry] = prot_seq
    
    return nuc_dict, prot_dict




def align_seqs(exec_, options, prot_dict, infile, outfile):
    with open(infile, 'w') as inf:
        for rec in prot_dict:
            inf.write('>' + str(rec) + '\n' + str(prot_dict[rec]) + '\n')
    run_align = subprocess.call( exec_ + ' ' + options + ' ' + infile + ' > ' + outfile, shell=True)
    assert(run_align == 0), "Could not perform mafft alignment."
    return list(SeqIO.parse(outfile, 'fasta'))
    
 
def back_translate(protseq, nucseq):
    '''Back translate an individual sequence''' 
    nucaln = ''
    start = 0; end = 3;
    for amino in protseq:
        if amino == '-':
            codon = '---'
        else:
            codon = nucseq[start:end]
            start += 3
            end += 3
        nucaln += codon
    assert(len(protseq)*3 == len(nucaln)), "Back-translation failed."
    return nucaln
  

def pal_to_nal(aln_records, nuc_dict, prot_outfile, nuc_outfile):
    protf = open(prot_outfile, 'w')
    nucf  = open(nuc_outfile,  'w')
    
    for protrec in aln_records:
        id = str(protrec.id)
        aln_nuc = back_translate( str(protrec.seq), nuc_dict[id] )
        protf.write('>' + id + '\n' + str(protrec.seq) + '\n')
        nucf.write('>' + id + '\n' + aln_nuc + '\n')
    protf.close()
    nucf.close()
    
   
def main():
    args = parse_args()
    while args.infile is None:
        args.infile = raw_input("\nYou gotta specify an input file: ")
        if not os.path.exists(args.infile):
            args.infile = None

    # Read in sequences and translate to protein
    print("Reading input sequences")
    nuc_dict, prot_dict = parse_translate_infile(args.infile)

    # Align protein sequences
    print("Creating protein alignment")
    aln_records = align_seqs(args.aln_exec, args.aln_opt, prot_dict, 'in.fasta', 'out.fasta')

    # Back-translate protein alignment to nucleotide alignment
    print "Back-translating protein alignment to a nucleotide alignment"
    pal_to_nal(aln_records, nuc_dict, args.prot_outfile, args.nuc_outfile)
    
    # clean up temp files
    os.remove('out.fasta')
    os.remove('in.fasta')
    
    outstring = "\nComplete! Your final AA alignment is in " + str( args.prot_outfile) +  " and your final nuc alignment is in " +  str(args.nuc_outfile) +  "\n"
    print(outstring)

main()









