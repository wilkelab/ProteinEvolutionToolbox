from Bio import AlignIO
import sys, os

def convert_seqs(infile):
    in_handle = open(infile, "rU")
    print infile
    out_handle = open(infile + ".phy", "w")
    alignment = AlignIO.parse(in_handle, "fasta")
    AlignIO.write(alignment, out_handle, "phylip")

    out_handle.close()
    in_handle.close()


def main():
    infile = sys.argv[1]
    convert_seqs(infile)
main()
