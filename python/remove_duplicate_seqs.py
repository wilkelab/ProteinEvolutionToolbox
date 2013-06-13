from Bio import SeqIO, SeqRecord
import sys, os

def collapse_unsorted_seqs(infile):
    unique = {}
    in_handle = open(infile, "r")
    records = list (SeqIO.parse(in_handle, "fasta"))
    records.sort(cmp=lambda x,y: cmp(len(y),len(x)))
    dupe_count = 0;
    for record in records:
        if str(record.seq) in unique:
            dupe_count+=1
        else:
            found = False
            
        if found == False:
                unique[str(record.seq)] = str(record.id)

        with open(infile + ".no_duplicates", "w") as collapsed_file:
            for key in sorted(unique, cmp=lambda x,y: cmp(len(y),len(x)), key=unique.get):
                ##Format the line and print it to the file as follows:
                ##>OriginalReadFastaIdentifier Read1Identifier Read2Identifier...ReadNIdentifer
                ##ATGAACGCCTAA(The sequence of the original read)
                my_line='>%s\n%s\n' %(str(unique[key]), str(key))
                collapsed_file.write(my_line)
    return 0

def main():
    infile = sys.argv[1]
    collapse_unsorted_seqs(infile)
main()
