Please add brief descriptions for every script you add to this directory.


fasta_to_phylip.py: Converts a standard fasta file to phylip sequential format (e.g. from a MAFFT alignment to a suitable input file for RAxML).
NOTE: This script does not trim sequence identifiers. Phylip prefers identifiers
of less than 10 characters.

parse_FB_output.py: parses the output of FUBAR. Prints the site count (starting with one), then the omega at that site as calculated from alpha and beta.

remove_duplicates.py: removes duplicates from a set of sequences. Removes both identical sequences and identical FASTA identifiers.
