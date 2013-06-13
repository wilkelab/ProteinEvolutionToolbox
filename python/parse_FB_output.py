from decimal import *
def find_omega(file_handle):
    outfile = "omegaFB.txt"
    out_handle = open(outfile, "w")
    count = 0
    out_handle.write("Site\tOmega\n")
    for a_line in file_handle.readlines():
        splits = a_line.split(",")
        codon = splits[0]
        alpha = splits[1]
        beta = splits[2]
        
        ##print alpha
        ##print beta
        ##a = Decimal((alpha))
        ##b = Decimal((beta))
        if alpha.isalpha():
            continue
        else:
            a = float(alpha)
            b = float(beta)
            c = a / b
            out_handle.write(codon + "\t" + str(c)  + "\n")
            print codon

def main():
    infile = "tree.tre.fubar.csv"
    file_handle = open(infile, "r")
    find_omega(file_handle)
main()
