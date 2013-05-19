# script for creating a list of sequences with specific taxonomy
# input required: output from QIIME assign_taxonomy.py

# name of output file from QIIME assign_taxonomy.py
taxonomyFile = "/Users/nathandavidolson/Desktop/rdp22_assigned_taxonomy/SRR072223_tax_assignments.txt"

# opening taxonomy file
taxa = open(taxonomyFile, 'r')

# creating output files for E. coli and L. monocytogenes
Ecoli = open("EcoliSeq.txt", 'w')
Lmono = open("LmonoSeq.txt", 'w')


# reading taxa file one line at a time
taxaLines = taxa.readlines()
for line in taxaLines:
    # determining if the sequence is for E. coli or L. mono
    # and printing the seqeunce name to a list
    if "Escherichia" in line:
        sline = line.split(" ")
        Ecoli.write(sline[0] + "\n")
    elif "Listeria" in line:
        sline = line.split(" ")
        Lmono.write(sline[0] + "\n")    
taxa.close()
Ecoli.close()
Lmono.close()
