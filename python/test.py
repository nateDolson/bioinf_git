################################################################################
## Script for parsing SAM file for input creating a csv file for statistical
## Analysis of sequencing errors in R
##
## Created by: Nate Olson
## October 18, 2012
##
## parsing based on SAM format specifications v1.4-r985, September 7, 2011
################################################################################

#importing required modules
import re
from Bio import SeqIO

# reading sam input file
sam = open('/Users/nathandavidolson/Desktop/SequencingErrors/SRR072238_10_10_2012.sam', 'r')

#reading reference fasta file
ref = {}
for seq_record in SeqIO.parse("/Users/nathandavidolson/Desktop/SequencingErrors/ECandLM_HMP_mock_16SRef.fasta","fasta"):
    ref[seq_record.id] = seq_record.seq

# creating output file
samTOcsv = open('/Users/nathandavidolson/Documents/SAMparse_10_28_2012.csv','w')

# write to file column headers
samTOcsv.write('read_name,refseq_name,ref_value, map_qual,cigar_value,base_call,quality_score,read_position, ref_position\n')


# iterating through each line of the sam file
for line in sam:
    if not line.startswith("@"): # skip elements that start with @
        lineSplit = line.split('\t')
        # skip elements with a flag score (second element in the array) of 4 (reads not mapped)
        if not lineSplit[2] == 4 and int(lineSplit[4]) > 100:
            # expand cigar (element 5 of sub arrays)
            # array of letters in the cigar
            letters = re.findall('[A-Z]',lineSplit[5])
            # array of numbers in the cigar as intigers 
            numbers = [int(x) for x in re.findall(r'\d+', lineSplit[5])]
            expCigar = ''
            for i in range(0,len(letters)):
                expCigar = expCigar + numbers[i]*letters[i]

            # writing to output file individual line for each base
            # split cigar, sequence, and qual into arrays
            cigar = list(expCigar) # cigar
            sequence = list(lineSplit[9]) # sequence
            qual = list(lineSplit[10]) # qual
            
            #reference sequence
            ref_seq = ref[lineSplit[2]]
            # interate through each position of each of the reads creating a line for each
            # setting initial ref position value
            # for offsetting sequence and qual looping position for deletion in the read
            seq_qual_pos = 0
            ref_pos = int(lineSplit[3]) - 1
            
            #iterating through the cigar
            for i in range(0,len(cigar)):
                # setting values and increasing counts based on cigar value
                if cigar[i] == 'D':
                    # setting values
                    if seq_value % 1 != 0:
                        seq_value = cigar[i] - 0.95
                        qual_value = cigar[i] - 0.95
                    else:
                        seq_value = seq_value + 0.05
                        qual_value = qual_value + 0.05
                    ref_value = ref_seq[ref_pos]
                    ref_position = ref_pos

                    # increasing counts
                    ref_pos = ref_pos + 1 
                	               	
                elif cigar[i] == "I":
                	# setting values
                	seq_value = sequence[seq_qual_pos]
                	qual_value = ord(qual[seq_qual_pos]) - 33
                	ref_value = '*'
                	if ref_position % 1 != 0: #incrementing by 0.05 to incorporate insertions, approach only works for 20 instertions
                	    ref_position = ref_position + 0.05
                	else:
                	    ref_position = ref_pos - 0.95
			    # increasing counts
                            seq_qual_pos = seq_qual_pos + 1
                	
                elif cigar[i] == "M":
                    # setting values 
                    seq_value = sequence[seq_qual_pos]
                    qual_value = ord(qual[seq_qual_pos]) - 33
                    ref_value = ref_seq[ref_pos]
                    ref_position = ref_pos
                    # increasing counts
                    ref_pos = ref_pos + 1
                    seq_qual_pos = seq_qual_pos + 1
                	
                elif cigar[i] == "S":
                    seq_qual_pos = seq_qual_pos + 1
                    continue
                
                else:
                    continue
                
                # writing values to output file     	
                if not ref_value == seq_value and seq_qual_pos < 244: #cut off set at 244 based on the trim for Kunin et al 2010, only printing error positions
                	samTOcsv.write(lineSplit[0] + ',' + lineSplit[2]  + ',' + ref_value + ',' + lineSplit[4] + ',' +  cigar[i] + ',' + seq_value + ',' + str(qual_value) + ',' + str(seq_qual_pos) + ',' + str(ref_position) + '\n')
                int(seq_qual_pos)
        
samTOcsv.close()
sam.close()


''' Notes for cleaning up code, stops at SRR.....40553- the qual vlaue is out of reange
this is the last sequence in the file.  There is no quality string for the sequence?'''
''' M in the cigar can be match or mismatch
    Need to incorporate the reference base-
    open refernce fasta sequence file
    create a dictionary with the key at hte sequence name and the value the sequence as an array
    add the reference base to the csv using the refernece position - 1 for the array position

'''
