#!/usr/bin/env python
from __future__ import division

__doc__ = """
SYNOPSIS

    pileup_parse.py input.pileup output.tsv [output filter]

DESCRIPTION

    Converts pileup file as input to a tab seperated
    file with individual bases as rows in the document
    for use in analysis of sequencing errors.

USAGE
    Input file must be labeled .pileup for script to function

    Large datsets must be filtered as output files have as many
    rows as bases in the pileup file
    
    output filter - ***need to incorporate this function***
        -ft Filter type - use to state whether to retain or exclude
                          specified filters
        -F filter - specific filters to keep or remove can be any
                    character in a pileup base or qual file, range of
                    base positions in relation to the reference as well
                    as any of the following:
                        location: start middle end
                        read direction: forward reverse

### add in this function
snp_summary     Parse samtools consensus pileup and output penetrance summary of snps
                        Tab delimited file with: seq_id, position, ref_base, consensus_base, fraction A, C, G, T
AUTHOR

        Nate Olson <nolson@nist.gov>

"""

__author__ = "Nate Olson"
__version__ = '0.1'

import sys
import csv

def main():
    if len(sys.argv) > 1:
        if sys.argv[1] == "help" or sys.argv[1] == "-h":
            print help()
        elif len(sys.argv) < 2:
            print("Input and Output files not specified" % help())
        else:
            if "pileup" in sys.argv[1]:
                input_file = sys.argv[1]
                output_file = sys.argv[2]
            elif "pileup" in sys.argv[2]:
                input_file = sys.argv[2]
                output_file = sys.argv[1]
            if len(sys.argv) > 2:
                filters = sys.argv[3:]
            else:
                filters = "none"
            expand(input_file, output_file, filters)
        
    else:
        print("No arguments were provided see documentation below:\n\n %s"  %help())

def help():
    return globals()['__doc__']

def expand(input_file, output_file, filters):
    print "reading files ..."

    pileup = csv.reader(open(input_file, 'rb'), delimiter = "\t")

    expanded_pileup = open(output_file, 'wb')

    wr = csv.writer(expanded_pileup   , quoting=csv.QUOTE_ALL, delimiter = "\t")

    print "expanding pileup ..."
    #interate through each line expanding the table
    for position in pileup:
        quals = map(lambda x: ord(x) - 33, list(position[-1]))
        bases = list(position[-2])
        for qual in quals:
            location = 'did not work'
            direction = 'did not work'
            base = 'did not work'
            call = 'did not work'
            # Determining base position in read
            # Note that the ascii character
            # indicating mapping quality at the start
            # of the read is not incorporated

            if bases[0] == '^':
                location = 'start'
                del(bases[0:2])
            elif len(bases) > 1:
                if bases[1] == '$':
                    location = 'end'
                    del(bases[1])
                else:
                    location = "middle"
            else:
                location = 'middle'

            # Determining read direction
            if bases[0] in ['.','A','C','G','T','N']:
                direction = "forward"
            else:
                direction = "reverse"

            # Setting base
            base = bases[0] #will be reset for indels
            
            #taking into consideration indels
            if bases[0] == '+': #insertion
                call = "insertion"
                base = "".join(bases[2: 2 + int(bases[1])])
                if base.isupper():
                    direction == "forward"
                del(bases[1: 2 + int(bases[1])])
                
            elif bases[0] == '-': #deletion
                call = "deletion"
                base = "".join(bases[2: 2 + int(bases[1])])
                if base.isupper():
                    direction == "forward"
                del(bases[1: 2 + int(bases[1])])
            elif bases[0] == "*":
                call = "deletion"
                base = bases[0]
            elif bases[0] in ['.',',']:
                call = "match"
            else:
                call = "mismatch"
            del(bases[0:1])
            mylist = position[0:-3] + [location] + [direction] + [call] + [base] + [qual] 
            wr.writerow(mylist) 

def snp_summary(args):
    '''
    Parse samtools consensus pileup and output penetrance summary of snps
        Tab delimited file with: seq_id, position, ref_base, consensus_base, fraction A, C, G, T
    '''
    usage = "Usage: %prog snp_summary pileup_variations_file"
    parser = OptionParser(usage=usage)
    (opts, args) = parser.parse_args(args)
    
    if len(args) < 2:
        print ("ERROR: Must specify a pileup varations file")
        parser.print_help()
        sys.exit()
    
    f = open(args[1])
    print 'seq_id\tpostion\tref base\tconsensus base\tcoverage\tcount A\tcount C\tcount G\tcount T\tPrimary Allele\tPrimary Allele Frequency\tSecondary Allele\tSecondary Allele Frequency'
    for line in f:
        var = Variation(line.strip().split("\t"))
        #bases = ('A', 'C', 'G', 'T')
        count = {}
        count['A'] = var.base_count('A')
        count['C'] = var.base_count('C')
        count['G'] = var.base_count('G')
        count['T'] = var.base_count('T')

        if (var.coverage != sum(count.values())):
            warnings.warn("Counts do not add up to coverage.  Sum of counts: %i  Coverage: %i\n" % (sum(count.values()), var.coverage))

        count_sorted = count.keys()
        count_sorted.sort(key=count.__getitem__, reverse=True)
        primary_allele = count_sorted[0]
        primary_allele_freq = count[count_sorted[0]] / var.coverage
        secondary_allele = count_sorted[1]
        secondary_allele_freq = count[count_sorted[1]] / var.coverage
        print "%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%s\t%f\t%s\t%f" % (var.seq_id,
                                                      var.position,
                                                      var.reference_base,
                                                      var.consensus_base,
                                                      var.coverage,
                                                      count['A'],
                                                      count['C'],
                                                      count['G'],
                                                      count['T'],
                                                      primary_allele,
                                                      primary_allele_freq,
                                                      secondary_allele,
                                                      secondary_allele_freq,)
class Variation:
    '''
    store data on genomic variations (snps and small indels)
    '''
    '''
    chromosome, 1-based coordinate, reference base, consensus base, consensus quality, SNP quality, maximum mapping quality, the number of reads covering the site, read bases and base qualities.
    chromosome, 1-based coordinate, a star, the genotype, consensus quality, SNP quality, RMS mapping quality, # covering reads, the first alllele, the second allele, # reads supporting the first allele, # reads supporting the second allele and # reads containing indels different from the top two alleles
    '''

    def __init__(self, pileup_data):
        '''
        Constructor
        '''
        self.seq_id = pileup_data[0]
        self.position = int(pileup_data[1])
        self.reference_base = pileup_data[2]
        self.consensus_base = pileup_data[3]
        self.consensus_quality = int(pileup_data[4])
        self.snp_quality = int(pileup_data[5])
        self.mapping_quality = int(pileup_data[6])
        self.coverage = int(pileup_data[7])
        
        if self.reference_base == '*':
            self.type = 'indel'
            self.read_bases = None
            self.read_qualities = None
            self.first_allele = pileup_data[8]
            self.second_allele = pileup_data[9]
            self.first_allele_reads = int(pileup_data[10])
            self.second_allele_reads = int(pileup_data[11])
            self.other_allele_reads = int(pileup_data[12])
            self.unknown1 = int(pileup_data[13])
            self.unknown2 = int(pileup_data[14])
        else:
            self.type = 'snv'
            self.read_bases = pileup_data[8]
            self.read_qualities = pileup_data[9]
            self.first_allele = None
            self.second_allele = None
            self.first_allele_reads = None
            self.second_allele_reads = None
            self.other_allele_reads = None
            
    def base_count(self, base):
        count = self.read_bases.count(base.lower()) + self.read_bases.count(base.upper())
        if (base.upper() == self.reference_base) or (base.lower() == self.reference_base):
            count += self.read_bases.count('.')
            count += self.read_bases.count(',')
        return count

    def __repr__(self):
        return repr((self.seq_id, self.position, self.reference_base, self.consensus_base, self.consensus_quality, self.snp_quality,
                     self.max_mapping_quality, self.coverage, self.read_bases, self.read_qualities))
    
    def key(self):
        return "%s%s%s%s" % (self.seq_id, str(self.position), self.reference_base, self.consensus_base)

    def pileup_line(self):
        if self.type == 'snv':
            return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.seq_id, str(self.position), self.reference_base, self.consensus_base,
                                                           str(self.consensus_quality), str(self.snp_quality), str(self.mapping_quality),
                                                           str(self.coverage), self.read_bases, self.read_qualities)
        elif self.type == 'indel':
            return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.seq_id, str(self.position), self.reference_base, self.consensus_base,
                                                           str(self.consensus_quality), str(self.snp_quality), str(self.mapping_quality),
                                                           str(self.coverage), self.first_allele, self.second_allele, self.first_allele_reads,
                                                           self.second_allele_reads, self.other_allele_reads, self.unknown1, self.unknown2)
        else:
            return ""
        #return '\t'.join([self.seq_id, str(self.position), self.reference_base, self.consensus_base, str(self.consensus_quality),
        #                  str(self.snp_quality), str(self.max_mapping_quality), str(self.coverage), self.read_bases, self.read_qualities])
            
                
if __name__ == '__main__':
    main()
    sys.exit()

'''
Notes on pileup formating 
. match in forward strand
, match in reverse strand
[ATCGN] mismatch in forward strand
[atcgn] mismatch in reverse strand
\+[0-9]+[ACGTNacgtn]+ insertion between this reference position and next
-[0-9]+[ACGTNacgtn]+ delection from the reference
^ marks the start of a read segment
'''
