#!/usr/bin/python 

import sys
import argparse
import numpy as np

def parse_options():
    parser = argparse.ArgumentParser(description="Takes an SNP file"
				     "that has position, reference"
				     "allele, and alternative allele(s)"
				     "along with a pileup-file, and is able"
				     "to develop 'Starter-Position-Files'"
				     "to be input into bootstrapping programs"
				     "that can robustly test allele-specific "
				     "expression and translation. The SNP file"
				     "and pileup-file are for one chromosome"
                                     )
    parser.add_argument('--pileup_file', dest='pileup_file', type=argparse.FileType('r'),
                       help="A valid pileup_file that comes in this format: "
                       "chr1\tposition\tN\tNumBasePairs\tBasePairs"
                       )
    parser.add_argument('--snp_file', dest='snp_file', type=argparse.FileType('r'),
                       help="A valid snp_file that comes in this format: "
                       "Position\tRef\tAlt\tHaplotype")
    parser.add_argument('--read_length', '-r', dest='read_length', type=int,
                       help="A valid integer that is the read length of a read or "
                       "any integer greater than zero which represents the number of lines in the starter position files you would like.")
    parser.add_argument('--unit', '-u', dest='chromosome', type=str, 
                        help="A valid string that either represents the name of the chromosome or the gene for which you are generating files")
    parser.add_argument('--min_reads', '-m', dest='min_reads', type = int,
                        help="A valid integer that marks the minimum number of unskipped base pairs in the pileup for a particular SNP position")
    parser.add_argument('--output_dir', dest='output_dir', default=".", type=str,
                        help= "The output directory you'd like to have")
    options = parser.parse_args()

    return options


def add_pop_dictionary(pileup_file, pileup_dictionary):
    line = pileup_file.readline()
    if len(line) < 1: 
        pileup_dictionary = 0
        keys = 0
        return pileup_dictionary, keys
    thelist = line.split('\t')
    pileup_dictionary[thelist[1]] = thelist[4]
    keys = list(pileup_dictionary.keys())
    keys = list(map(int, keys))
    pop = str(min(keys))
    pileup_dictionary.pop(pop)
    return pileup_dictionary, keys


def fill_dict(pileup_file, read_length, pileup_dictionary):
    for x in range(0, read_length):
        line = pileup_file.readline()
        print(line)
        if len(line) < 1:
            print('Not enough lines in pileup file')
            pileup_dictionary = 0
            keys = 0
            return pileup_dictionary, keys
        thelist = line.split('\t')
        pileup_dictionary[thelist[1]] = thelist[4]
        keys = list(pileup_dictionary.keys())
        keys = list(map(int, keys))
        
    return pileup_dictionary, keys


def initialize_or_update(snp_file):
    line = options.snp_file.readline()
    if len(line) < 1:
        print('Not enough lines in snp file')
        snp_pos = 0
        hap = 0
        ref = 0
        alt = 0
        return snp_pos, hap, ref, alt
    thelist = line.split('\t')
    snp_pos = int(thelist[0])
    hap = 'NA\n'
    if len(thelist) > 3:
        hap = thelist[3]
    ref = thelist[1]
    alt = thelist[2]
    return snp_pos, hap, ref, alt



def parse_ref_alt(ref, alt, current_position_bp):
    ref_count = len(current_position_bp[np.where(current_position_bp ==ref)])
    alt_1_count = 0
    alt_2_count = 0
    alt_3_count = 0
    
    if len(alt) > 1:
        alt_np = alt.split(',')
        alt_1_count = len(current_position_bp[np.where(current_position_bp == alt_np[0])])
        if len(alt_np) > 1:
            alt_2_count = len(current_position_bp[np.where(current_position_bp == alt_np[1])])
        if len(alt_np) == 3:
            alt_3_count = len(current_position_bp[np.where(current_position_bp == alt_np[2])])
    else:
        alt_1_count = len(current_position_bp[np.where(current_position_bp == alt)])

    return ref_count, alt_1_count, alt_2_count, alt_3_count


def create_output(pileup_dictionary, ref, alt, hap, snp_pos, chromosome, read_length, min_reads, output_directory):
    
    if hap == '1|1' or hap == '0|0':
        return
    filename = chromosome + "." + str(snp_pos)
    basepair_and_numbases = pileup_dictionary[str(snp_pos)]
    truncated_string= ''.join([i for i in basepair_and_numbases if i.isalpha()])

    if len(truncated_string) < min_reads:
        print('Reads at ' + str(snp_pos) + ' do not mean minimum requirement\n')
        return
    if len(set(filter(str.isalpha,truncated_string)))==1:
        print('Homozygous at supposedly heterozygous position at ' + str(snp_pos) + "\n")
        return


    if len(ref) > 1:
        print('Invalid Reference Allele(s) at ' + str(snp_pos) + "\n")
        return
    if ref == '-':
        print('Invalid Reference Allele(s) at ' + str(snp_pos) + "\n")
        return
    if len(alt.split(',')[0]) > 1:
        print(alt.split(',')[0])
        print('Invalid Alternate Allele(s) at ' + str(snp_pos) + "\n")
        return
    
    np_allbasepairs = basepair_and_numbases.upper()
    np_allbasepairs = np.array(list(np_allbasepairs))

    numstart = len(np.where(np_allbasepairs == "^")[0])
    index = numstart*3+1

    current_positions_bp = np_allbasepairs[-index:]
    np_allbasepairs = np_allbasepairs[:-index]
    dollar_index = np.where(np_allbasepairs == "$")[0]
    permanent_basepairs = np.delete(np_allbasepairs, dollar_index)
    ref_count, alt_1_count, alt_2_count, alt_3_count = parse_ref_alt(ref, alt, current_positions_bp)
    filename = output_directory + "/" + filename
    outfile = open(filename, "w")
    outfile.write("SNP_POS\tREF\tAALT\tREF_COUNT\tALT_1_COUNT\tALT_2_COUNT\tALT_3_COUNT\tHAP\n")
    outfile.write(str(snp_pos)+ "\t"+ ref + "\t" + alt+ "\t"+ str(ref_count) +"\t" +str(alt_1_count)+ "\t" +str(alt_2_count)+ "\t"+ str(alt_3_count) + "\t"+hap)

    for x in range(1, read_length):
        snp_pos = snp_pos - 1
        if str(snp_pos) in pileup_dictionary:
            basepair_and_numbases = pileup_dictionary[str(snp_pos)].split('\t')
            np_allbasepairs = np.array(list(basepair_and_numbases[0]))
            numstart = len(np.where(np_allbasepairs == "^")[0])
            if numstart == 0:
                outfile.write(str(snp_pos) + "\t" + ref + "\t" + alt + "\t0\t0\t0\t0\t" + hap)
                continue
            current_positions_bp = permanent_basepairs[-numstart:]
            permanent_basepairs = permanent_basepairs[:-numstart]
            ref_count, alt_1_count, alt_2_count, alt_3_count = parse_ref_alt(ref, alt, current_positions_bp)
            outfile.write(str(snp_pos)+"\t"+ref+"\t"+alt+"\t"+ str(ref_count)+"\t"+str(alt_1_count)+"\t"+str(alt_2_count)+"\t"+str(alt_3_count)+"\t"+hap)
    
    outfile.close()

    


def main(options):

    options = parse_options()
    pileup_dictionary = dict()
    keys = []
    pileup_dictionary, keys = fill_dict(options.pileup_file, options.read_length, pileup_dictionary)
    snp_pos, hap, ref, alt = initialize_or_update(options.snp_file)

    while True:
        if pileup_dictionary == 0 or snp_pos == 0: 
            break
    
        if snp_pos not in pileup_dictionary and snp_pos > max(keys):
            pileup_dictionary, keys = add_pop_dictionary(options.pileup_file, pileup_dictionary)
            if pileup_dictionary == 0:
                break
        
        if snp_pos not in pileup_dictionary and snp_pos < max(keys):
            snp_pos, hap, ref, alt = initialize_or_update(options.snp_file)
            if snp_pos == 0:
                break

        if str(snp_pos) in pileup_dictionary:
            create_output(pileup_dictionary, ref, alt, hap, snp_pos, options.chromosome, options.read_length, options.min_reads)
            pileup_dictionary, keys = add_pop_dictionary(options.pileup_file, pileup_dictionary)
            snp_pos, hap, ref, alt = initialize_or_update(options.snp_file)
            if pileup_dictionary == 0 or snp_pos == 0:
                break
    

if __name__ == '__main__':
    options = parse_options()
    main(options)

