#!/usr/bin/python 

import sys
import argparse
import numpy as np
import pileuptable

##########################################################################
##########################################################################
##########################################################################
# This program takes in an SNP file and a PILEUP file and pumps out a ####
# 'Starter Position File' for each SNP in the SNP file that reaches the ##
# minimum read count. Each SNP's 'Starter Position File', will contain a #
# a list of x number of consecutive positions starting from the position #
# that is SNP - read_length all the way to the SNP. At each position #####
# leading up to the SNP of interest, we will have a count of how many ####
# reads from each haplotype began exactly at that position. Please refer##
# to the README on the github page for more information on what a ########
# 'Starter Position File' looks like. ####################################
##########################################################################
##########################################################################
##########################################################################

def main(options):

    options = parse_options()

    #######################################################################
    #### the pileuptable is a class derived from pileuptable.py that ######
    #### keeps track of SNPs from the SNP file and positions from the #####
    #### pileup file. See class for details. ##############################
    #######################################################################

    the_snp_pileup_table = pileuptable.SNP_PILEUP_TABLE( options.read_length)

    # Initialize our snp_pileup_table class by reading in the first how ever many lines ##
    # from the pileup file, and storing the necessary information ######

    the_snp_pileup_table.fill(options.pileup_file, options.read_length)

    # get our current snp # 
    the_snp_pileup_table.update_snp(options.snp_file)

    while True:

       ## if we've come to a point where either there are no pileup lines or 
       ## snps we stop and break.
        if the_snp_pileup_table.number == 0: 
            break
    
       ### if the snp isn't in our current list, but the SNP position is greater than the maximum we need to read more lines from the pileup table ####
        if the_snp_pileup_table.snp_pos not in the_snp_pileup_table.nppos and the_snp_pileup_table.snp_pos > max(the_snp_pileup_table.nppos):
            the_snp_pileup_table.update_pileup_table(options.pileup_file)
            
        ##### if the snp isn't in our current list but our pileup table starts with positions after the SNP, we need to update our SNP #####
        if the_snp_pileup_table.snp_pos not in the_snp_pileup_table.nppos and the_snp_pileup_table.snp_pos < min(the_snp_pileup_table.nppos):
            the_snp_pileup_table.update_snp(options.snp_file)
 
        ##### if the snp is in our current list, but it doesn't have enough positions before it, then we need to update the SNP #####
        if the_snp_pileup_table.snp_pos in the_snp_pileup_table.nppos and the_snp_pileup_table.snp_pos_minus_read_length not in the_snp_pileup_table.nppos:
            the_snp_pileup_table.update_snp(options.snp_file)
        
        ##### if the snp is in our current list, and all other metrics check out, we proceed to writing its starter position file ###### 
        if the_snp_pileup_table.snp_pos in the_snp_pileup_table.nppos:
            the_snp_pileup_table.create_output(options.chromosome,  options.min_reads, options.output_dir)
            the_snp_pileup_table.update_pileup_table(options.pileup_file)
            the_snp_pileup_table.update_snp(options.snp_file)

if __name__ == '__main__':
    options = parse_options()
    main(options)

