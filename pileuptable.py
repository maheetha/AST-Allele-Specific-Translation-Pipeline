import sys
import numpy as np
import bases


######### Parsing options ############
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
    parser.add_argument('--read_length', '-r', dest='read_length', default = 30, type=int,
                       help="A valid integer that is the read length of a read or "
                       "any integer greater than zero which represents the" 
                       "number of lines in the starter position files you would like.")
    parser.add_argument('--unit', '-u', dest='chromosome', type=str, 
                        help="A valid string that either represents the name of the chromosome or the gene for which you are generating files")
    parser.add_argument('--min_reads', '-m', dest='min_reads', default = 20, type = int,
                        help="A valid integer that marks the minimum number of unskipped base pairs in the pileup for a particular SNP position")
    parser.add_argument('--output_dir', dest='output_dir', default=".", type=str,
                        help= "The output directory you'd like to have")
    options = parser.parse_args()


    if not (options.snp_file or options.pileup_file):
        parser.error("expected --snp_file and --pileup_file "
                         "at least one of these files is missing")
    
    if not (options.unit):
        parser.error("you must give a name for the prefix of output files "
                        )
   
    return options

######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### 
######### The Class SNP_PILEUP_TABLE keeps track of many ######### 
######### variables that are necessary to create 'Starter ######### 
######### 'Position Files' for each SNP. ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### 

class SNP_PILEUP_TABLE:

	## Initializing the class object ###
    def __init__(self, read_length):

        self.initialize(read_length)

	#### signals to end the program ####
    def clear(self):
       
        self.number = 0

######### ######### ######### ######### ######### ######### ######### ######### 
######### The initialize method initializes all variables needed to keep ######### 
######### track of the relationship between the current snp of interest ######### 
######### and the positions and reads in the pileup files ######### ######### #####
######### ######### ######### ######### ######### ######### ######### ######### ###
######### ######### ######### ######### ######### ######### ######### ######### 


    def initialize(self, read_length):
        ### keeps track of readlength. This is important to make the arrays 
        ### that will keep track of the position and how many beginning reads
        self.readlength = read_length 

        ### array that keeps track of which position ###
        self.nppos = np.empty([self.readlength]) 

        ### array that keeps track of how many "starter reads" at each position #####
        self.npreads = np.empty([self.readlength]) 

        ### keeps track of the snp of interest
        self.snp_pos = 0 

	### keeps track of the snp_pos - 30, so we know we have enough reads to create a SPF
        self.snp_pos_minus_read_length = 0 

	### keeps track of haplotype, reference, alternative allele info for our SNP of interets
        self.hap = "" 
        self.ref = "" 
        self.alt1 = ""
        self.alt2 = ""
        self.alt3 = ""

	### everytime, we need to count the number of ref alleles ###
	### and alt alleles in a particular array of bases ####
        self.refcount = 0  
        self.alt1count = 0
        self.alt2count = 0
        self.alt3count = 0

	### keeps track of the current reads and current position
        self.currentreads = "" 
        self.currentpos = 0 

	### this signals that we still have SNPs and lines in the pileup table ###
        self.number = 1 


#################################################################################
#################################################################################
######### update_pileup_table simply reads another line from the pileup #########
######### file, parses the fields, and updates various variables in #########
######### our class. ######################################################
#################################################################################

    def update_pileup_table(self, pileup_file):
        pileupline = pileup_file.readline()
        
        ### parse the pileup file line and get the position 
	### and the alleles at that location
        pileuplist = pileupline.split('\t')

        ## signifies that pileup line is either too short, or we've come to the end of the file
        if pileuplist.size != 6: 
            sys.stderr.write("WARNING: either pileup file doesn't have "
                         "anymore lines, or pileup file line has different number of fields."
			 " Potential file corruption. Aborting." )
            self.clear()
            return

        

        ### a neat method to simply rearrange the the elements in 
	### the array so that the first position becomes the last position ###
        np.roll(self.nppos,-1) 
        np.roll(self.npreads, -1)

        #### store necessary parts of pileup line. 
	#### namely the position and the number of ^ characters its reads have
        basepairs = np.array(list(pileuplist[4]))
        self.nppos[self.readlength-1] = pileuplist[1]
        self.npreads[self.readlength-1] = len(np.where(basepairs == "^")[0])
        self.currentreads = pileuplist[4]
        self.currentpos = pileuplist[1]
        
#################################################################################
#################################################################################
######### the fill function simply fills up the nppos and npreads with the ######### 
######### first X lines (where X is equivalent to the average read length). ######### 
######### It's used once at the beginning of the script, and never used again #######
#################################################################################
#################################################################################


    def fill(self, pileup_file, read_length):

	#### read the lines from a pileup file ####
        for x in range(0, read_length):
            pileupline = pileup_file.readline()
            pileuplist = pileupline.split('\t')

	#### a metric to show a corrupted pileup file ####
            if pileuplist.size != 6:
                sys.stderr.write("WARNING: either pileup file doesn't have "
                         "anymore lines, or pileup file line has different number of fields."
			 " Potential file corruption. Aborting." )
                self.clear()
                return

       ###### store all necessary information (potentially make this a submethod)
	###### let's do it. 

            basepairs = np.array(list(pileuplist[4]))
            self.nppos[x] = pileuplist[1]
            self.npreads[x] = len(np.where(basepairs == "^")[0])
            self.currentreads = pileuplist[4]
            self.currentpos = pileuplist[1] 
            

#################################################################################
#################################################################################
######### update_snp is self-explanatory. It reads another line from the SNP ######### 
######### file and updates all of our SNP information  #########  #########  ######### 
#################################################################################
#################################################################################


    def update_snp(self, snp_file):
        snpline = snp_file.readline()
        snplist = snpline.split('\t')

        if snplist.size < 3:
            sys.stderr.write("WARNING: either snp file doesn't have "
                         "anymore lines, or snp file line has different number of fields."
			 " Potential file corruption. Aborting." )
            self.clear()
            return 

	##### storing snp information ######
        self.snp_pos = int(snplist[0])
        self.snp_pos_minus_read_length = self.snp_pos - self.readlength
	
	### making sure haplotype information exists
        if len(snplist) > 3:
            self.hap = snplist[3]

	### making sure all the refs and alts have been stored properly
        self.ref = snplist[1]

	### this tries to find multiple alleles (for example if alt has C,G) 
        if len(self.alt1) > 1:
            alt = snplist[2].split(',')
            self.alt1 = alt[0]

            if len(alt) > 1:
                self.alt2 = alt[1]

            if len(alt) == 3:
                self.alt3 = alt[2]
        else:
            self.alt1 = snplist[2]


#### counting bases in our reads ######
    def count(self, bases):
        
        self.refcount = len(bases[np.where(bases==self.ref)])
        self.alt1count = len(bases[np.where(bases==self.alt1)])
        self.alt2count = len(bases[np.where(bases==self.alt2)])
        self.alt3count = len(bases[np.where(bases==self.alt3)])



#################################################################################
#################################################################################
######### create_output creates the actual starter position file for a ######### 
######### particular SNP. ######### ######### ######### ######### ######### ######### 
#################################################################################
################################################################################# 



    def create_output(self, chromosome, min_reads, output_directory):
        
        counter = self.readlength-1
        self.currentpos = self.snp_pos        
        
        
        if self.hap == '1|1' or self.hap == '0|0':
            return
        
        
        filename = chromosome + "." + str(self.snp_pos)
        filename = output_directory + "/" + filename
        outfile = open(filename, "w")
        outfile.write("SNP_POS\tREF\tAALT\tREF_COUNT\tALT_1_COUNT\tALT_2_COUNT\tALT_3_COUNT\tHAP\n")
        
        #### this class keeps track of the bases that are shrunk and come at each position ####
        the_bases = bases.BASES(self.currentreads)

        if len(the_bases.basepairs) < min_reads:
            print('Reads at ' + str(self.snp_pos) + ' do not mean minimum requirement\n')
            return
    
        if len(set(filter(str.isalpha, the_bases.basepairs)))==1:
            print('Homozygous at supposedly heterozygous position at ' + str(self.snp_pos) + "\n")
            return

        if len(self.ref) > 1:
            print('Invalid Reference Allele(s) at ' + str(self.snp_pos) + "\n")
            return
        if self.ref == '-':
            print('Invalid Reference Allele(s) at ' + str(self.snp_pos) + "\n")
            return
        if len(self.alt1) > 1 or len(self.alt2) > 1 or len(self.alt3) > 1:
            print('Invalid Alternate Allele(s) at ' + str(self.snp_pos) + "\n")
            return
  
	#### the algorithm for figuring out which bases belong to which haplotype ####
        for index in range(0, self.readlength):
            the_bases.shrink(counter, self.npreads)
            self.count(the_bases.currentbases)
            outfile.write(str(self.currentpos)+ "\t"+ self.ref + "\t" + self.alt1 + "\t"+ str(self.refcount) +"\t" +str(self.alt1count)+ "\t" +str(self.alt2count)+ "\t"+ str(self.alt3count) + "\t"+ self.hap)
            counter = counter - 1
            self.currentpos = self.currentpos - 1
 
        outfile.close()
