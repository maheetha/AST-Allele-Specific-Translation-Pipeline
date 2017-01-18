import sys
import numpy as np

######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### 
######### the BASES class keeps track of the bases that are #########
######### at a particular heterozygous position and shrinks #########
######### them accordingly based on how many "^" starter ############
######### characters exist in the previous positions ################
#####################################################################
#####################################################################


class BASES:

    def __init__(self, bases):
        self.initialize(bases)

    ### takes out all unnecessary characters from the read string ####
    ### turns the string into an array of characters #####
    def initialize(self, bases):
        self.array = np.array(list(bases))
        ind1 = np.where(self.array == '^')[0]

        if ind1.size != 0:
            ind2 = ind1 + 1
        else:
            ind2 = []

        self.ind = np.append(ind1, ind2)
        del ind1
        del ind2
        np.delete(self.array, self.ind)

     #### store array of characters ####
        self.basepairs = ''.join([i for i in self.array if i.isalpha()])
        self.basepairs = np.array(list(self.basepairs))
        self.currentbases = []


    #### shrinks the basepairs and adds to currentbases depending on the ####
    #### number of "^" at the index of npreads #####

    def shrink(self, index, npreads):
        
        num = int(npreads[index])
       
        if num == 0:
            string = ''
            self.currentbases = np.array(list(string))
            return

        self.currentbases = self.basepairs[-num:]
        self.basepairs = self.basepairs[:-num]
       
