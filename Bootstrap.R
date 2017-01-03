#!/srv/gs1/software/R/3.2.0/bin/Rscript

## The purpose of this script is to take a SPF output, and develop a bootstrapped ratio of haplotype 1 / haplotype 2 ###
library(data.table)
library(plyr)


#### First argument is the SPF, and the second argument is how many times you want the bootstrap function to sample ####
args = commandArgs(trailingOnly = TRUE)
spf = read.table(args[1], header = T)
bootstrap_interval = args[2]


### core bootstrap function: simple####
bootstrap <- function(spf) {

	ratios = c()
	
	### for a certain number of times ####
	for (count in 1:bootstrap_interval){
	
		adjustment_count = length(mribo$REF_COUNT)

		### randomly sample x number of positions
		sampleribo = sample.int(adjustment_count, size = adjustment_count, replace = TRUE)

		### sum up the counts
		matribocounts = spf$REF_COUNT[c(sampleribo)]
		patribocounts = spf$ALT_1_COUNT[c(sampleribo)]
	
		### ratio
		ratio = (sum(matribocounts) + adjustment_count) / (sum(matribocounts)+ sum(patribocounts) + 2*adjustment_count)

		### add to the vector of ratios
		ratios = c(ratios, ratio)

	}


	return ratios

	}
	
#### get the lower bound and upper bound of a 99% confidence intreval and print
vector_of_ratios = bootstrap(spf)
upperbound = quantile(vectorribo, 0.99, na.rm = TRUE)
lowerbound = quantile(vectorribo, 0.01, na.rm = TRUE)
medianbound = median(vectorribo)
CI = paste(lowerbound, upperbound, medianbound, sep = "\t")	
print(CI)
