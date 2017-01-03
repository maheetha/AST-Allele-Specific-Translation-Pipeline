#!/srv/gs1/software/R/3.2.0/bin/Rscript

args = commandArgs(trailingOnly = TRUE)
table = read.table(args[1], header = T
gene_list = table$V1
gene_position = table$V2

riboconfidenceintervals = c()
rnaconfidenceintervals = c()
significance = c()
pvalues = c()

names = seq(1, 100, 1)
pvalue_het = c()

for (x in 1:length(table$V1)){
	maternalribo = paste(gene_list[x], gene_position[x], "maternalRIBO", sep = ".")
	paternalribo = paste(gene_list[x], gene_position[x], "paternalRIBO", sep = ".")
	maternalrna = paste(gene_list[x], gene_position[x], "maternalRNA", sep = ".")
	paternalrna = paste(gene_list[x], gene_position[x], "paternalRNA", sep = ".")

	maternalribotable = read.table(maternalribo)
	paternalribotable = read.table(paternalribo)
	maternalrnatable = read.table(maternalrna)
	paternalrnatable = read.table(paternalrna)

	vectorribo = c()
	vectorrna = c()
	vectorOfPvalues = c()

	for (o in 1:1000){
	
	sampleribo = sample.int(length(maternalribotable[,3]), size = length(maternalribotable[,3]), replace = TRUE)
	samplerna = sample.int(length(maternalrnatable[,3]), size = length(maternalrnatable[,3]), replace = TRUE)
	
	matribocounts = maternalribotable$V3[c(sampleribo)]
	patribocounts = paternalribotable$V3[c(sampleribo)]
	matrnacounts = maternalrnatable$V3[c(samplerna)]
	patrnacounts = paternalrnatable$V3[c(samplerna)]

	ratioribo = (sum(matribocounts) + length(maternalribotable[,3])) / (sum(matribocounts)+sum(patribocounts) + length(maternalribotable[,3]) + length(maternalribotable[,3]))
	ratiorna = (sum(matrnacounts) + length(maternalrnatable[,3])) / (sum(matrnacounts)+sum(patrnacounts) + (2*length(maternalrnatable[,3])))
	vectorribo = c(vectorribo, ratioribo)
	vectorrna = c(vectorrna, ratiorna)

	}
	
	meanribo = mean(vectorribo, na.rm = TRUE)
	meanrna = mean(vectorrna, na.rm = TRUE)


	ribolength = 30
	rnalength = 75
	upperboundribo = quantile(vectorribo, 0.99, na.rm = TRUE)
	lowerboundribo = quantile(vectorribo, 0.01, na.rm = TRUE)
	RIB = median(vectorribo)
	upperboundrna = quantile(vectorrna, 0.99, na.rm = TRUE)
	lowerboundrna = quantile(vectorrna, 0.01, na.rm = TRUE)
	RNA = median(vectorrna)
	print(upperboundribo)
	print(lowerboundribo)
	print(upperboundrna)
	print(lowerboundrna)

	riboCI = paste(lowerboundribo, upperboundribo, RIB, sep = "\t")
	rnaCI = paste(lowerboundrna, upperboundrna, RNA, sep = "\t")
	riboconfidenceintervals = c(riboconfidenceintervals, riboCI)
	rnaconfidenceintervals = c(rnaconfidenceintervals, rnaCI)

	if (	meanribo > meanrna){
		if (lowerboundribo > upperboundrna){
			significance = c(significance, "yes")		
		} else {
			significance = c(significance, "no")		
		}
	} else if (meanrna > meanribo){
		if (lowerboundrna > upperboundribo){
			significance = c(significance, "yes")		
		} else {
			significance = c(significance, "no")		
		}
	} else if (meanrna == meanribo){
		significance = c(significance, "no")
	}

}

finaltable = data.frame(gene_list, gene_position, riboconfidenceintervals, rnaconfidenceintervals, significance)

tablename = paste(args[2], "FINAL_1000", sep = ".")
write.table(finaltable, tablename, row.names = F, quote = F)
