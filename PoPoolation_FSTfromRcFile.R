#!/usr/bin/env Rscript

#############################################################################################
###	Option, Arguments declarations and usage
#############################################################################################

library("optparse")
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Popoolation2 rc file ", metavar="character"),
  make_option(c("-l", "--list"), action="store_true", default=FALSE,
              help="file provided is a list of file "),
  make_option(c("-g", "--graphical"), action="store_true", default=FALSE,
              help="print fancy dendrogram !"),
  make_option(c("-i", "--legend"), type="character", default=NULL,
              help="use tsv file with first col as ordered name of population !"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output suffix", metavar="character")
  );
usage = "This R script will compute global FST on pairs of population from rc file of PoPoolation2 run.\nIt will write two files with a FST (Fu and Fw) matrix oz size [n, n] (n: number of populations)."
opt_parser = OptionParser(option_list=option_list, description=usage);
opt = parse_args(opt_parser);

if ( is.null(opt$file) & is.null(opt$out) ){
  print_help(opt_parser)
  stop("Two arguments must be supplied (input and output files names).", call.=FALSE)
}

#############################################################################################
###	Function to compute FST values by pop across SNP sites (From Popoolation2 rc output)
#############################################################################################

#################################################
### Main Function
#From a set of SNP between two population : Compute Fst by SNP, FST global (Weir and Cockerham, 1983 and Bhatia et al., 2013) see Jackson et al. (2014) for literal development of equations.
# According Jackson et al. (2014) :
# ---------------------------------
# PIs is the expected divergence between a pair of allele from two different populations
# PIb is the expected within poputation diversity
# Jackson (2014) used the following estimation of PIs and PIb from Charlesworth (1998) :
compute_PIs = function(Pa,Pb){
	return( Pa*(1-Pa) + Pb*(1-Pb) )
}
compute_PIb = function(Pa,Pb){
	return( Pa*(1-Pb) + Pb*(1-Pa) )
}
# with P the frequence of an allelic state common between population A and B.
# Pa is the frequence of this state in the population A
# Pb is the frequence of this state in the population B


# According Weir and Cockerham (1984) and Hudson et al. (1992), FST which is abbreviated as F can be written :
# this estimator is used for a single allele between two population
compute_F = function(PIs, PIb){
	return( (PIb-PIs)/PIb )
}

# In case of multiple allele (for desambiguation multiple locus with bi-allelic state) and according reccomendation of Weir and Cockerham (1984); Bhatia et al. (2013) :
# We can define to cross-loci FST, the first one Fu defined as the average of the FST and a second (Fw) take into account the cross-loci PIs and PIb to compute a FST statistics. With the words of Chen et al. (2016) in G3 :
# "we took the average of individual SNP FST values, or we calculated the ratio of the expectations of the denominator and numerator across all SNPs.

# With Fvec, a vector of the F for S sites :
compute_Fu = function(Fvec){
	S = length(Fvec)
	return( sum(Fvec)/S )
}
#
compute_Fw = function(PIs, PIb){
	S = length(PIs)
	stopifnot(S == length(PIb))
	return( sum(PIb - PIs) / sum(PIb) )
}

# According the result of Bhatia et al. (2013). The Fu mesures give equal weigth for every SNPs and Fw with higher expected levels of polymorphism

#################################################
### Formating Function
# Precedent functions need Pa and Pb (frequence of an allelic state common between population A and B)
# Because Major is less prone to error from sequencing error (and more frequently shared by population). This will be the used allele. But first we need to know how many population are stored in the file.
NumberOfPop = function(rcDf){
	DimFile=dim(rcDf)
	Npop=(DimFile[2]-9)/2
	return(Npop)
}

Determination_CommonAlleleFrequency = function(vectorInfosSite){
# 	Accept a vector containing (for one site) : MaaA, MaaB, MiaA, MiaB, MaafA, MaafB, MiafA, MiafB
# 	return a vector of lenth 2 with the correct composition c(maa,maa), c(mia, mia), c(mia,maa), c(maa, mia), c(nc,nc)
# 	Reformat the input vector :
	ReconstituedVector = unlist(strsplit(vectorInfosSite ,' '))
	#if MaaA == MaaB
	if(ReconstituedVector[1] == ReconstituedVector[2]){
		return( c(ReconstituedVector[5],ReconstituedVector[6]) )
	}else if(ReconstituedVector[1] == ReconstituedVector[4]){ #if MaaA == MiaB
		return( c(ReconstituedVector[5],ReconstituedVector[8]) ) 
	}else if(ReconstituedVector[3] == ReconstituedVector[2]){ #if MiaA == MaaB
		return( c(ReconstituedVector[7],ReconstituedVector[6]) ) 
	}else if(ReconstituedVector[3] == ReconstituedVector[4]){ #if MiaA == MiaB
		return( c(ReconstituedVector[7],ReconstituedVector[8]) )
	}else{
		return( c(NA, NA))
	}
}


# Extract informations needs fo a pair of population
PrepareAndCompute_FuFw = function(rcDf, A, B){
	npop=NumberOfPop(rcDf)
# 	A should be an integer from 1 to n (number of population)
# 	B should be an integer from 1 to n (number of population)
	stopifnot(A %in% seq(1,npop) & B %in% seq(1,npop))
	
# 	Look for Allelic state
	Maa = t(as.data.frame(strsplit(as.character(rcDf[,8]),''), row.names=NULL, col.names=NULL, fix.empty.names=F))
	Mia = t(as.data.frame(strsplit(as.character(rcDf[,9]),''), row.names=NULL, col.names=NULL, fix.empty.names=F))
	MaaA = as.character(Maa[,A])
	MaaB = as.character(Maa[,B])
	MiaA = as.character(Mia[,A])
	MiaB = as.character(Mia[,B])
	
#	store allele Count for minor and major (vector with two position: i. allele count, and ii. sites count )
	MaacA = t(as.data.frame(strsplit(as.character(rcDf[,c(9+A)]),'/'), row.names=NULL, col.names=NULL, fix.empty.names=F))
	MiacA = t(as.data.frame(strsplit(as.character(rcDf[,c(9+npop+A)]),'/'), row.names=NULL, col.names=NULL, fix.empty.names=F))
	MaacB = t(as.data.frame(strsplit(as.character(rcDf[,c(9+B)]),'/'), row.names=NULL, col.names=NULL, fix.empty.names=F))
	MiacB = t(as.data.frame(strsplit(as.character(rcDf[,c(9+npop+B)]),'/'), row.names=NULL, col.names=NULL, fix.empty.names=F))

#	store allele frequency for minor and major 
	MiafA = as.numeric(MiacA[,1])/as.numeric(MiacA[,2])
	MaafA = as.numeric(MaacA[,1])/as.numeric(MaacA[,2])
	MiafB = as.numeric(MiacB[,1])/as.numeric(MiacB[,2])
	MaafB = as.numeric(MaacB[,1])/as.numeric(MaacB[,2])

	
# 	Return the Common Allele frequency
	AlleleStateFrequency = paste(MaaA, MaaB, MiaA, MiaB, MaafA, MaafB, MiafA, MiafB)
	CommonAlleleFrequencyVector = t(as.data.frame(lapply(c(AlleleStateFrequency), Determination_CommonAlleleFrequency), row.names=NULL, col.names=NULL, fix.empty.names=F))
	CommonAlleleFrequencyReady = na.omit(CommonAlleleFrequencyVector)
	PiS = compute_PIs(as.numeric(CommonAlleleFrequencyReady[,1]), as.numeric(CommonAlleleFrequencyReady[,2]))
	PiB = compute_PIb(as.numeric(CommonAlleleFrequencyReady[,1]), as.numeric(CommonAlleleFrequencyReady[,2]))
	Fvector = compute_F(PiS, PiB)
	Fu = compute_Fu(na.omit(Fvector)) # Ffector return NaN when both allele are fixed, or identical
	Fw = compute_Fw(PiS, PiB)
	return(c(Fu, Fw, length(PiS)))
# given a df with MajorA MinorA MajorB MinorB CountMajA/CountA
}

main = function(fileIN){
	data = read.csv(fileIN, h=T, sep="\t" )
	Npop = NumberOfPop(data)
	Infos = matrix(nrow=Npop, ncol=Npop)
	FuI = matrix(nrow=Npop, ncol=Npop)
	FwI = matrix(nrow=Npop, ncol=Npop)
	for(i in seq(1, Npop)){
		for(j in seq(1, Npop)){
			cat("Working on", fileIN, "pair of pop:",i,"-",j,"\n" )
			if(i==j){
				FuI[i,j] = 0
				FwI[i,j] = 0
			}else if(is.na(FuI[i,j])==FALSE){ # to prevent computation of miror Fst (pair[i,j] == pair[j,i])
				next
			}else{
				result = PrepareAndCompute_FuFw(data, i, j)
				FuIpop = result[1]
				FwIpop = result[2]
				Nsites = result[3]
				FuI[i,j] = FuIpop
				FuI[j,i] = FuIpop
				FwI[i,j] = FwIpop
				FwI[j,i] = FwIpop
				Infos[i,j] = Nsites
				Infos[j,i] = Nsites
			}
		}
	}
	write.table(file=paste(fileIN, opt$out, "_fw",sep=""),FwI, sep="\t", quote=F, row.names=F, col.names=F)
	write.table(file=paste(fileIN, opt$out, "_fu",sep=""),FuI, sep="\t", quote=F, row.names=F, col.names=F)
	write.table(file=paste(fileIN, opt$out, "_Nsites",sep=""),Infos, sep="\t", quote=F, row.names=F, col.names=F)
	
	# Produce fancy graph if -g option
	if(opt$graphical){
		if(is.null(opt$legend)){
			PopName = seq(1, length(FwI[,1]), by=1 )
		} else{
			infos = read.csv(sep="\t", h=F, file=opt$legend)
			PopName = infos$V1
		}
		names(FwI) = PopName
		row.names(FwI) = PopName
		fw_r = cbind(FwI, PopName) 
		names(FuI) = PopName
		row.names(FuI) = PopName
		fu_r = cbind(FuI, PopName) 
		
		pdf(file=paste(fileIN, opt$out, "_fst.pdf",sep=""))
		par(mfrow=c(1,2))
		plot(hclust(as.dist(fu_r[,1:length(fu_r[,1])]), method="average"), xlab='UPGMA clustering of the 12 pops using Fst_u')
		plot(hclust(as.dist(fw_r[,1:length(fw_r[,1])]), method="average"), xlab='UPGMA clustering of the 12 pops using Fst_w')
		dev.off()
	}
	
	
}

#======================================================================================
#	Main
#======================================================================================
# logic of this statement is inverted, I can't store_true from optparse option...
if(opt$list){
	dataList = readLines(opt$file)
	for( i in seq(1, length(dataList))){
		main(dataList[i])
	}
	
} else {
	main(opt$file)
}
# testFile = "12pops_rc-test"
# main(testFile)
