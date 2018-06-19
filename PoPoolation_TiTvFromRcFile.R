#!/usr/bin/env Rscript

library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Popoolation2 rc file ", metavar="character"),
  make_option(c("-l", "--list"), action="store_true", default=FALSE,
              help="file provided is a list of file "),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output suffix", metavar="character")
  );
usage = "This R script will compute Transition/transversion from rc file of PoPoolation2 run.\nIt will return a line by file with the name of the file, (for n populations) the n col. folowing the filename are the number of fixed sites, then n col. with the number of ti and then tv."
opt_parser = OptionParser(option_list=option_list, description=usage);
opt = parse_args(opt_parser);

if ( is.null(opt$file) & is.null(opt$out) ){
  print_help(opt_parser)
  stop("Two arguments must be supplied (input and output files names).", call.=FALSE)
}
#======================================================================================
#	Functions
#======================================================================================
MaaOverMia_PopVector = function(VecMaaOverVecMia){
	MaaAndMia = unlist(strsplit(VecMaaOverVecMia, ","))
	# return a list of allel state (maj/min)
	return(paste(unlist(strsplit(as.character(MaaAndMia[1]), '')), unlist(strsplit(as.character(MaaAndMia[2]), '')), sep='/') )
}

tiOrtv = function(MaaOverMia){
	Transition = c("A/G", "G/A", "C/T", "T/C")
	Transversion = c("A/C", "A/T", "T/G", "T/A", "C/A", "C/G", "G/C", "G/T")
	majmin = MaaOverMia
	state = c()
	if(majmin %in% Transition){state = "ti"
	}else if(majmin %in% Transversion){state = "tv"
	}else{ state = "fixed"}
	#return the state of the sustitution (ti,tv, or unknown)
	return(state)
}

CountTiTv = function(vectorAlleleState){
	# Need a vector of allele states (maj/min) FUN: MaaOverMia_PopVector()
	TiTv_vector = unlist(lapply(vectorAlleleState, tiOrtv))
	Infos = summary(as.factor(TiTv_vector))
	Nsites = length(TiTv_vector)
	
	# return a vector with the number of transition, the number of transversion and the number of position examined
	return(c(Infos,Nsites))
}

ratioTiTv = function(CountTiTvVector){
	#Needs: a vector with count of Ti, Tv, fixed, Nsites
	TiTvRatio=CountTiTvVector[2]/CountTiTvVector[3]
	Pfixed=CountTiTvVector[1]/CountTiTvVector[3]
	#return: TiTvRatio and Nfixed
	return(TiTvRatio,Nfixed)
	
}

main = function(fileIN){
	cat("Working on ", fileIN, "\n" )
	data = read.csv(fileIN, h=T, sep="\t" )
	DimFile=dim(data)
	Npop=(DimFile[2]-9)/2
	Npos=DimFile[1]
	# I'm only interested by Maa:Col8 ; Mia:Col9... I will rebuild a matrix with maa/mia Mat(Npop, Npos) :
	vectorOfAllele = paste(data[,8], data[,9], sep=",")
	MatAlleleState=as.matrix(t(as.data.frame(lapply(vectorOfAllele,MaaOverMia_PopVector), row.names=NULL, col.names=NULL, fix.empty.names=F)))
	#Count Ti, Tv and fixed sites
	Bigresult = apply(MatAlleleState, 2, CountTiTv)
	rownames(Bigresult)=c("Fixed", "Ti", "Tv", "Total")
	#return fixed, ti, tv
	# print result
	write.table(file=paste(fileIN, opt$out, sep=""),Bigresult, sep="\t", quote=F, row.names=T, col.names=F)

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


