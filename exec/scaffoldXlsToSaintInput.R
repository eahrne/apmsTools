#!/usr/bin/Rscript

# Converts a Scaffold Spectral Count export to Saint (interaction.txt, prey.txt, bait.txt)  or CRAPome input files   
# 
# Author: erikahrne
###############################################################################

kVersion <- "0.1"

# FUNCTIONS

getUserOptions <- function(version=version){
	
	# get command line options, if help option encountered print help and exit,
	# otherwise if options not found on command line then set defaults,
	cmdOpt <- parse_args(OptionParser( prog=paste("scaffoldXlsToSaintInput",kVersion), option_list=optionList))
	
	### CMD OPTIONS END						
	
	### SET USER OPTIONS
	userOptions <- list()
	
	### VERBOSE
	#VERBOSE: verbose
	userOptions$verbose <- cmdOpt$verbose
	### VERBOSE	END
	
# I/O
	#I/O: targetInputFile
	userOptions$targetInputFile <- cmdOpt$targetInputFile
	if( userOptions$targetInputFile == "" ){
			stop("Please specify a TARGET input file.")
	}else if(!file.exists(userOptions$targetInputFile)){
		stop("Input file ", userOptions$targetInputFile, " not found!")
	}
	
	#I/O: controlInputFile
	userOptions$controlInputFile <- cmdOpt$controlInputFile
	if( userOptions$controlInputFile == "" ){
		stop("Please specify a CONTROL input file.")
	}else if(!file.exists(userOptions$controlInputFile)){
		stop("Input file ", userOptions$controlInputFile, " not found!")
	}
	
	#I/O: outputDirPath
	userOptions$outputDirPath <- cmdOpt$outputDir
	if(!file.exists(userOptions$outputDirPath) & userOptions$outputDirPath != "" ){
		stop("No such directory: ",userOptions$outputDirPath)		
	}else if((substr(userOptions$outputDirPath,nchar(userOptions$outputDirPath),nchar(userOptions$outputDirPath)) != "/")
			||
			!(userOptions$outputDirPath== "")
			){
		userOptions$outputDirPath <- paste(userOptions$outputDirPath,"/", sep="")
	}
	
	
	
# I/O END

# .FASTA FILE	
	userOptions$fastaFile <- cmdOpt$fastaFile 
	if(!is.na(userOptions$fastaFile) && !file.exists(userOptions$fastaFile)){
		stop("Fasta file ", userOptions$fastaFile, " not found!")
	}
# .FASTA FILE END		

# CRAPome EXPORT	
	userOptions$crapomeExport <- cmdOpt$crapomeExport	
# CRAPome EXPORT END	
	
	if(userOptions$verbose){
		cat("### SPECIFIED OPTIONS ###\n\n")
		for(opt in names(userOptions)){
			cat(opt,":",userOptions[[opt]],"\n")
		}
		cat("\n### SPECIFIED OPTIONS ###\n\n")
	}
	
	return(userOptions)	
}

# FUNCTIONS END

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(apmsTools))
### LOAD EXT LIBRARIES END
### CMD OPTIONS
optionList <- list(
		### I/O
		make_option(c("-t", "--targetInputFile"), type="character", default="",
				help="I/O: Scaffold .xls target file (REQUIRED)"
		),
		make_option(c("-c", "--controlInputFile"), type="character", default="",
				help="I/O: Scaffold .xls control file (REQUIRED)"
		),	
		make_option(c("-o", "--outputDir"), type="character", default=".",
				help="I/O:  Results Output Directory [default ./]"
		),
		### I/O
		
		### .FASTA FILE
		make_option(c("-f","--fastaFile"), type="character", default=NA,
				help="Path to protein sequnece .fasta file [default None]
		Used to update protein lengths.",
				metavar=".fasta file path"),
		### .FASTA FILE END
		
		make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
			help="Print progress information [default %default]"
		),
		make_option(c("-r", "--crapomeExport"), action="store_true", default=FALSE,
				help="Create CRAPome input file [default %default]"
		)
	
)

### I/O END

#get User Options
userOptions <- getUserOptions(kVersion)

aExp <- apmsExp(readScaffoldSpecCountFile(userOptions$targetInputFile)
		,readScaffoldSpecCountFile(userOptions$controlInputFile, isControl=T))

# extract and set protein lengths from user specified .fasta
if(!is.na(userOptions$fastaFile)){
	if(userOptions$verbose){
		cat("\nUpdating protein lengths \n")
	}
	aExp <- setProteinLength(aExp,userOptions$fastaFile)
}

if(userOptions$crapomeExport){
	writeCRAPomeFile(aExp,paste(userOptions$outputDirPath,"CRAPomeInput.txt",sep="/"),verbose=T)
}else{
	writeSaintInteractionFile(aExp,paste(userOptions$outputDirPath,"interaction.txt",sep="/"),verbose=T)
	writeSaintPreyFile(aExp,paste(userOptions$outputDirPath,"prey.txt",sep="/"),verbose=T)
	writeSaintBaitFile(aExp,paste(userOptions$outputDirPath,"bait.txt",sep="/"),verbose=T)
}



quit(status=1)
