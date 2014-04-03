# TODO: Add comment
# 
# Author: erikahrne
###############################################################################
library(RUnit)

### INIT
kBaseDir <- "/Users/erikahrne/dev/R/workspace/apmsTools/"
source(paste(kBaseDir,"/R/apmsExp.R",sep=""))

scaffoldTargetFile <- paste(kBaseDir,"data/PP2A-p53.xls",sep="/")
scaffoldControlFile <- paste(kBaseDir,"data/GFP5.xls",sep="/")
scaffoldControlNonRedundantFile <- paste(kBaseDir,"data/GFP_non-redundant_protein_cluster_analysis.xls",sep="/")
fastaFile <-  paste(kBaseDir,"data/c_human_TG.decoy_Dehio_Ficproteins_04112013.fasta",sep="")

aExp <- apmsExp(readScaffoldSpecCountFile(scaffoldTargetFile)
		,readScaffoldSpecCountFile(scaffoldControlNonRedundantFile, isControl=T))

### INIT END


### TESTS

testReadScaffoldSpecCountFile <-function(){
	
	cat("--- testReadScaffoldSpecCountFile ---" ,checkEqualsNumeric(nrow(readScaffoldSpecCountFile(scaffoldTargetFile)),565),"\n")
	cat("--- testReadScaffoldSpecCountFile ---" ,checkEqualsNumeric(nrow(readScaffoldSpecCountFile(scaffoldControlFile)),450),"\n")
	cat("--- testReadScaffoldSpecCountFile ---" ,checkEqualsNumeric(nrow(readScaffoldSpecCountFile(scaffoldControlNonRedundantFile)),655),"\n")
	
} 

testApmsExp <-function(){
	cat("--- testApmsExp ---", checkEqualsNumeric(
								ncol(apmsExp(readScaffoldSpecCountFile(scaffoldTargetFile)
										   ,readScaffoldSpecCountFile(scaffoldControlNonRedundantFile))$specCountDf)
							   ,11),"\n")
}

testWriteSaintBaitFile <- function(){
	
	cat("--- testWriteSaintBaitFile --- ")
	writeSaintBaitFile(aExp,paste(tempdir(),"bait.txt",sep="/"),verbose=T)
	
}

testWriteSaintInteractionFile <- function(){
	
	cat("--- testWriteSaintInteractionFile --- ")
	writeSaintInteractionFile(aExp,paste(tempdir(),"interaction.txt",sep="/"),verbose=T)
	
}


testWriteSaintPreyFile <- function(){
		
	cat("--- writeSaintPreyFile --- ")
	writeSaintPreyFile(aExp,paste(tempdir(),"prey.txt",sep="/"),verbose=T)
	
}

testSetProteinLength <- function(){
	cat("--- testSetProteinLength ---" ,checkEqualsNumeric(setProteinLength(aExp,fastaFile)$specCountDf$proteinLength[565],222),"\n")
}


### TESTS END

### RUN TESTS

testReadScaffoldSpecCountFile()
testApmsExp()
testWriteSaintBaitFile()
testWriteSaintInteractionFile()
testWriteSaintPreyFile()
testSetProteinLength()

### RUN TESTS END




	