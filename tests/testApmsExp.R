# TODO: Add comment
# 
# Author: erikahrne
###############################################################################


### INIT
if(file.exists("/Users/erikahrne/dev/R/workspace/apmsTools/inst/testData/")){
	kBaseDir <- "/Users/erikahrne/dev/R/workspace/apmsTools/inst/"
	source("/Users/erikahrne/dev/R/workspace/apmsTools/R/apmsExp.R")
}else{
	library(apmsTools)
	kBaseDir <- path.package("apmsTools")
}

scaffoldTargetFile <- paste(kBaseDir,"/testData/PP2A-p53.xls",sep="/")
scaffoldControlFile <- paste(kBaseDir,"/testData/GFP5.xls",sep="/")
scaffoldControlNonRedundantFile <- paste(kBaseDir,"/testData/GFP_non-redundant_protein_cluster_analysis.xls",sep="/")
fastaFile <-  paste(kBaseDir,"/testData/c_human_TG.decoy_Dehio_Ficproteins_04112013.fasta",sep="")

aExp <- apmsExp(readScaffoldSpecCountFile(scaffoldTargetFile)
		,readScaffoldSpecCountFile(scaffoldControlNonRedundantFile, isControl=T))

### INIT END


### TESTS

testReadScaffoldSpecCountFile <-function(){
	
	cat("--- testReadScaffoldSpecCountFile: --- \n")
	stopifnot(all.equal(nrow(readScaffoldSpecCountFile(scaffoldTargetFile)),565))
	stopifnot(all.equal(nrow(readScaffoldSpecCountFile(scaffoldControlFile)),450))
	stopifnot(all.equal(nrow(readScaffoldSpecCountFile(scaffoldControlNonRedundantFile)),655))
	cat("--- testReadScaffoldSpecCountFile: PASS ALL TEST --- \n")
} 

testApmsExp <-function(){
	
	cat("--- testApmsExp: --- \n")
	stopifnot(all.equal(ncol(apmsExp(readScaffoldSpecCountFile(scaffoldTargetFile)
					       ,readScaffoldSpecCountFile(scaffoldControlNonRedundantFile))$specCountDf)
						   ,11))	
	cat("--- testApmsExp: PASS ALL TEST --- \n")
}

testWriteSaintBaitFile <- function(){
	
	cat("--- testWriteSaintBaitFile --- ")
	writeSaintBaitFile(aExp,paste(tempdir(),"bait.txt",sep="/"),verbose=T)
	cat("--- testWriteSaintBaitFile: PASS ALL TEST --- \n")
	
}

testWriteSaintInteractionFile <- function(){
	
	cat("--- testWriteSaintInteractionFile: --- ")
	writeSaintInteractionFile(aExp,paste(tempdir(),"interaction.txt",sep="/"),verbose=T)
	cat("--- testWriteSaintInteractionFile: PASS ALL TEST --- \n")
}

testWriteSaintPreyFile <- function(){
		
	cat("--- writeSaintPreyFile: --- ")
	writeSaintPreyFile(aExp,paste(tempdir(),"prey.txt",sep="/"),verbose=T)
	cat("--- testWriteSaintPreyFile: PASS ALL TEST --- \n")
	
}

testSetProteinLength <- function(){
	cat("--- testSetProteinLength: ---" ,stopifnot(all.equal(setProteinLength(aExp,fastaFile)$specCountDf$proteinLength[565],222)),"\n")
	cat("--- testSetProteinLength: PASS ALL TEST --- \n")
}

testWriteCRAPomeFile <-function(){
	
	cat("--- testWriteCRAPomeFile: --- ")
	writeCRAPomeFile(aExp,paste(tempdir(),"CRAPome_input.tsv",sep="/"),verbose=T)
	cat("--- testWriteCRAPomeFile: PASS ALL TEST --- \n")
	
}

### TESTS END

### RUN TESTS

testReadScaffoldSpecCountFile()
testApmsExp()
testWriteSaintBaitFile()
testWriteSaintInteractionFile()
testWriteSaintPreyFile()
testSetProteinLength()
testWriteCRAPomeFile()

### RUN TESTS END







	