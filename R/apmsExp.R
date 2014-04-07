# "apmsExp", S3 class including functions for parsing, storing and exporting apms spectral count data
# 
# Author: erikahrne
###############################################################################
apmsExp <- function( targetSpecCountDf,controlSpecCountDf){
	
	# Creates an apms experiment data container
    # class: "apmsExp"
	#
    # Args:
    #   targetSpecCount: target (bait-prey) spec. count data.frame(). 
	#		rownames: prey proteinACs,
    #		#df[,1:3]: "protein description","gene name","protein length" 
    #		df[,4:n] bait-prey spec counts	
    #   controlSpecCount: GFP control spec. count data.frame() 
    #		#df[,1:3]: "protein description","gene name","protein length" 
    #		df[,4:n] GFP-prey spec counts	
    #
    # Returns: apmsExp object
	#   desc

	
	### merge dataframes
	uniqueProteinACs <- unique(rownames(targetSpecCountDf),rownames(controlSpecCountDf))
	merge <- cbind(targetSpecCountDf[uniqueProteinACs,4:ncol(targetSpecCountDf)]
			,controlSpecCountDf[uniqueProteinACs,4:ncol(controlSpecCountDf)]
	)
	### replace NA's with 0
	merge[is.na(merge)] <- 0
	merge <- cbind(rbind(targetSpecCountDf[,1:3],controlSpecCountDf[,1:3])[uniqueProteinACs,],merge)
	
	
	out <- list()
	
	out$specCountDf <- merge
	
	class(out) <- "apmsExp"
	
	return(out)
	
}

readScaffoldSpecCountFile <- function(file,sep="\t", isControl=F, isUniprot=T){
	
	# Reads scaffold spectral count export and creates 
	#
    # Args:
    #   file: path to scaffold .xls spectrum count export
    #   sep: deliminator
	#	isControl: is the data from GFP control experiments
	#	isUniprot: Were protein identified searching a uniprot ONLY database?
	
    #
    # Returns:
    #   A data frame 
	# rownames: protein accession numbers
    # df[,1:3]: "protein description","gene name","'approx.' protein length" 
	#
	# Commments: 
	# 	Annotate replicates in scaffold file by adding replicate number to bait gene name.
	# 	IP-name PPP2R1A-1 (as given in scaffold)
	# 	bait-name PPP2R1A
	
	# read file
	specCountDf <- read.csv(file,skip=2,sep=sep, head=T)
	# remove END OF FILE
	specCountDf <-  specCountDf[1:(nrow(specCountDf)-1),]
	
	
	### manage protein cluster analysis export. 1) Discard leadinf protein, 2) Remove protein description indentation
	# Example:
	#Cluster of Heat shock protein HSP 90-alpha OS=Homo sapiens GN=HSP90AA1 PE=1 SV=5 (sp|P07900|HS90A_HUMAN)                                          
	#     Heat shock protein HSP 90-alpha OS=Homo sapiens GN=HSP90AA1 PE=1 SV=5  
	#	  Putative heat shock protein HSP 90-beta 2 OS=Homo sapiens GN=HSP90AB2P PE=1 SV=2 
	
	# 1) Discard leadinf protein
	proteinClusterSel <- 	regexpr("^Cluster of ",specCountDf[,3]) == -1
	specCountDf <- 	specCountDf[proteinClusterSel,]
	
	# 2) Remove protein description indentation
	proteinDesc <- gsub("^ *","",specCountDf[,3])
	
	# discard grouping info sp|P60709|ACTB_HUMAN (+1) ->  sp|P60709|ACTB_HUMAN
	proteinAC <- gsub(" .*","",specCountDf[,4])

	# approx. assuming avg. a.a mass 100Da  .Call setProteinLength(file) to get protein lenghts from .fasta file
	proteinLength <- as.numeric(gsub(" kDa","",specCountDf[,5]))*10
	
	specCounts <- specCountDf[,10:(ncol(specCountDf)-1)]
	# label baits as GFP1 ... GFPn
	if(isControl){
		names(specCounts) <- paste("GFP",1:ncol(specCounts),sep="_")
	}
	
	### parse gene name, 
	geneName <- rep(NA,nrow(specCountDf))
	# E.g. Keratin, type II cytoskeletal 1 OS=Homo sapiens GN=KRT1 PE=1 SV=6
	# geneName -> KRT1
	if(isUniprot){ 
		geneName <- gsub(" PE\\=.*","",gsub(".*GN\\=","",proteinDesc))
		# check if failed to parse gene name
		nasSel <- regexpr("GN\\=",proteinDesc) == -1
		geneName[nasSel] <- NA
		
		if(sum(nasSel) > 0 ){
			warning("Found ", sum(nasSel)," non Uniprot entries.", call. = FALSE)
		}
	}
	
	out <- data.frame(proteinDesc,geneName,proteinLength,specCounts)
	rownames(out) <- proteinAC
	names(specCountDf)
	
	return(out)
	
}

writeSaintBaitFile <- function(apmsExp,file,verbose=F){
	
	# Write bait SAINT input file  
	# Example:
	#	PPP2R1A	PPP2R1A	T
	#	PPP2CB	PPP2CB	T
	#	PPME1	PPME1	T
	#	TP53	TP53	T
	#	GFP1	GFP1	C
	#	GFP2	GFP2	C
	#	GFP3	GFP3	C
	#	GFP4	GFP4	C
	#	
	#    # class "apmsExp"
    #
    # Args:
    #   apmsExp: apmsExp object
    #   file: path to bait output file
    #
	
	#@ How to label replicates?
	ipName <- names(apmsExp$specCountDf)[4:ncol(apmsExp$specCountDf)]
	baitName <- gsub("\\-[0-9]{1,}$","",ipName) # strip-off replicate number
	isControlStr <- c("T","C")[(regexpr("^GFP",ipName) > -1)+1]
	
	baitOut <- data.frame(ipName,baitName,isControlStr)
	write.table(baitOut,file=file,row.names=F,col.names=F,sep="\t", quote=F)
	
	if(verbose){
		cat("Created File ",file,"\n")
	}
	
}

writeSaintInteractionFile <- function(apmsExp,file,verbose=F){
	
	# Write interaction SAINT input file  
	# Example:
	#	PPP2R1A	PPP2R1A	sp|P30153|2AAA_HUMAN	64
	#	PPP2R1A	PPP2R1A	sp|Q9BUJ2|HNRL1_HUMAN	62
	#	PPP2R1A	PPP2R1A	sp|P08107|HSP71_HUMAN	29
	#	
	#    # class "apmsExp"
    #
    # Args:
    #   apmsExp: apmsExp object
    #   file: path to interaction output file
    #
	
	nbPrey <- nrow(apmsExp$specCountDf)
	ipName <- names(apmsExp$specCountDf)[4:ncol(apmsExp$specCountDf)]
	baitName <- gsub("\\-[0-9]{1,}$","",ipName) # strip-off replicate number
	
	interactionIpName <- c()
	interactionBaitName <- c()
	interactionPreyName <- c()
	interactionSpecCount <- c()
	
	for(i in 1:length(ipName)){
		
		sel <- apmsExp$specCountDf[,ipName[i]] > 0
		selLength <- sum(sel)
		
		interactionIpName <- c(interactionIpName,rep(ipName[i],selLength))
		interactionBaitName <- c(interactionBaitName,rep(baitName[i],selLength))
		interactionPreyName <- c(interactionPreyName,rownames(apmsExp$specCountDf)[sel])
		interactionSpecCount <- c(interactionSpecCount,apmsExp$specCountDf[sel,ipName[i]])
		
	}
	
	interactionOut <- data.frame(interactionIpName
			,interactionBaitName
			,interactionPreyName
			,interactionSpecCount	
	)
	
	write.table(file=file,interactionOut,row.names=F,col.names=F,sep="\t",  quote=F)
	
	if(verbose){
		cat("Created File ",file,"\n")
	}
	
}

writeSaintPreyFile <- function(apmsExp,file,verbose=F){
	
	# Write prey SAINT input file  
	# Example:
	# 	sp|P30153|2AAA_HUMAN	650	PPP2R1A
	#	sp|Q9BUJ2|HNRL1_HUMAN	960	HNRNPUL1
	#	sp|P08107|HSP71_HUMAN	700	HSPA1A
	#	
	#   
    # class "apmsExp"
    #
    # Args:
    #   apmsExp: apmsExp object
    #   file: path to prey output file
    #
	
	preyAC <- rownames(apmsExp$specCountDf)
	preySeqLength <- apmsExp$specCountDf$proteinLength
	
	preyGeneName <- apmsExp$specCountDf$geneName
	preyOut <- data.frame(preyAC,preySeqLength,preyGeneName)
	
	write.table(file=file,preyOut,row.names=F,col.names=F,sep="\t", quote=F)
	
	if(verbose){
		cat("Created File ",file,"\n")
	}
}

print.apmsExp <- function(apmsExp){
	print(unclass(apmsExp))
}

plot.apmsExp <- function(apmsExp,...){
	pairs(apmsExp$specCountDf[,4:ncol(apmsExp$specCountDf)],...)
}

summary.apmsExp <- function(apmsExp){
	
	specCountOnlyDf <- apmsExp$specCountDf[,4:ncol(apmsExp$specCountDf)]
	isControl <- regexpr("^GFP\\_",names(specCountOnlyDf)) > -1
	
	ret <-list()
	
	ret$targetBaits <- names(specCountOnlyDf)[!isControl]
	ret$controlBaits <- names(specCountOnlyDf)[isControl]
	ret$nbPrey <- nrow(apmsExp$specCountDf)
	ret$totSpecCountPerBait <- apply(specCountOnlyDf,2,sum)

	class(ret) <- c("summary.apmsExp") 

	return(ret)
}

print.summary.apmsExp <-function(apmsExpSummary){
	cat("\n")
	cat("Target Baits: ",  apmsExpSummary$targetBaits,"\n")
	cat("Control Baits: ",  apmsExpSummary$controlBaits,"\n")
	cat("Nb. Prey: ",  apmsExpSummary$nbPrey,"\n")
	cat("Total Spectral Counts per Bait:\n")
	print(apmsExpSummary$totSpecCountPerBait)
	cat("\n")
}


setProteinLength <- function(apmsExp,fastaFile){
	
	# Updates protein lengths (extracted from .fasta) in apmsExp object. 
	# class: "apmsExp"
    # Args:
    #   apmsExp: apmsExp object
    #   fastaFile: full path to .fasta file
	#
    # Returns:
    #   An updated apmsExp object
	
	library(seqinr) # package for parsing .fasta file
	proteinDB <- read.fasta(fastaFile,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
	
	preyAC <- rownames(apmsExp$specCountDf)
	preySeqLength <- lapply(proteinDB[preyAC],nchar)
	
	### Warn but don't creash if incorrect .fasta was provided
	isNA <- is.na(names(preySeqLength))
	names(preySeqLength)[isNA] <- preyAC[isNA]
	preySeqLength[isNA] <- round(mean(unlist(preySeqLength)))
	preySeqLength <- as.vector(unlist(preySeqLength)) 

	if(sum(isNA)){
		warning("Missing protein lengths, Check .fasta file version and AC formatting.", call. = FALSE)
		print(preyAC[isNA])
	}
	
	apmsExp$specCountDf$proteinLength <- preySeqLength 
	detach(package:seqinr)
	
	return(apmsExp)
}

writeCRAPomeFile <- function(apmsExp,file,verbose=F){
	
	# Write CRAPome input file  
	# Example:
	# PROTID	GENEID	PROTLEN	C1_SPC	C2_SPC	HDAC5_1_SPC	HDAC5_2_SPC
	# PROTID	GENEID	PROTLEN	C	C	HDAC5	HDAC5
	# Q9Y6M1	IGF2BP2	599	9	2	2	0
	# Q9Y697	NFS1	457	0	0	14	8
	# Q9Y618	NCOR2	2525	0	0	19	16
	#	
	#   
    # class "apmsExp"
    #
    # Args:
    #   apmsExp: apmsExp object
    #   file: path to CRAPome output file
	#
	
	# create headers
	h1 <- c("PROTID","GENEID","PROTLEN",paste(names(apmsExp$specCountDf)[4:ncol(apmsExp$specCountDf)],"SPC",sep="_"))
	h2 <- gsub("^GFP_.*","C",c(h1[1:3],names(apmsExp$specCountDf)[4:ncol(apmsExp$specCountDf)]))
	
	# create output df
	out <- cbind(rownames(apmsExp$specCountDf),apmsExp$specCountDf[,-1])
	names(out) <- h2
 	
	# write h1 to file
	cat(h1,file=file,sep="\t")
	cat("\n",file=file,append=T)
	# write h2 to file
	cat(h2,file=file,sep="\t",append=T)
	cat("\n",file=file,append=T)
	# append data.frame to file
	write.table(out,file=file,append=T,sep="\t",row.names=F,col.names=F)
	
	if(verbose){
		cat("Created File ",file,"\n")
	}
}
