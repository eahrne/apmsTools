\name{readScaffoldSpecCountFile}
\alias{readScaffoldSpecCountFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parse scaffold spectral count export (.xls).
}
\description{
Parses scaffold .xls and creates data.frame object.

Annotate replicates in scaffold file by adding replicate number to bait gene name.

IP-name: PPP2R1A-1 (as given in scaffold)

bait-name: PPP2R1A

}
\usage{
readScaffoldSpecCountFile(file, sep = "\t", isControl = F, isUniprot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
Path to Scaffold .xls File.
}
  \item{sep}{
the field separator character.
}
  \item{isControl}{
Does the file contain GFP control baits only 
}
  \item{isUniprot}{
Uniprot Protein Accession Numbers TRUE/FALSE
}
}
\details{
Not yet available
}
\value{
	Returns a data frame with the following structure:
	 
	rownames: protein accession numbers
	
    df[,1:3]: "protein description","gene name","'approx.' protein length" 
    
    df[,4:]: bait-prey spectral count.

}
\references{
Not yet available
}
\author{
Erik Ahrne

erik.ahrne@unibas.ch
}
\note{
Not yet available
}

\seealso{
Not yet available
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
