\name{setProteinLength}
\alias{setProteinLength}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update protein lengths 
}
\description{
Updates protein lengths (extracted from .fasta) in apmsExp object, 
}
\usage{
setProteinLength(apmsExp, fastaFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{apmsExp}{
An object of class apmsExp.
}
  \item{fastaFile}{
Full path to .fasta file
}
}
\details{
Not yet available
}
\value{
An updated apmsExp object

\item{specCountDf }{data.frame, df  
	
			df[,1:3]: "protein description","gene name","protein length" 
			
			df[,4:n] bait-prey spec counts 
		 }

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
