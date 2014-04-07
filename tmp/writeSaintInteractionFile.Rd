\name{writeSaintInteractionFile}
\alias{writeSaintInteractionFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeSaintInteractionFile(apmsExp, file, verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{apmsExp}{
%%     ~~Describe \code{apmsExp} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (apmsExp, file, verbose = F) 
{
    nbPrey <- nrow(apmsExp$specCountDf)
    ipName <- names(apmsExp$specCountDf)[4:ncol(apmsExp$specCountDf)]
    baitName <- gsub("\\-[0-9]{1,}$", "", ipName)
    interactionIpName <- c()
    interactionBaitName <- c()
    interactionPreyName <- c()
    interactionSpecCount <- c()
    for (i in 1:length(ipName)) {
        sel <- apmsExp$specCountDf[, ipName[i]] > 0
        selLength <- sum(sel)
        interactionIpName <- c(interactionIpName, rep(ipName[i], 
            selLength))
        interactionBaitName <- c(interactionBaitName, rep(baitName[i], 
            selLength))
        interactionPreyName <- c(interactionPreyName, rownames(apmsExp$specCountDf)[sel])
        interactionSpecCount <- c(interactionSpecCount, apmsExp$specCountDf[sel, 
            ipName[i]])
    }
    interactionOut <- data.frame(interactionIpName, interactionBaitName, 
        interactionPreyName, interactionSpecCount)
    write.table(file = file, interactionOut, row.names = F, col.names = F, 
        sep = "\t", quote = F)
    if (verbose) {
        cat("Created File ", file, "\n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
