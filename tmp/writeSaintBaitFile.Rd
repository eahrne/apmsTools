\name{writeSaintBaitFile}
\alias{writeSaintBaitFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeSaintBaitFile(apmsExp, file, verbose = F)
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
    ipName <- names(apmsExp$specCountDf)[4:ncol(apmsExp$specCountDf)]
    baitName <- gsub("\\-[0-9]{1,}$", "", ipName)
    isControlStr <- c("T", "C")[(regexpr("^GFP", ipName) > -1) + 
        1]
    baitOut <- data.frame(ipName, baitName, isControlStr)
    write.table(baitOut, file = file, row.names = F, col.names = F, 
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
