\name{compute_possibleCovariates_milkType}
\alias{compute_possibleCovariates_milkType}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute p-values for ANOVA testing differences in means of milk consumption type group in general chacteristics of cohorts
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_possibleCovariates_milkType()


}
%- maybe also 'usage' for other objects documented here.
\arguments{
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output will be a data frame with the following columns:
\item{allOutcomes }{Name of covariate}
\item{bfList }{Mean (Standard Deviation) of breastfed group}
\item{tfList }{Mean (Standard Deviation) of traditional formula fed group}
\item{asfList }{Mean (Standard Deviation) of added-sugar formula fed group}
\item{pValSugarList }{P-value of  ANOVA testing differences in means of milk consumption type groups}
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
compute_possibleCovariates_milkType()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
