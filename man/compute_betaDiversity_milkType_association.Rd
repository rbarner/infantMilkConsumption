\name{compute_betaDiversity_milkType_association}
\alias{compute_betaDiversity_milkType_association}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute association between Beta-Diversity and milk consumption type
%%  ~~function to do ... ~~
}
\description{
This function will output the results of multiple statistical models determining association between beta diversity of gut microbiota and milk consumption type at 6-months of age.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_betaDiversity_milkType_association()

}
%- maybe also 'usage' for other objects documented here.
\arguments{

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output will be a data frame with the following columns:
\item{allOutcomes }{Beta Diversity axes}
\item{bfList }{Mean (Standard Deviation) of beta diversity axes in breastfed group}
\item{tfList }{Mean (Standard Deviation) of beta diversity axes in traditional formula fed group}
\item{asfList }{Mean (Standard Deviation) of beta diversity axes in added-sugar formula fed group}
\item{difftfbfList }{Mean difference (Lower confidence level, upper confidence level) between beta diversity axes of traditional formula and breastfed groups}
\item{pValTFBF_list }{P-value of  difference  between traditional formula group and breastfed group}
\item{diffasfbfList }{Mean difference (Lower confidence level, upper confidence level) between beta diversity axes of added-sugar formula group and breastfed group}
\item{pValTFBF_list }{P-value of  difference  between added-sugar formula group and breastfed group}
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
compute_betaDiversity_milkType_association()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
