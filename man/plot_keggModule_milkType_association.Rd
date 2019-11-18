\name{plot_keggModule_milkType_association}
\alias{plot_keggModule_milkType_association}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate boxplots of reads mapping to predicted KEGG modules of gut microbiota by milk consumption type
%%  ~~function to do ... ~~
}
\description{
This function will output a boxplot with means of a specified KEGG (Kyoto Encyclopedia of Genes and Genomes) module and milk consumption type at 6-months of age.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_keggModule_milkType_association(variable = "M00221", variableName = "M00221: Putative simple sugar transport system")


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variable}{
 KEGG module to be plotted.
%%     ~~Describe \code{x} here~~
}
  \item{variableName}{
  Name of KEGG module to be plotted.
%%     ~~Describe \code{x} here~~
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
Roshonda B. Jones
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
plot_keggModule_milkType_association()
plot_keggModule_milkType_association(variable="M00233",variableName="M00233: Glutamate transport system")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
