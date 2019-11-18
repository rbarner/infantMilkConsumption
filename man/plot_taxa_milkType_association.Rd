\name{plot_taxa_milkType_association}
\alias{plot_taxa_milkType_association}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create boxplots of specified taxonomic group by milk consumption type
%%  ~~function to do ... ~~
}
\description{
This function will output the a boxplot of means of a specified taxonomic group across milk consumption types at 6-months of age.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_taxa_milkType_association <- function(taxaFile = phylumTaxa,variable = "k__Bacteria.p__Actinobacteria",variableName = "Logged Abundance of Actinobacteria")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxaFile}{
  R Data file containing read counts of taxa at specified level with samples as columns and taxonomic groups as rows.
%%     ~~Describe \code{x} here~~
}
  \item{variable}{
  Taxa to be plotted.
%%     ~~Describe \code{x} here~~
}
  \item{variableName}{
  Name of taxa to be plotted. Will be used as the title of the y-axis.
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
plot_taxa_milkType_association()
plot_taxa_milkType_association(taxaFile = genusTaxa,variable ="p__Firmicutes.c__Clostridia.o__Clostridiales.f__Lachnospiraceae.g__Blautia",variableName="Logged Abundance of Blautia")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
