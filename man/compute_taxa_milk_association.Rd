\name{compute_taxa_milkType_association}
\alias{compute_taxa_milkType_association}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute association between reads mapping to taxanomic groups and milk consumption type
%%  ~~function to do ... ~~
}
\description{
This function will output the results of multiple statistical models determining association between taxonomic groups of gut microbiota and milk consumption type at 6-months of age.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_taxa_milkType_association(asfbf_fdrFile=fdrCI_lefse_phylum_level_6m_asfbf,tfbf_fdrFile = fdrCI_lefse_phylum_level_6m_tfbf, taxaFile = phylumTaxa)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{asfbf_fdrFile}{
  R Data file containing the permuttated false discover rates for taxa differences between added-sugar formula group and breastfed groups.
%%     ~~Describe \code{x} here~~
}
  \item{tfbf_fdrFile}{
  R Data file containing the permuttated false discover rates for taxa differences between traditional formula group and breastfed groups.
%%     ~~Describe \code{x} here~~
}
  \item{taxaFile}{
  R Data file containing read counts of taxa at specified level with samples as columns and taxonomic groups as rows.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output will be a data frame with the following columns:
\item{allOutcomes }{Taxonomic group}
\item{bfList }{Mean (Standard Deviation) of taxonomic group in breastfed group}
\item{tfList }{Mean (Standard Deviation) of taxonomic group in traditional formula fed group}
\item{asfList }{Mean (Standard Deviation) of taxonomic group in added-sugar formula fed group}
\item{difftfbfList }{Mean difference (Lower confidence level, upper confidence level) between diversity metric in traditional formula and breastfed groups}
\item{pValTFBF_list }{P-value of  difference  between traditional formula group and breastfed group}
\item{diffasfbfList }{Mean difference (Lower confidence level, upper confidence level) between added-sugar formula group and breastfed group}
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
Roshonda B. Jones
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
compute_taxa_milkType_association()
compute_taxa_milkType_association(asfbf_fdrFile=fdrCI_lefse_genus_level_6m_asfbf,tfbf_fdrFile = fdrCI_lefse_genus_level_6m_tfbf, taxaFile = genusTaxa)
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
