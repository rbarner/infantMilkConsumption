\name{plot_alphaDiversity_milkType_association}
\alias{plot_alphaDiversity_milkType_association}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create boxplots of alpha diversity and milk consumption type
%%  ~~function to do ... ~~
}
\description{
This function will output the a boxplot of means of  alpha diversity of gut microbiota across milk consumption types at 6-months of age.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_alphaDiversity_milkType_association(variable = "shannonDiversity",variableName = "Shannon Diversity")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variable}{
  Diversity index to be plotted. Options are "shannonDiversity","simpsonDiversity","richness","evenness","numSeqPerSample"
%%     ~~Describe \code{x} here~~
}
  \item{variableName}{
  Name of diversity index to be plotted. Options are "Shannon Diversity","Simpson Diversity","Richness","Evenness","Number of Sequences Per Sample"
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\examples{
plot_alphaDiversity_milkType_association()
plot_alphaDiversity_milkType_association(variable = "simpsonDiversity",variableName = "Simpson Diversity")
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
