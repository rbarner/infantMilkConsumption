\name{plot_somaticGrowth_milkType_association}
\alias{plot_somaticGrowth_milkType_association}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create boxplots of child growth at 24 months by milk consumption type
%%  ~~function to do ... ~~
}
\description{
This function will output a boxplot with means of child growth measure at 24-months of age by milk consumption type at 6-months of age.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_somaticGrowth_milkType_association(variable = "inf_weight_kg",variableName = "Infant Weight")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variable}{
  Growth measure to be plotted. Options are "inf_weight_kg","zbmi","zwfl","zwei","zlen","skinf_tricep_mm","skinf_subscap_mm","skinf_supra_mm","skinf_midthigh_mm","circ_umb_cm"
%%     ~~Describe \code{x} here~~
}
  \item{variableName}{
  Name of growth measure to be plotted. Options are "Infant Weight","BMI Z-score","Weight-for-Length Z-score","Weight Z-score","Length Z-score","Tricep Skinfold (mm)","Subscapular Skinfold (mm)","Suprailiac Skinfold (mm)","Midthigh Skinfold (mm)","Abdominal Circumference (cm)"
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
plot_somaticGrowth_milkType_association()
plot_somaticGrowth_milkType_association(variable = "zbmi",variableName = "BMI Z-score")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
