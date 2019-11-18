\name{plot_bayley_milkType_association}
\alias{plot_bayley_milkType_association}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create boxplots of child development at 24 months and milk consumption type
%%  ~~function to do ... ~~
}
\description{
This function will output a boxplot with means of child development scaled score and at 24-months of age by milk consumption type at 6-months of age.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_bayley_milkType_association(variable = "bsid_cog_ss",variableName = "Cognitive Scaled Score")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variable}{
  Development outcome to be plotted. Options are "bsid_cog_ss","bsid_lang_rc_ss","bsid_lang_ec_ss","bsid_lang_ss","bsid_mot_fm_ss","bsid_mot_gm_ss","bsid_mot_ss","bsid_se_ss","bsid_ab_com_ss","bsid_ab_cu_ss","bsid_ab_fa_ss","bsid_ab_hl_ss","bsid_ab_hs_ss","bsid_ab_ls_ss","bsid_ab_sc_ss","bsid_ab_sd_ss","bsid_ab_soc_ss","bsid_ab_mo_ss","bsid_ab_ss"
%%     ~~Describe \code{x} here~~
}
  \item{variableName}{
  Name of devlopmental outcome to be plotted. Options are "Cognitive Scaled Score","Receptive Communication (RC) Total Scaled Score","Expressive Communication (EC) Total Scaled Score","Language Sum Scaled Score","Fine Motor (FM) Scaled Score","Gross Motor (GM) Scaled Score","Motor Sum Scaled Score","Social-Emotional Scaled Score","Communication (Com) Scaled Score","Community Use (CU) Scaled Score","Functional Pre-Academics (FA) Scaled Score","Home Living (HL) Scaled Score","Health and Safety (HS) Scaled Score","Leisure (LS) Scaled Score","Self-Care (SC) Scaled Score","Self-Direction (SD) Scaled Score","Social (Soc) Scaled Score","Motor (MO) Scaled Score","Adaptive Behavior Sum Scaled Score"
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
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
plot_bayley_milkType_association()
plot_bayley_milkType_association(variable = "bsid_mot_fm_ss",variableName = "Fine Motor (FM) Scaled Score")
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
