\name{PRFplot}
\alias{PRFplot}
\title{Person response function (PRF)}
\description{Plot the nonparametric person response function with variability bands.}
\usage{
PRFplot(x, respID, h=.09, N.FPts=101, 
        VarBands=FALSE, VarBands.area=FALSE, alpha=.05,
        Xlabel=NA, Xcex=1.5, Ylabel=NA, Ycex=1.5, title=NA, Tcex=1.5)
}

\arguments{
  \item{x}{Object of class "PerFit".}
  \item{respID}{Vector specifying the respondents for whom PRFs are to be computed.}
  \item{h}{Bandwidth value. Default is 0.09.}
  \item{N.FPts}{Number of (equidistant) focal points in the [0,1] interval. Default is 101.}
  \item{VarBands}{Logical: Draw the \code{(1-alpha)} variability bands? Default is \code{FALSE}.}
  \item{VarBands.area}{Logical: Draw the area between the \code{(1-alpha)} variability bands? Default is \code{FALSE}.}
  \item{alpha}{Significance level for the variability bands. Default is 0.05.}
  \item{Xlabel}{Define label of x-axis, otherwise a default label is shown.}
  \item{Xcex}{Font size of the label of x-axis.}
  \item{Ylabel}{Define label of y-axis, otherwise a default label is shown.}
  \item{Ycex}{Font size of the label of y-axis.}
  \item{title}{Define the title of the plot, otherwise a default title is shown.}
  \item{Tcex}{Font size of the title of the plot.}
}

\details{
Function \code{PRFplot} displays the so-called nonparametric person response functions (PRFs; Emons, Sijtsma, and Meijer, 2004; Sijtsma and Meijer, 2001),  for the respondents identified in \code{respID}. The PRF relates item difficulty (0-1 range on the x-axis) with the associated probability of correct response (on the y-axis). The PRF is typically nonincreasing, implying that the probability of answering increasingly difficult items should (at least) not increase. The code is based on nonparametric kernel smoothing (Emons et al., 2004). The value of the PRF at each focal point (representing a difficulty parameter between 0 and 1) is estimated as a weighted sum score, where scores pertaining to items with difficulty close to the focal point are given the largest weights. The weights are functions of the Gaussian kernel function. It is necessary to specify a bandwidth value (\code{h}) in order to compute the weights. The \code{h} value controls the trade-off between bias and sampling variation (Emons et al., 2004). Small \code{h} values reduce bias but increase variance, leading to PRFs that capture too much measurement error. Large \code{h} values, on the other hand, increase bias which renders PRFs that are often too flat, thus missing potentially relevant misfitting response behavior. Therefore, it is important to carefuly specify the value \code{h}. Emons et al. (2004, pp. 10-13), after a simulation study, advised that "\code{h} values between 0.07 and 0.11 are reasonable choices".

Moreover, variability bands of level \code{alpha} (0.05 by default) can also be added to the plot. These bands are computed following the jackknife procedure explained in Emons et al. (2004). 

The PRFs and variability bands for each respondent are approximated by means of functional data objects (e.g., Ramsay, Hooker, and Graves, 2009), with the help of the \code{fda} package. This procedure follows two steps:
\itemize{
\item Compute a B-splines basis system. This basis consists of a set of (thirteen) piecewise polinomials, all of degree three/order four (i.e., cubic polinomial segments), with one knot per break point. This allows any two consecutive splines, sp1 and sp2, with common break point BP, verifying sp1(BP) = sp2(BP), sp1'(BP) = sp2'(BP), and sp1''(BP) = sp2''(BP). At 0 and 1 (extremes of the x-range), four (= order) knots are used.
\item Specify coefficients c for the B-splines basis system computed above and then create functional data objects. Based on smoothing using regression analysis (Ramsay et al., 2009, section 4.3).
}

\code{PRFplot} returns three functional data objects (for the PRFs, lower-bound of the variability bands, and upper-band of the variability bands) for all respondents in the sample. 
}

\value{
The output is a list with three functional data objects of class \code{fd} (see package \code{fda}):
  \item{\code{PRF.FDO}}{Functional data object of the PRFs for the entire sample.}
  \item{\code{VarBandsLow.FDO}}{Functional data object of the lower-bound of the variability bands for the entire sample.}
  \item{\code{VarBandsHigh.FDO}}{Functional data object of the upper-bound of the variability bands for the entire sample.}
}

\references{
Emons, W. M., Sijtsma, K., and Meijer, R. R. (2004) Testing hypotheses about the person-response function in person-fit analysis. \emph{Multivariate Behavioral Research}, \bold{39(1)}, 1--35.

Ramsay, J. O., Hooker, G., and Graves, S. (2009) \emph{Functional data analysis with R and MATLAB}. New York: US.

Sijtsma, K., and Meijer, R. R. (2001) The person response function as a tool in person-fit research. \emph{Psychometrika}, \bold{66(2)}, 191--207.
}

\author{Jorge N. Tendeiro \email{j.n.tendeiro@rug.nl}}
\seealso{\code{\link{cutoff}}, \code{\link{plot.PerFit}}, \code{\link{flagged.resp}}}
\examples{
# Load the inadequacy scale data (dichotomous item scores):
data(InadequacyData);

# As an example, compute the Ht person-fit scores:
Ht.out <- Ht(InadequacyData)
Ht.out$PFscores

# Determine which respondents were flagged by Ht at 1% level:
Ht.flagged <- flagged.resp(Ht.out, Blvl=.01, scores=FALSE)
Ht.flagged <- Ht.flagged$PFSscores[,1]
# Flagged respondents: 30  37  46  49 137 158 165 181 216 243 286 303 
#         397 401 531 549 559 582 591 619 750

# Plot the PRFs of the first three flagged respondents:
Flagged    <- PRFplot(Ht.out, respID=Ht.flagged[1:3])
# Plot the person response function of respondent 35 (not flagged as aberrant):
NotFlagged <- PRFplot(Ht.out, respID=35)
# Plot the PRFs of all respondents:
plot(Flagged$PRF.FDO)
}
\keyword{univar}

