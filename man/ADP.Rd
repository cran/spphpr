\name{ADP}
\alias{ADP}
\title{
Function for Implementing the Accumulated Developmental Progress Method 
}
\description{
Estimates the starting date (\eqn{S} in day of year) and the parameters in 
a developmental rate model in the accumulated developmental progress (ADP) 
method using mean daily air temperatures (Wagner et al., 1984; Shi et al., 2017a, 2017b).
}
\usage{
ADP( S.arr, expr, ini.val, Year1, Time, Year2, DOY, Temp, DOY.ul = 120,
     fig.opt = TRUE, control = list(), verbose = TRUE )
}
\arguments{
  \item{S.arr}{the candidate starting dates for thermal accumulation (in day of year)}
  \item{expr}{a user-defined model that is used in the accumulated developmental progress (ADP) method}
  \item{ini.val}{a vector or a list that saves the initial values of the parameters in \code{expr}}
  \item{Year1}{the vector of the years recording a particular phenological event}
  \item{Time}{the vector of the occurence times (in day of year) of a particular phenological event across many years}
  \item{Year2}{the vector of the years recording the climate data corresponding to the occurrence times}
  \item{DOY}{the vector of the dates (in day of year) when the climate data exist}
  \item{Temp}{the mean daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{DOY.ul}{the upper limit of \code{DOY} used to predict the occurrence time}
  \item{fig.opt}{an optional argument of drawing the figures assoicated with the temperature-dependent developmental rate curve, 
the mean daily temperatures versus years, and a comparision between the predicted and observed occurrence times}
  \item{control}{the list of control parameters for using the \code{\link{optim}} function in package \pkg{stats}}
  \item{verbose}{an optional argument of allowing users to suppress printing of computation progress}
}
\details{

It is better not to set too much candiate starting dates, which will be time-consuming. If \code{expr} 
is selected as Arrhenius' equation, \code{S.arr} can be selected as \code{S} obtained from the output of 
carrying out the \code{\link{ADTS}} function. Here, \code{expr} can be other nonlinear temperature-dependent 
developmental rate functions (see Shi et al. [2017b] for details). Here, \code{expr} can be any an arbitrary 
user-defined temperature-dependent developmental rate function, e.g., a function named \code{myfun}, 
but it needs to take the following form of \code{myfun <- function(P, x){...}}, 
where \code{P} is the vector of the model parameter(s), and \code{x} is the vector of the 
predictor variable, i.e., the temperature variable.

\eqn{\qquad}The function does not require that \code{Year1} is the same as the unique of \code{Year2}, 
and the intersection of two years will be finally kept. The unused years that have phenological 
records but lack the climate data will be showed in \code{unused.years} in the returned list. 

\eqn{\qquad}The numerical value of \code{DOY.ul} should be larger than or equal to the maximum \code{Time}.

\eqn{\qquad} Let \eqn{r} represent the temperature-dependent developmental rate, i.e., 
the reciprocal of the developmental 
duration at a constant temperature required for completing a particular phenological event. 
In the accumulated developmental progress (ADP) method, when the annual accumulated developmental 
progress (AADP) reaches 100\%, the phenological event is predicted to occurr for each year. 
Let \eqn{\mathrm{AADP}_{i}} denote the AADP of the \eqn{i}th year, which equals

\deqn{\mathrm{AADP}_{i} = \sum_{j=S}^{E_{i}}r_{ij}\left(\mathrm{\mathbf{P}}; T_{ij}\right),}

where \eqn{S} represents the starting date (in day of year), \eqn{E_{i}} represents the ending date 
(in day of year), i.e., the occurrence time of a pariticular phenological event in the \eqn{i}th year, 
\eqn{\mathrm{\mathbf{P}}} is the vector of the model parameters in \code{expr}, 
and \eqn{T_{ij}} represents the mean daily temperature of the \eqn{j}th day of the \eqn{i}th 
year (in \eqn{{}^{\circ}}C or K). In theory, \eqn{\mathrm{AADP}_{i} = 100\%}, 
i.e., the AADP values of different years are a constant of 100\%. However, in practice, there is 
a certain deviation of \eqn{\mathrm{AADP}_{i}} from 100\%. The following approach 
is used to determine the predicted occurrence time. 
When \eqn{\sum_{j=S}^{F}r_{ij} = 100\%} (where \eqn{F \geq S}), it follows that \eqn{F} is 
the predicted occurrence time; when \eqn{\sum_{j=S}^{F}r_{ij} < 100\%} and 
\eqn{\sum_{j=S}^{F+1}r_{ij} > 100\%}, the trapezoid method (Ring and Harris, 1983) 
is used to determine the predicted occurrence time. Assume that there are \eqn{n}-year phenological records. 
When the starting date \eqn{S} and the temperature-dependent developmental rate model are known, 
the model parameters can be estimated using the Nelder-Mead optiminization method 
(Nelder and Mead, 1965) to minimize the root-mean-square error (RMSE) between the observed and predicted 
occurrence times, i.e.,

\deqn{\mathrm{\mathbf{\hat{P}}} = \mathrm{arg}\,\mathop{\mathrm{min}}\limits_{
\mathrm{\mathbf{P}}}\left\{\mathrm{RMSE}\right\} = \mathrm{arg}\,\mathop{\mathrm{min}}\limits_{
\mathrm{\mathbf{P}}}\sqrt{\frac{\sum_{i=1}^{n}\left(E_{i}-\hat{E}_i\right){}^{2}}{n}}.}

Because \eqn{S} is not determined, a group of candidate \eqn{S} values (in day of year) need to 
be provided. Assume that there are \eqn{m} candidate \eqn{S} values, i.e., \eqn{S_{1}, S_{2}, S_{3}, 
\cdots, S_{m}}. For each \eqn{S_{q}} (where \eqn{q} ranges between 1 and \eqn{m}), we can obtain 
a vector of the estimated model parameters, \eqn{\mathrm{\mathbf{\hat{P}}}_{q}}, by minimizing 
\eqn{\mathrm{RMSE}_{q}} using the Nelder-Mead optiminization method. Then we finally selected 
\eqn{\mathrm{\mathbf{\hat{P}}}} associated with \eqn{\mathrm{min}\left\{\mathrm{RMSE}_{1}, 
\mathrm{RMSE}_{2}, \mathrm{RMSE}_{3}, \cdots, \mathrm{RMSE}_{m}\right\}} 
as the target parameter vector.      
}
\value{
  \item{TDDR}{the temperature-dependent developmental rate matrix consisting of the year, 
day of year, mean daily temperature and developmental rate columns}
  \item{MAT}{a matrix consisting of the candidate starting dates, the estimates of 
candidate model parameters with the corresponding RMSEs}
  \item{Dev.accum}{the calculated annual accumulated developmental progresses in different years} 
  \item{Year}{The intersected years between \code{Year1} and \code{Year2}}
  \item{Time}{The observed occurence times (day of year) in the intersected years 
between \code{Year1} and \code{Year2}}
  \item{Time.pred}{the predicted occurence times in different years}
  \item{S}{the determined starting date (day of year)}
  \item{par}{the estiamtes of model parameters}
  \item{RMSE}{the RMSE (in days) between the observed and predicted occurrence times}
  \item{unused.years}{the years that have phenological records but lack the climate data}
}
\references{
Nelder, J.A., Mead, R. (1965) A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Ring, D.R., Harris, M.K. (1983) Predicting pecan nut casebearer (Lepidoptera: Pyralidae) activity 
at College Station, Texas. \emph{Environmental Entomology} 12, 482\eqn{-}486. \doi{10.1093/ee/12.2.482}

Shi, P., Chen, Z., Reddy, G.V.P., Hui, C., Huang, J., Xiao, M. (2017a) Timing of cherry tree blooming: 
Contrasting effects of rising winter low temperatures and early spring temperatures. 
\emph{Agricultural and Forest Meteorology} 240\eqn{-}241, 78\eqn{-}89. \doi{10.1016/j.agrformet.2017.04.001}

Shi, P., Fan, M., Reddy, G.V.P. (2017b) Comparison of thermal performance equations in describing 
temperature-dependent developmental rates of insects: (III) Phenological applications.  
\emph{Annals of the Entomological Society of America} 110, 558\eqn{-}564. \doi{10.1093/aesa/sax063}

Wagner, T.L., Wu, H.-I., Sharpe, P.J.H., Shcoolfield, R.M., Coulson, R.N. (1984) Modelling insect 
development rates: a literature review and application of a biophysical model. 
\emph{Annals of the Entomological Society of America} 77, 208\eqn{-}225. \doi{10.1093/aesa/77.2.208}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Zhenghong Chen \email{chenzh64@126.com}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\note{
The entire mean daily temperature data in the spring of each year should be provided. 
In \code{TDDR}, the first column of \code{Year} saves the years, the second column 
of \code{DOY} saves the day of year values, the third column of \code{Temperature} 
saves the mean daily air temperatures between the starting date to the occurrence times, 
and the fourth column of \code{Rate} saves the calculated developmental rates 
corresponding to the mean daily temperatures. 
}
\seealso{
\code{\link{predADP}}
}
\examples{

data(apricotFFD)
data(BJMDT)
X1 <- apricotFFD
X2 <- BJMDT

Year1.val  <- X1$Year
Time.val   <- X1$Time
Year2.val  <- X2$Year
DOY.val    <- X2$DOY
Temp.val   <- X2$MDT
DOY.ul.val <- 120
S.arr0     <- 47

#### Defines a re-parameterized Arrhenius' equation ########################### 
Arrhenius.eqn <- function(P, x){
  B  <- P[1]
  Ea <- P[2]
  R  <- 1.987 * 10^(-3)
  x  <- x + 273.15
  10^12*exp(B-Ea/(R*x))
}
##############################################################################

#### Provides the initial values of the parameter of Arrhenius' equation #####
ini.val0 <- list( B = 20, Ea = 14 )
##############################################################################

\donttest{
  res5 <- ADP( S.arr = S.arr0, expr = Arrhenius.eqn, ini.val = ini.val0, Year1 = Year1.val, 
               Time = Time.val, Year2 = Year2.val, DOY = DOY.val, Temp = Temp.val, 
               DOY.ul = DOY.ul.val, fig.opt = TRUE, control = list(trace = FALSE, 
               reltol = 1e-12, maxit = 5000), verbose = TRUE )
  res5

  TDDR <- res5$TDDR
  T    <- TDDR$Temperature
  r    <- TDDR$Rate
  Y    <- res5$Year
  DP   <- res5$Dev.accum

  dev.new()
  par1 <- par(family="serif")
  par2 <- par(mar=c(5, 5, 2, 2))
  par3 <- par(mgp=c(3, 1, 0))
  Ind <- sort(T, index.return=TRUE)$ix 
  T1  <- T[Ind]
  r1  <- r[Ind]
  plot( T1, r1, cex.lab = 1.5, cex.axis = 1.5, pch = 1, cex = 1.5, col = 2, type = "l", 
        xlab = expression(paste("Mean daily temperature (", degree, "C)", sep = "")), 
        ylab = expression(paste("Calculated developmental rate (", {day}^{"-1"}, ")", sep = "")) ) 
  par(par1)
  par(par2)
  par(par3)

  dev.new()
  par1 <- par(family="serif")
  par2 <- par(mar=c(5, 5, 2, 2))
  par3 <- par(mgp=c(3, 1, 0))
  plot( Y, DP * 100, xlab = "Year", 
        ylab = "Accumulated developmental progress (\%)", 
        ylim = c(50, 150), cex.lab=1.5, cex.axis = 1.5, cex = 1.5 )
  abline( h = 1 * 100, lwd = 1, col = 4, lty = 2 ) 
  par(par1)
  par(par2)
  par(par3)

  # graphics.off()
}
}