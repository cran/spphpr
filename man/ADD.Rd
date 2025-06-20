\name{ADD}
\alias{ADD}
\title{
Function for Implementing the Accumulated Degree Days Method Using Mean Daily Temperatures
}
\description{
Estimates the starting date (\eqn{S}, in day-of-year) and base tempeature (\eqn{T_{0}}, in \eqn{{}^{\circ}}C) 
in the accumulated degree days method using mean daily air temperatures (Aono, 1993; Shi et al., 2017a, b).
}
\usage{
ADD(S.pd = NULL, T0.arr, Year1, Time, Year2, DOY, Temp, DOY.ul = 120, 
    fig.opt = TRUE, S.def = 54, verbose = TRUE)
}
\arguments{
  \item{S.pd}{the pre-determined starting date for thermal accumulation (in day-of-year)}
  \item{T0.arr}{the candidate base temperatures (in \eqn{{}^{\circ}}C)}
  \item{Year1}{the vector of the years in which a particular phenological event was recorded}
  \item{Time}{the vector of the occurrence times (in day-of-year) of a particular phenological event across many years}
  \item{Year2}{the vector of the years recording the climate data corresponding to the occurrence times}
  \item{DOY}{the vector of the dates (in day-of-year) for which climate data exist}
  \item{Temp}{the mean daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{DOY.ul}{the upper limit of \code{DOY} used to predict the occurrence time}
  \item{fig.opt}{an optional argument to draw the figures associated with the determinations 
of the starting date and base temperature, and a comparison between the predicted and observed occurrence times}
  \item{S.def}{a mandatory defintion of the starting date when (i) \code{S.pd} is NULL 
and (ii) the minimum correlation coefficient method fails to find a suitable starting date}
  \item{verbose}{an optional argument allowing users to suppress the printing of computation progress}
}
\details{

The default of \code{S.pd} is NULL. In this case, the date associated with the minimum 
correlation coefficient [between the mean of the mean daily temperatures (from a candidate starting 
date to the observed occurrence time) and the observed occurrence time] 
will be determined to be the starting date on the condition that it is smaller than 
the mininum phenological occurrence time. If the determined date associated with the minimum 
correlation coefficient is greater than the mininum phenological occurrence time, \code{S.def} 
will be used as the starting date. If \code{S.pd} is not NULL, the starting date will be 
directly set as \code{S.pd} irrespective of the minimum correlation coefficient method 
and the value of \code{S.def}. This means that \code{S.pd} is superior to \code{S.def} in determining 
the starting date.

\eqn{\qquad}The function does not require that \code{Year1} is the same as \code{unique(Year2)}, 
and the intersection of the two vectors of years will be kept. The unused years that have phenological 
records but lack climate data will be showed in \code{unused.years} in the returned list. 

\eqn{\qquad}The numerical value of \code{DOY.ul} should be greater than or equal to the maximum \code{Time}.

}
\value{
  \item{S.arr}{the candidate starting dates (in day-of-year), whose default ranges from 
the minimum \code{DOY} to min(\code{DOY.ul}, the maximum \code{DOY})}
  \item{cor.coef.arr}{the candidate correlation coefficients between the mean of the mean daily tempertures 
(from a candidate starting date to the observed occurrence time) and the observed occurrence time}
  \item{cor.coef}{the minimum correlation coefficient, i.e., min(\code{cor.coef.arr})}
  \item{search.failure}{a value of 0 or 1 of showing whether the starting date is successfully determined 
by the minimum correlation coefficient method when \code{S.pd = NULL}, where 0 represents success and 1 represents failure}
  \item{mAADD.arr}{an vector saving the interannual mean of the annual acccumulated degree days 
(AADD) values for each of the candidate base temperatures}
  \item{RMSE.arr}{a vector saving the candidate root-mean-square errors (in days) between 
the observed and predicted occurrence times for each of the candidate base temperatures}
  \item{AADD.arr}{the annual accumulated degree days (AADD) values in different years}
  \item{Year}{The overlapping years between \code{Year1} and \code{Year2}}
  \item{Time}{The observed occurrence times (day-of-year) in the overlapping years 
between \code{Year1} and \code{Year2}}
  \item{Time.pred}{the predicted occurrence times in different years}
  \item{S}{the determined starting date (day-of-year)}
  \item{T0}{the determined base temperature (in \eqn{{}^{\circ}}C)}
  \item{AADD}{the expected annual accumulated degree days}
  \item{RMSE}{the smallest RMSE (in days) from the different candidate base temperatures}
  \item{unused.years}{the years that have phenological records but lack climate data}
}
\references{
Aono, Y. (1993) Climatological studies on blooming of cherry tree (\emph{Prunus yedoensis}) by means 
of DTS method. \emph{Bulletin of the University of Osaka Prefecture. Ser. B, Agriculture and life sciences} 
45, 155\eqn{-}192 (in Japanese with English abstract).

Shi, P., Chen, Z., Reddy, G.V.P., Hui, C., Huang, J., Xiao, M. (2017a) Timing of cherry tree blooming: 
Contrasting effects of rising winter low temperatures and early spring temperatures. 
\emph{Agricultural and Forest Meteorology} 240\eqn{-}241, 78\eqn{-}89. \doi{10.1016/j.agrformet.2017.04.001}

Shi, P., Fan, M., Reddy, G.V.P. (2017b) Comparison of thermal performance equations in describing 
temperature-dependent developmental rates of insects: (III) Phenological applications.  
\emph{Annals of the Entomological Society of America} 110, 558\eqn{-}564. \doi{10.1093/aesa/sax063}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Zhenghong Chen \email{chenzh64@126.com},
Jing Tan \email{jmjwyb@163.com}, Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\note{
The entire mean daily temperature data set for the spring of each year should be provided. 
\code{AADD} is represented by the mean of \code{AADD.arr} in the output. 
When the argument of \code{S.pd} is not \code{NULL}, the returned value of \code{search.failure} will be \code{NA}. 
When the argument of \code{S.pd} is \code{NULL}, and the minimum correlation coefficient method fails 
to find a suitable starting date, the argument of \code{S.def} is then defined as the determined starting 
date, i.e., the returned value of \code{S}. At the same time, the returned value of \code{cor.coef} is defined as \code{NA}.   

}
\seealso{
\code{\link{predADD}}
}
\examples{

data(apricotFFD)
data(BJDAT)
X1 <- apricotFFD
X2 <- BJDAT

Year1.val  <- X1$Year
Time.val   <- X1$Time
Year2.val  <- X2$Year
DOY.val    <- X2$DOY
Temp.val   <- X2$MDT
DOY.ul.val <- 120
T0.arr0    <- seq(-5, 5, by = 0.1) 
S.pd0      <- NULL

\donttest{
  res1 <- ADD( S.pd = S.pd0, T0.arr = T0.arr0, Year1 = Year1.val, Time = Time.val, 
               Year2 = Year2.val, DOY = DOY.val, Temp = Temp.val, 
               DOY.ul = DOY.ul.val, fig.opt = TRUE, S.def=54, verbose = TRUE )
  res1

  S0 <- res1$S.arr 
  r0 <- res1$cor.coef.arr
  
  dev.new()
  par1 <- par(family="serif")
  par2 <- par(mar=c(5, 5, 2, 2))
  par3 <- par(mgp=c(3, 1, 0))
  plot( S0, r0, cex.lab = 1.5, cex.axis = 1.5, xlab = "Candidate starting date (day-of-year)", 
        ylab="Correlation coefficient between the mean temperature and FFD", type="l" )
  ind <- which.min(r0)
  points(S0[ind], r0[ind], cex = 1.5, pch = 16)
  text(S0[ind], r0[ind] + 0.1, bquote(paste(italic(S), " = ", .(S0[ind]), sep = "")), cex = 1.5)
  par(par1)
  par(par2)
  par(par3)

  resu1 <- ADD( S.pd = 65, T0.arr = seq(-10, 0, by = 0.1), Year1 = Year1.val, Time = Time.val, 
                Year2 = Year2.val, DOY = DOY.val, Temp = Temp.val, 
                DOY.ul = DOY.ul.val, fig.opt = TRUE, S.def = 54, verbose = TRUE )
  resu1

  # graphics.off()

}
}