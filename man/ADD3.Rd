 \name{ADD3}
\alias{ADD3}
\title{
Function for Implementing the Accumulated Degree Days Method Using Mean Daily Temperatures for 
the Combinations of the Starting Date and Base Temperature
}
\description{
Estimates the starting date (\eqn{S}, in day-of-year) and base temperature (\eqn{T_{0}}, 
in \eqn{{}^{\circ}}C) in the accumulated degree days (ADD) method using mean daily air temperatures 
(Konno and Sugihara, 1986; Aono, 1993; Shi et al., 2017a, b).
}
\usage{
ADD3( S.arr, T0.arr, Year1, Time, Year2, DOY, Temp, DOY.ul = 120, 
      fig.opt = TRUE, verbose = TRUE )
}
\arguments{
  \item{S.arr}{the candidate starting dates for thermal accumulation (in day-of-year)}
  \item{T0.arr}{the candidate base temperatures (in \eqn{{}^{\circ}}C)}
  \item{Year1}{the vector of the years in which a particular phenological event was recorded}
  \item{Time}{the vector of the occurrence times (in day-of-year) of a particular phenological event across many years}
  \item{Year2}{the vector of the years recording the climate data corresponding to the occurrence times}
  \item{DOY}{the vector of the dates (in day-of-year) for which climate data exist}
  \item{Temp}{the mean daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{DOY.ul}{the upper limit of \code{DOY} used to predict the occurrence time}
  \item{fig.opt}{an optional argument to draw the figures associated with the determination 
of the combination the starting date and base temperature, 
and a comparison between the predicted and observed occurrence times}
  \item{verbose}{an optional argument allowing users to suppress the printing of computation progress}
}
\details{

When \code{fig.opt} is equal to \code{TRUE}, it will show the contours of the root-mean-square 
errors (RMSEs) based on different combinations of \eqn{S} and \eqn{T_{0}}.

\eqn{\qquad}The function does not require that \code{Year1} is the same as \code{unique(Year2)}, 
and the intersection of the two vectors of years will be kept. The unused years that have phenological 
records but lack climate data will be showed in \code{unused.years} in the returned list. 

\eqn{\qquad}The numerical value of \code{DOY.ul} should be greater than or equal to the maximum \code{Time}.

}
\value{
  \item{mAADD.mat}{a matrix consisting of the means of the annual accumulated degree days (AADD) 
values from the combinations of \eqn{S} and \eqn{T_{0}}}
  \item{RMSE.mat}{the matrix consisting of the RMSEs (in days) from different 
combinations of \eqn{S} and \eqn{T_{0}}}
  \item{AADD.arr}{the AADD values in different years 
associated with the smallest value in \code{RMSE.mat}} 
  \item{Year}{The overlapping years between \code{Year1} and \code{Year2}}
  \item{Time}{The observed occurrence times (day-of-year) in the overlapping years 
between \code{Year1} and \code{Year2}}
  \item{Time.pred}{the predicted occurrence times in different years}
  \item{S}{the determined starting date (day-of-year)}
  \item{T0}{the determined base temperature (in \eqn{{}^{\circ}}C)} 
  \item{AADD}{the expected AADD}
  \item{RMSE}{the smallest RMSE (in days) in \code{RMSE.mat} from 
different combinations of \eqn{S} and \eqn{T_{0}}}
  \item{unused.years}{the years that have phenological records but lack climate data}
}
\references{
Aono, Y. (1993) Climatological studies on blooming of cherry tree (\emph{Prunus yedoensis}) by means 
of DTS method. \emph{Bulletin of the University of Osaka Prefecture. Ser. B, Agriculture and life sciences} 
45, 155\eqn{-}192 (in Japanese with English abstract).

Konno, T., Sugihara, S. (1986) Temperature index for characterizing biological activity in soil and 
its application to decomposition of soil organic matter. \emph{Bulletin of National Institute for 
Agro-Environmental Sciences} 1, 51\eqn{-}68 (in Japanese with English abstract).

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
S.arr0     <- seq(60, 70, by = 1)
T0.arr0    <- seq(-2, 5, by = 1)

\donttest{
  RES1 <- ADD3( S.arr = S.arr0, T0.arr = T0.arr0, Year1 = Year1.val, Time = Time.val, 
                Year2 = Year2.val, DOY = DOY.val, Temp = Temp.val, DOY.ul = DOY.ul.val, 
                fig.opt = TRUE, verbose = TRUE)
  RES1

  RMSE.mat0  <- RES1$RMSE.mat
  RMSE.range <- range(RMSE.mat0)

  dev.new()
  par1 <- par(family="serif")
  par2 <- par(mar=c(5, 5, 2, 2))
  par3 <- par(mgp=c(3, 1, 0))
  image( S.arr0, T0.arr0, RMSE.mat0, col = terrain.colors(200), axes = TRUE, 
         cex.axis = 1.5, cex.lab = 1.5, xlab = "Starting date (day-of-year)", 
         ylab = expression(paste("Base temperature (", degree, "C)", sep = "")))
  points( RES1$S, RES1$T0, cex = 1.5, pch = 16, col = 2 )
  contour( S.arr0, T0.arr0, RMSE.mat0, levels = round(seq(RMSE.range[1], 
           RMSE.range[2], len = 20), 4), add = TRUE, cex = 1.5, col = "#696969", labcex = 1.5)
  par(par1)
  par(par2)
  par(par3)

  # graphics.off()

}
}