\name{BJDAT}
\docType{data}
\alias{BJDAT}
\title{
Daily Air Temperature Data of Beijing from 1952 to 2012. 
}
\description{
The data include the mean, minimum, and maximum daily temperatures (in \eqn{{}^{\circ}}C) 
of Beijing between 1952 and 2012. 
\bold{Data source}: China Meteorological Data Service Centre (\url{https://data.cma.cn/en}).
}
\usage{
data(BJDAT)
}
\details{
In the data set, there are seven columns of vectors: \code{Year}, \code{Month}, 
\code{Day}, \code{DOY}, \code{MDT}, \code{MinDT}, and \code{MaxDT}. 
\code{Year} saves the recording years; 
\code{Month} saves the recording months; 
\code{Day} saves the recording days; 
\code{DOY} saves the dates in day-of-year; 
\code{MDT} saves the mean daily temperatures (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}; 
\code{MinDT} saves the minimum daily temperatures (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}; 
\code{MaxDT} saves the maximum daily temperatures (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}.
}
\references{
Guo, L., Xu, J., Dai, J., Cheng, J., Wu, H., Luedeling, E. (2015) Statistical identification of chilling 
and heat requirements for apricotflower buds in Beijing, China. 
\emph{Scientia Horticulturae} 195, 138\eqn{-}144. \doi{10.1016/j.scienta.2015.09.006}
}
\examples{
data(BJDAT)
attach(BJDAT)

x    <- as.numeric( tapply(DOY, DOY, mean) )
y    <- as.numeric( tapply(MDT, DOY, mean) )
y.sd <- as.numeric( tapply(MDT, DOY, sd) )

dev.new()
par1 <- par(family="serif")
par2 <- par(mar=c(5, 5, 2, 2))
par3 <- par(mgp=c(3, 1, 0))
plot( x, y, cex = 1.5, xlim = c(0, 367), ylim = c(-10, 30), 
      cex.lab = 1.5, cex.axis = 1.5, type = "n", xlab = "Day-of-year", 
      ylab = expression(paste("Mean daily temperature (", degree, "C)", sep="")) )
for(i in 1:length(x)){
  lines(c(x[i], x[i]), c(y[i]-y.sd[i], y[i]+y.sd[i]), col=4)
}
points(x, y, cex = 1.5)
par(par1)
par(par2)
par(par3)

# graphics.off()
}
