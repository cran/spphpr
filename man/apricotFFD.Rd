\name{apricotFFD}
\docType{data}
\alias{apricotFFD}
\title{
First flowering date records of \emph{Prunus armeniaca}
}
\description{
The data consist of the first flowering date records of \emph{Prunus armeniaca} 
at the Summer Palace (39\eqn{{}^{\circ}}54\eqn{'}38\eqn{''} N, 
116\eqn{{}^{\circ}}8\eqn{'}28\eqn{''} E, 50 m a.s.l.) in Beijing of China 
between 1963 and 2010. \bold{Data source}: Chinese Phenological Observation Network (Guo et al., 2015).

}
\usage{
data(apricotFFD)
}
\details{
In the data set, there are two columns of vectors: \code{Year}, and \code{Time}. 
\code{Code} saves the recording years; and 
\code{Time} saves 1963\eqn{-}2010 first flowering dates of \emph{Prunus armeniaca} in day of year.
}
\references{
Guo, L., Xu, J., Dai, J., Cheng, J., Wu, H., Luedeling, E. (2015) Statistical identification of chilling 
and heat requirements for apricotflower buds in Beijing, China. 
\emph{Scientia Horticulturae} 195, 138\eqn{-}144. \doi{10.1016/j.scienta.2015.09.006}
}
\examples{
data(apricotFFD)
attach(apricotFFD)

dev.new()
par1 <- par(family="serif")
par2 <- par(mar=c(5, 5, 2, 2))
par3 <- par(mgp=c(3, 1, 0))
plot( Year, Time, asp = 1, cex.lab = 1.5, cex.axis = 1.5, 
      xlab = "Year", ylab = "First flowering date (day of year)" ) 
par(par1)
par(par2)
par(par3)

# graphics.off()
}
