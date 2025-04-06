\name{toDOY}
\alias{toDOY}
\title{
Function for Transferring a Date to the Value of Day-of-Year
}
\description{
Transfers the date (from year, month and day) to the value of day-of-year.
}
\usage{
toDOY(Year, Month, Day)
}
\arguments{
  \item{Year}{the vector of years}
  \item{Month}{the vector of months}
  \item{Day}{the vector of days}
}
\details{
The user needs to provide the three separate vectors of \code{Year}, \code{Month} and \code{Day}, 
rather than providing a single date vector. The arguments can be numerical vectors or character vectors.
}
\value{
  The returned value is a vector of transferred dates in day-of-year.
}
\references{

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
The returned vector, \code{DOY}, usually matches with the year vector and the mean daily 
temperature vector as arguments in other functions, e.g., the \code{\link{ADD}} function.
}
\seealso{
\code{\link{BJDAT}}
}
\examples{

data(BJDAT)
X2   <- BJDAT
DOY2 <- toDOY(X2$Year, X2$Month, X2$Day)
# cbind(X2$DOY, DOY2)

}