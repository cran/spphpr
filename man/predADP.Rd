\name{predADP}
\alias{predADP}
\title{
Prediction Function of the Accumulated Developmental Progress Method
}
\description{
Predicts the occurrence times using the accumulated developmental progress (ADP) 
method based on observed or predicted mean daily air temperatures 
(Wagner et al., 1984; Shi et al., 2017a, 2017b).
}
\usage{
predADP(S, expr, theta, Year2, DOY, Temp, DOY.ul = 120)
}
\arguments{
  \item{S}{the starting date for thermal accumulation (in day of year)}
  \item{expr}{a user-defined model that is used in the accumulated developmental progress (ADP) method}
  \item{theta}{a vector saves the numerical values of the parameters in \code{expr}}
  \item{Year2}{the vector of the years recording the climate data for predicting the occurrence times}
  \item{DOY}{the vector of the dates (in day of year) when the climate data exist}
  \item{Temp}{the mean daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{DOY.ul}{the upper limit of \code{DOY} used to predict the occurrence time}
}
\details{
Organisms showing phenological events in early spring often experience several cold days 
during the development. In this case, Arrhenius' equation (Shi et al., 2017a, 2017b, 
and references therein) has been recommended to describe the effect of the absolute temperature 
(\eqn{T} in Kelvin [K]) on the developmental rate (\eqn{r}):
\deqn{r = \mathrm{exp}\left(B - \frac{E_{a}}{R\,T}\right),}
where \eqn{E_{a}} represents the activation free energy (in kcal \eqn{\cdot} mol\eqn{{}^{-1}}); 
\eqn{R} is the universal gas constant (= 1.987 cal \eqn{\cdot} mol\eqn{{}^{-1}} \eqn{\cdot} K\eqn{{}^{-1}}); 
\eqn{B} is a constant. To keep the consistence of the unit used in \eqn{E_{a}} and \eqn{R}, we need to 
re-assign \eqn{R} to be 1.987\eqn{\times {10}^{-3}} to make its unit 1.987\eqn{\times {10}^{-3}} 
kcal \eqn{\cdot} mol\eqn{{}^{-1}} \eqn{\cdot} K\eqn{{}^{-1}} in the above formula. 

\eqn{\qquad}In the accumulated developmental progress (ADP) method, when the annual accumulated developmental 
progress (AADP) reaches 100\%, the phenological event is predicted to occurr for each year. 
Let \eqn{\mathrm{AADP}_{i}} denote the AADP of the \eqn{i}th year, which equals

\deqn{\mathrm{AADP}_{i} = \sum_{j=S}^{E_{i}}r_{ij},}

where \eqn{E_{i}} represents the ending date (in day of year), i.e., the occurrence time of a pariticular 
phenological event in the \eqn{i}th year. If the temperature-dependent developmental rate follows 
Arrhenius' equation, the AADP of the \eqn{i}th year is equal to  

\deqn{\mathrm{AADP}_{i} = \sum_{j=S}^{E_{i}}\mathrm{exp}\left(B - \frac{E_{a}}{R\,T_{ij}}\right),}

where \eqn{T_{ij}} represents the mean daily temperature of the 
\eqn{j}th day of the \eqn{i}th year (in K). In theory, \eqn{\mathrm{AADP}_{i} = 100\%}, 
i.e., the AADP values of different years are a constant of 100\%. However, in practice, there is 
a certain deviation of \eqn{\mathrm{AADP}_{i}} from 100\%. The following approach 
is used to determine the predicted occurrence time. 
When \eqn{\sum_{j=S}^{F}r_{ij} = 100\%} (where \eqn{F \geq S}), it follows that \eqn{F} is 
the predicted occurrence time; when \eqn{\sum_{j=S}^{F}r_{ij} < 100\%} and 
\eqn{\sum_{j=S}^{F+1}r_{ij} > 100\%}, the trapezoid method (Ring and Harris, 1983) 
is used to determine the predicted occurrence time. 

\eqn{\qquad}The argument of \code{expr} can be any an arbitrary user-defined temperature-dependent 
developmental rate function, e.g., a function named \code{myfun}, 
but it needs to take the following form of \code{myfun <- function(P, x){...}}, 
where \code{P} is the vector of the model parameter(s), and \code{x} is the vector of the 
predictor variable, i.e., the temperature variable.  
}
\value{
  \item{Year}{the years with climate data}
  \item{Time.pred}{the predicted occurence times (day of year) in different years}
}
\references{
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
There is a need to note that the unit of \code{Temp} in \bold{Arguments} is \eqn{{}^{\circ}}C, not K. 
In addition, when using Arrhenius' equation to describe \eqn{r}, to reduce the size of \eqn{B} 
in this equation, Arrhenius' equation is multiplied by \eqn{{10}^{12}} in calculating the 
AADP value for each year, i.e., 

\deqn{\mathrm{AADP}_{i} = \sum_{j=S}^{E_{i}}\left[{10}^{12} \cdot \mathrm{exp}\left(B - \frac{E_{a}}{R\,T_{ij}}\right)\right].}
}
\seealso{
\code{\link{ADP}}
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
S.val      <- 47

# Defines a re-parameterized Arrhenius' equation
Arrhenius.eqn <- function(P, x){
  B  <- P[1]
  Ea <- P[2]
  R  <- 1.987 * 10^(-3)
  x  <- x + 273.15
  10^12*exp(B-Ea/(R*x))
}

P0 <- c(-4.7823, 14.8198)
T2 <- seq(-10, 20, len = 2000)
r2 <- Arrhenius.eqn(P = P0, x = T2)

dev.new()
par1 <- par(family="serif")
par2 <- par(mar=c(5, 5, 2, 2))
par3 <- par(mgp=c(3, 1, 0))
plot( T2, r2, cex.lab = 1.5, cex.axis = 1.5, pch = 1, cex = 1.5, col = 2, type = "l", 
      xlab = expression(paste("Temperature (", degree, "C)", sep = "")), 
      ylab = expression(paste("Developmental rate (", {day}^{"-1"}, ")", sep="")) ) 
par(par1)
par(par2)
par(par3)

res6 <- predADP( S = S.val, expr = Arrhenius.eqn, theta = P0, Year2 = Year2.val, 
                 DOY = DOY.val, Temp = Temp.val, DOY.ul = DOY.ul.val )
res6

ind5 <- res6$Year \%in\% intersect(res6$Year, Year1.val)
ind6  <- Year1.val \%in\% intersect(res6$Year, Year1.val)
RMSE3 <- sqrt( sum((Time.val[ind6]-res6$Time.pred[ind5])^2) / length(Time.val[ind6]) ) 
RMSE3 

}