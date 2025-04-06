\name{predADP2}
\alias{predADP2}
\title{
Prediction Function of the Accumulated Developmental Progress Method Using Minimum and Maximum Daily Temperatures
}
\description{
Predicts the occurrence times using the accumulated developmental progress (ADP) 
method based on observed or predicted minimum and maximum daily air temperatures 
(Wagner et al., 1984; Shi et al., 2017a, b).
}
\usage{
predADP2(S, expr, theta, Year2, DOY, Tmin, Tmax, DOY.ul = 120)
}
\arguments{
  \item{S}{the starting date for thermal accumulation (in day-of-year)}
  \item{expr}{a user-defined model that is used in the accumulated developmental progress (ADP) method}
  \item{theta}{a vector saves the numerical values of the parameters in \code{expr}}
  \item{Year2}{the vector of the years recording the climate data for predicting the occurrence times}
  \item{DOY}{the vector of the dates (in day-of-year) for which climate data exist}
  \item{Tmin}{the minimum daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{Tmax}{the maximum daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{DOY.ul}{the upper limit of \code{DOY} used to predict the occurrence time}
}
\details{
Organisms exhibiting phenological events in early spring often experience several cold days 
during their development. In this case, Arrhenius' equation (Shi et al., 2017a, b, 
and references therein) has been recommended to describe the effect of the absolute temperature 
(\eqn{T} in Kelvin [K]) on the developmental rate (\eqn{r}):
\deqn{r = \mathrm{exp}\left(B - \frac{E_{a}}{R\,T}\right),}
where \eqn{E_{a}} represents the activation free energy (in kcal \eqn{\cdot} mol\eqn{{}^{-1}}); 
\eqn{R} is the universal gas constant (= 1.987 cal \eqn{\cdot} mol\eqn{{}^{-1}} \eqn{\cdot} K\eqn{{}^{-1}}); 
\eqn{B} is a constant. To maintain consistency between the units used for \eqn{E_{a}} and \eqn{R}, we need to 
re-assign \eqn{R} to be 1.987\eqn{\times {10}^{-3}}, making its unit 1.987\eqn{\times {10}^{-3}} 
kcal \eqn{\cdot} mol\eqn{{}^{-1}} \eqn{\cdot} K\eqn{{}^{-1}} in the above formula. 

\eqn{\qquad}In the accumulated developmental progress (ADP) method, when the annual accumulated developmental 
progress (AADP) reaches 100\%, the phenological event is predicted to occur for each year. 
Let \eqn{\mathrm{AADP}_{i}} denote the AADP of the \eqn{i}th year, which equals

\deqn{\mathrm{AADP}_{i} = \sum_{j=S}^{E_{i}}\sum_{w=1}^{24}\frac{r_{ijw}}{24},}

where \eqn{E_{i}} represents the ending date (in day-of-year), i.e., the occurrence time of a pariticular 
phenological event in the \eqn{i}th year. \eqn{r_{ijw}} is the developmental rate (per hour), which is 
transferred to \eqn{r_{ij}} (per day) by dividing 24. If the temperature-dependent developmental rate follows 
Arrhenius' equation, the AADP of the \eqn{i}th year is equal to  

\deqn{\mathrm{AADP}_{i} = \sum_{j=S}^{E_{i}}\sum_{w=1}^{24}\left\{\frac{1}{24}\mathrm{exp}\left(B - \frac{E_{a}}{R\,T_{ijw}}\right)\right\},}

where \eqn{T_{ijw}} represents the estimated mean hourly temperature of the \eqn{w}th hour of the  
\eqn{j}th day of the \eqn{i}th year (in K). This packages takes the method proposed by 
Zohner et al. (2020) to estimate the mean hourly temperature (\eqn{T_{w}}) for each of 24 hours:

\deqn{T_{w} = \frac{T_{\mathrm{max}} - T_{\mathrm{min}}}{2}\, \mathrm{sin}\left(\frac{w\pi}{12}-
\frac{\pi}{2}\right)+\frac{T_{\mathrm{max}} + T_{\mathrm{min}}}{2},}

where \eqn{w} represents the \eqn{w}th hour of a day, and \eqn{T_{\mathrm{min}}} and \eqn{T_{\mathrm{max}}} 
represent the minimum and maximum temperatures of the day, respectively. 

\eqn{\qquad}In theory, \eqn{\mathrm{AADP}_{i} = 100\%}, 
i.e., the AADP values of different years are a constant 100\%. However, in practice, there is 
a certain deviation of \eqn{\mathrm{AADP}_{i}} from 100\%. The following approach 
is used to determine the predicted occurrence time. 
When \eqn{\sum_{j=S}^{F}\sum_{w=1}^{24}\left(r_{ijw}/24\right) = 100\%} (where \eqn{F \geq S}), it follows that \eqn{F} is 
the predicted occurrence time; when \eqn{\sum_{j=S}^{F}\sum_{w=1}^{24}\left(r_{ijw}/24\right) < 100\%} and 
\eqn{\sum_{j=S}^{F+1}\sum_{w=1}^{24}\left(r_{ijw}/24\right) > 100\%}, the trapezoid method (Ring and Harris, 1983) 
is used to determine the predicted occurrence time. 

\eqn{\qquad}The argument of \code{expr} can be any an arbitrary user-defined temperature-dependent 
developmental rate function, e.g., a function named \code{myfun}, 
but it needs to take the form of \code{myfun <- function(P, x){...}}, 
where \code{P} is the vector of the model parameter(s), and \code{x} is the vector of the 
predictor variable, i.e., the temperature variable.  
}
\value{
  \item{Year}{the years with climate data}
  \item{Time.pred}{the predicted occurrence times (day-of-year) in different years}
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

Zohner, C.M., Mo, L., Sebald, V., Renner, S.S. (2020) Leaf-out in northern ecotypes of wide-ranging 
trees requires less spring warming, enhancing the risk of spring frost damage at cold limits. 
\emph{Global Ecology and Biogeography} 29, 1056\eqn{-}1072. \doi{10.1111/geb.13088}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Zhenghong Chen \email{chenzh64@126.com},
Jing Tan \email{jmjwyb@163.com}, Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\note{
The entire minimum and maximum daily temperature data set for the spring of each year should be provided. 
It should be noted that the unit of \code{Tmin} and \code{Tmax} in \bold{Arguments} is \eqn{{}^{\circ}}C, not K. 
In addition, when using Arrhenius' equation to describe \eqn{r}, to reduce the size of \eqn{B} 
in this equation, Arrhenius' equation is multiplied by \eqn{{10}^{12}} in calculating the 
AADP value for each year, i.e., 

\deqn{\mathrm{AADP}_{i} = \sum_{j=S}^{E_{i}}\sum_{w=1}^{24}\left[{10}^{12} \cdot \frac{1}{24} \cdot \mathrm{exp}\left(B - \frac{E_{a}}{R\,T_{ijw}}\right)\right].}
}
\seealso{
\code{\link{ADP2}}
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
Tmin.val   <- X2$MinDT
Tmax.val   <- X2$MaxDT
DOY.ul.val <- 120
S.val      <- 46

# Defines a re-parameterized Arrhenius' equation
Arrhenius.eqn <- function(P, x){
  B  <- P[1]
  Ea <- P[2]
  R  <- 1.987 * 10^(-3)
  x  <- x + 273.15
  10^12*exp(B-Ea/(R*x))
}

P0 <- c(8.220327, 22.185942)
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

\donttest{
  cand.res6 <- predADP2( S = S.val, expr = Arrhenius.eqn, theta = P0, Year2 = Year2.val, 
                         DOY = DOY.val, Tmin = Tmin.val, Tmax = Tmax.val, DOY.ul = DOY.ul.val )
  cand.res6

  ind5  <- cand.res6$Year \%in\% intersect(cand.res6$Year, Year1.val)
  ind6  <- Year1.val \%in\% intersect(cand.res6$Year, Year1.val)
  RMSE3 <- sqrt( sum((Time.val[ind6]-cand.res6$Time.pred[ind5])^2) / length(Time.val[ind6]) ) 
  RMSE3 
}
}