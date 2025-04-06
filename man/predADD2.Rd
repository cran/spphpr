\name{predADD2}
\alias{predADD2}
\title{
Prediction Function of the Accumulated Degree Days Method Using Minimum and Maximum Daily Temperatures
}
\description{
Predicts the occurrence times using the accumulated degree days method based on 
observed or predicted minimum and maximum daily air temperatures (Aono, 1993; Shi et al., 2017a, b).
}
\usage{
predADD2(S, T0, AADD, Year2, DOY, Tmin, Tmax, DOY.ul = 120)
}
\arguments{
  \item{S}{the starting date for thermal accumulation (in day-of-year)}
  \item{T0}{the base temperature (in \eqn{{}^{\circ}}C)}
  \item{AADD}{the expected annual accumulated degree days}
  \item{Year2}{the vector of the years recording the climate data for predicting the occurrence times}
  \item{DOY}{the vector of the dates (in day-of-year) for which climate data exist}
  \item{Tmin}{the minimum daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{Tmax}{the maximum daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{DOY.ul}{the upper limit of \code{DOY} used to predict the occurrence time}
}
\details{
In the accumulated degree days (ADD) method (Shi et al., 2017a, b), the starting date 
(\eqn{S}), the base temperature 
(\eqn{T_{0}}), and the annual accumulated degree days (AADD, which is denoted by \eqn{k}) 
are assumed to be constants across different years. Let \eqn{k_{i}} denote the AADD of 
the \eqn{i}th year, which equals

\deqn{k_{i} = \sum_{j=S}^{E_{i}}\sum_{w=1}^{24}\frac{\left(T_{ijw}-T_{0}\right)}{24},}

where \eqn{E_{i}} represents the ending date (in day-of-year), i.e., the occurrence time of a particular 
phenological event in the \eqn{i}th year, and \eqn{T_{ijw}} represents the estimated mean hourly temperature of 
the \eqn{w}th hour of the \eqn{j}th day of the \eqn{i}th year (in \eqn{{}^{\circ}}C). 
If \eqn{T_{ijw} \le T_{0}}, \eqn{T_{ijw} - T_{0}} is defined to be zero. This packages takes the method proposed by 
Zohner et al. (2020) to estimate the mean hourly temperature (\eqn{T_{w}}) for each of 24 hours:

\deqn{T_{w} = \frac{T_{\mathrm{max}} - T_{\mathrm{min}}}{2}\, \mathrm{sin}\left(\frac{w\pi}{12}-
\frac{\pi}{2}\right)+\frac{T_{\mathrm{max}} + T_{\mathrm{min}}}{2},}

where \eqn{w} represents the \eqn{w}th hour of a day, and \eqn{T_{\mathrm{min}}} and \eqn{T_{\mathrm{max}}} 
represent the minimum and maximum temperatures of the day, respectively. 

\eqn{\qquad}In theory, \eqn{k_{i} = k}, i.e., the AADD values of different years are a constant. 
However, in practice, there is a certain deviation of \eqn{k_{i}} from \eqn{k} that is estimated by \eqn{\overline{k}}  
(i.e., the mean of the \eqn{k_{i}} values). The following approach is used to determine the predicted occurrence time. 
When \eqn{\sum_{j=S}^{F}\sum_{w=1}^{24}\left(T_{ijw}-T_{0}\right)/24 = \overline{k}} (where \eqn{F \geq S}), 
it follows that \eqn{F} is 
the predicted occurrence time; when \eqn{\sum_{j=S}^{F}\sum_{w=1}^{24}\left(T_{ijw}-T_{0}\right)/24 < \overline{k}} and 
\eqn{\sum_{j=S}^{F+1}\sum_{w=1}^{24}\left(T_{ijw}-T_{0}\right)/24 > \overline{k}}, the trapezoid method (Ring and Harris, 1983) 
is used to determine the predicted occurrence time. 

}
\value{
  \item{Year}{the years with climate data}
  \item{Time.pred}{the predicted occurrence times (day-of-year) in different years}
}
\references{
Aono, Y. (1993) Climatological studies on blooming of cherry tree (\emph{Prunus yedoensis}) by means 
of DTS method. \emph{Bulletin of the University of Osaka Prefecture. Ser. B, Agriculture and life sciences} 
45, 155\eqn{-}192 (in Japanese with English abstract).

Ring, D.R., Harris, M.K. (1983) Predicting pecan nut casebearer (Lepidoptera: Pyralidae) activity 
at College Station, Texas. \emph{Environmental Entomology} 12, 482\eqn{-}486. \doi{10.1093/ee/12.2.482}

Shi, P., Chen, Z., Reddy, G.V.P., Hui, C., Huang, J., Xiao, M. (2017a) Timing of cherry tree blooming: 
Contrasting effects of rising winter low temperatures and early spring temperatures. 
\emph{Agricultural and Forest Meteorology} 240\eqn{-}241, 78\eqn{-}89. \doi{10.1016/j.agrformet.2017.04.001}


Shi, P., Fan, M., Reddy, G.V.P. (2017b) Comparison of thermal performance equations in describing 
temperature-dependent developmental rates of insects: (III) Phenological applications.  
\emph{Annals of the Entomological Society of America} 110, 558\eqn{-}564. \doi{10.1093/aesa/sax063}

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
}
\seealso{
\code{\link{ADD2}}
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
S.val      <- 65
T0.val     <- 3.8
AADD.val   <- 136.5805

\donttest{
  cand.res2 <- predADD2( S = S.val, T0 = T0.val, AADD = AADD.val, Year2 = Year2.val, 
                         DOY = DOY.val, Tmin = Tmin.val, Tmax = Tmax.val, 
                         DOY.ul = DOY.ul.val )
  cand.res2


  ind1  <- cand.res2$Year \%in\% intersect(cand.res2$Year, Year1.val)
  ind2  <- Year1.val \%in\% intersect(cand.res2$Year, Year1.val)
  RMSE1 <- sqrt( sum((Time.val[ind2]-cand.res2$Time.pred[ind1])^2) / length(Time.val[ind2]) ) 
  RMSE1 
}

}