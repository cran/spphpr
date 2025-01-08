\name{predADTS}
\alias{predADTS}
\title{
Prediction Function of the Accumulated Days Transferred 
to a Standardized Temperature Method
}
\description{
Predicts the occurrence times using the accumulated days transferred to a standardized 
temperature (ADTS) method based on observed or predicted mean daily air temperatures 
(Konno and Sugihara, 1986; Aono, 1993; Shi et al., 2017a, 2017b).
}
\usage{
predADTS(S, Ea, AADTS, Year2, DOY, Temp, DOY.ul = 120)
}
\arguments{
  \item{S}{the starting date for thermal accumulation (in day-of-year)}
  \item{Ea}{the activation free energy (in kcal \eqn{\cdot} mol\eqn{{}^{-1}})}
  \item{AADTS}{the expected annual accumulated days transferred to a standardized temperature}
  \item{Year2}{the vector of the years recording the climate data for predicting the occurrence times}
  \item{DOY}{the vector of the dates (in day-of-year) for which climate data exist}
  \item{Temp}{the mean daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{DOY.ul}{the upper limit of \code{DOY} used to predict the occurrence time}
}
\details{
Organisms exhibiting phenological events in early spring often experience several cold days 
during their development. In this case, Arrhenius' equation (Shi et al., 2017a, 2017b, 
and references therein) has been recommended to describe the effect of the absolute temperature 
(\eqn{T} in Kelvin [K]) on the developmental rate (\eqn{r}):
\deqn{r = \mathrm{exp}\left(B - \frac{E_{a}}{R\,T}\right),}
where \eqn{E_{a}} represents the activation free energy (in kcal \eqn{\cdot} mol\eqn{{}^{-1}}); 
\eqn{R} is the universal gas constant (= 1.987 cal \eqn{\cdot} mol\eqn{{}^{-1}} \eqn{\cdot} K\eqn{{}^{-1}}); 
\eqn{B} is a constant. To maintain consistence between the units used for \eqn{E_{a}} and \eqn{R}, we need to 
re-assign \eqn{R} to be 1.987\eqn{\times {10}^{-3}}, making its unit 1.987\eqn{\times {10}^{-3}} 
kcal \eqn{\cdot} mol\eqn{{}^{-1}} \eqn{\cdot} K\eqn{{}^{-1}} in the above formula. 

\eqn{\qquad}According to the definition of the developmental rate (\eqn{r}), 
it is the developmental progress per unit time (e.g., per day, per hour), 
which equals the reciprocal of the developmental duration {\eqn{D}}, i.e., \eqn{r = 1/D}. Let \eqn{T_{s}} 
represent the standard temperature (in K), and \eqn{r_{s}} represent the developmental rate at \eqn{T_{s}}. 
Let \eqn{r_{j}} represent the developmental rate at \eqn{T_{j}}, an arbitrary 
temperature (in K). It is apparent that \eqn{D_{s}r_{s} = D_{j}r_{j} = 1}. It follows that 

\deqn{\frac{D_{s}}{D_{j}} = \frac{r_{j}}{r_{s}} = 
\mathrm{exp}\left[\frac{E_{a}\left(T_{j}-T_{s}\right)}{R\,T_{j}\,T_{s}}\right],}

where \eqn{D_{s}/D_{j}} is referred to as the number of days transferred to a standardized temperature 
(DTS) (Konno and Sugihara, 1986; Aono, 1993). 

\eqn{\qquad}In the accumulated days transferred to a standardized temperature (ADTS) method, 
the annual accumulated days transferred to a standardized temperature (AADTS) is assumed to be a constant. 
Let \eqn{\mathrm{AADTS}_{i}} denote the AADTS of the \eqn{i}th year, which equals

\deqn{\mathrm{AADTS}_{i} = \sum_{j=S}^{E_{i}}\left\{\mathrm{exp}\left[\frac{E_{a}\left(T_{ij}-T_{s}\right)}{R\,T_{ij}\,T_{s}}\right]\right\},}

where \eqn{E_{i}} represents the ending date (in day-of-year), i.e., the occurrence time of a pariticular 
phenological event in the \eqn{i}th year, and \eqn{T_{ij}} represents the mean daily temperature of the 
\eqn{j}th day of the \eqn{i}th year (in K). In theory, \eqn{\mathrm{AADTS}_{i} = \mathrm{AADTS}}, 
i.e., the AADTS values of different years are a constant. However, in practice, there is 
a certain deviation of \eqn{\mathrm{AADTS}_{i}} from \eqn{\mathrm{AADTS}}. The following approach 
is used to determine the predicted occurrence time. 
When \eqn{\sum_{j=S}^{F}\left\{\mathrm{exp}\left[\frac{E_{a}\left(T_{ij}-T_{s}\right)}
{R\,T_{ij}\,T_{s}}\right]\right\} = \mathrm{AADTS}} (where \eqn{F \geq S}), it follows that \eqn{F} is 
the predicted occurrence time; when \eqn{\sum_{j=S}^{F}\left\{\mathrm{exp}\left[
\frac{E_{a}\left(T_{ij}-T_{s}\right)}{R\,T_{ij}\,T_{s}}\right]\right\} < \mathrm{AADTS}} and 
\eqn{\sum_{j=S}^{F+1}\left\{\mathrm{exp}\left[\frac{E_{a}\left(T_{ij}-T_{s}\right)}
{R\,T_{ij}\,T_{s}}\right]\right\} > \mathrm{AADTS}}, the trapezoid method (Ring and Harris, 1983) 
is used to determine the predicted occurrence time. 

}
\value{
  \item{Year}{the years with climate data}
  \item{Time.pred}{the predicted occurence times (day-of-year) in different years}
}
\references{
Aono, Y. (1993) Climatological studies on blooming of cherry tree (\emph{Prunus yedoensis}) by means 
of DTS method. \emph{Bulletin of the University of Osaka Prefecture. Ser. B, Agriculture and life sciences} 
45, 155\eqn{-}192 (in Japanese with English abstract).

Konno, T., Sugihara, S. (1986) Temperature index for characterizing biological activity in soil and 
its application to decomposition of soil organic matter. \emph{Bulletin of National Institute for 
Agro-Environmental Sciences} 1, 51\eqn{-}68 (in Japanese with English abstract).

Ring, D.R., Harris, M.K. (1983) Predicting pecan nut casebearer (Lepidoptera: Pyralidae) activity 
at College Station, Texas. \emph{Environmental Entomology} 12, 482\eqn{-}486. \doi{10.1093/ee/12.2.482}

Shi, P., Chen, Z., Reddy, G.V.P., Hui, C., Huang, J., Xiao, M. (2017a) Timing of cherry tree blooming: 
Contrasting effects of rising winter low temperatures and early spring temperatures. 
\emph{Agricultural and Forest Meteorology} 240\eqn{-}241, 78\eqn{-}89. \doi{10.1016/j.agrformet.2017.04.001}


Shi, P., Fan, M., Reddy, G.V.P. (2017b) Comparison of thermal performance equations in describing 
temperature-dependent developmental rates of insects: (III) Phenological applications.  
\emph{Annals of the Entomological Society of America} 110, 558\eqn{-}564. \doi{10.1093/aesa/sax063}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Zhenghong Chen \email{chenzh64@126.com}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\note{
The entire mean daily temperature data set for the spring of each year should be provided. 
It should be noted that the unit of \code{Temp} in \bold{Arguments} is \eqn{{}^{\circ}}C, not K.
}
\seealso{
\code{\link{ADTS}}
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
S.val      <- 47
Ea.val     <- 15 
AADTS.val  <- 8.5879

res4 <- predADTS( S = S.val, Ea = Ea.val, AADTS = AADTS.val, 
                  Year2 = Year2.val, DOY = DOY.val, Temp = Temp.val, 
                  DOY.ul = DOY.ul.val )
res4

ind3  <- res4$Year \%in\% intersect(res4$Year, Year1.val)
ind4  <- Year1.val \%in\% intersect(res4$Year, Year1.val)
RMSE2 <- sqrt( sum((Time.val[ind4]-res4$Time.pred[ind3])^2) / length(Time.val[ind4]) ) 
RMSE2 

}