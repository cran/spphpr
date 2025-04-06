\name{ADP2}
\alias{ADP2}
\title{
Function for Implementing the Accumulated Developmental Progress Method Using Minimum and Maximum Daily Temperatures
}
\description{
Estimates the starting date (\eqn{S}, in day-of-year) and the parameters of 
a developmental rate model in the accumulated developmental progress (ADP) 
method using minimum and maximum daily air temperatures (Wagner et al., 1984; Shi et al., 2017a, b).
}
\usage{
ADP2( S.arr, expr, ini.val, Year1, Time, Year2, DOY, Tmin, Tmax, 
      DOY.ul = 120, fig.opt = TRUE, control = list(), verbose = TRUE )
}
\arguments{
  \item{S.arr}{the candidate starting dates for thermal accumulation (in day-of-year)}
  \item{expr}{a user-defined model that is used in the accumulated developmental progress (ADP) method}
  \item{ini.val}{a vector or a list that saves the initial values of the parameters in \code{expr}}
  \item{Year1}{the vector of the years in which a particular phenological event was recorded}
  \item{Time}{the vector of the occurrence times (in day-of-year) of a particular phenological event across many years}
  \item{Year2}{the vector of the years recording the climate data corresponding to the occurrence times}
  \item{DOY}{the vector of the dates (in day-of-year) for which climate data exist}
  \item{Tmin}{the minimum daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{Tmax}{the maximum daily air temperature data (in \eqn{{}^{\circ}}C) corresponding to \code{DOY}}
  \item{DOY.ul}{the upper limit of \code{DOY} used to predict the occurrence time}
  \item{fig.opt}{an optional argument to draw the figures associated with the temperature-dependent developmental rate curve, 
the mean daily temperatures versus years, and a comparison between the predicted and observed occurrence times}
  \item{control}{the list of control parameters for using the \code{\link{optim}} function in the \pkg{stats} package}
  \item{verbose}{an optional argument allowing users to suppress the printing of computation progress}
}
\details{

It is better not to set too many candidate starting dates, as doing so will be time-consuming. If \code{expr} 
is selected as Arrhenius' equation, \code{S.arr} can be selected as the \code{S} obtained from the output of 
the \code{\link{ADTS2}} function. Here, \code{expr} can be other nonlinear temperature-dependent 
developmental rate functions (see Shi et al. [2017b] for details). Further, \code{expr} can be any an arbitrary 
user-defined temperature-dependent developmental rate function, e.g., a function named \code{myfun}, 
but it needs to take the form of \code{myfun <- function(P, x){...}}, 
where \code{P} is the vector of the model parameter(s), and \code{x} is the vector of the 
predictor variable, i.e., the temperature variable.

\eqn{\qquad}The function does not require that \code{Year1} is the same as \code{unique(Year2)}, 
and the intersection of the two vectors of years will be kept. The unused years that have phenological 
records but lack climate data will be showed in \code{unused.years} in the returned list. 

\eqn{\qquad}The numerical value of \code{DOY.ul} should be greater than or equal to the maximum \code{Time}.

\eqn{\qquad} Let \eqn{r} represent the temperature-dependent developmental rate, i.e., 
the reciprocal of the developmental 
duration required for completing a particular phenological event, at a constant temperature. 
In the accumulated developmental progress (ADP) method, when the annual accumulated developmental 
progress (AADP) reaches 100\%, the phenological event is predicted to occurr for each year. 
Let \eqn{\mathrm{AADP}_{i}} denote the AADP of the \eqn{i}th year, which equals

\deqn{\mathrm{AADP}_{i} = \sum_{j=S}^{E_{i}}\sum_{w=1}^{24}\frac{r_{ijw}\left(\mathrm{\mathbf{P}}; T_{ijw}\right)}{24},}

where \eqn{S} represents the starting date (in day-of-year), \eqn{E_{i}} represents the ending date 
(in day-of-year), i.e., the occurrence time of a pariticular phenological event in the \eqn{i}th year, 
\eqn{\mathrm{\mathbf{P}}} is the vector of the model parameters in \code{expr}, 
\eqn{T_{ijw}} represents the estimated mean hourly temperature of the \eqn{w}th hour of the \eqn{j}th day of the \eqn{i}th 
year (in \eqn{{}^{\circ}}C or K), and \eqn{r_{ijw}} represents the developmental rate (per hour) at \eqn{T_{ijw}}, 
which is transferred to \eqn{r_{ij}} (per day) by dividing 24. This packages takes the method proposed 
by Zohner et al. (2020) to estimate the mean hourly temperature for each of 24 hours:

\deqn{T_{w} = \frac{T_{\mathrm{max}} - T_{\mathrm{min}}}{2}\, \mathrm{sin}\left(\frac{w\pi}{12}-
\frac{\pi}{2}\right)+\frac{T_{\mathrm{max}} + T_{\mathrm{min}}}{2},}

where \eqn{w} represents the \eqn{w}th hour of a day, and \eqn{T_{\mathrm{min}}} and \eqn{T_{\mathrm{max}}} 
represent the minimum and maximum temperatures of the day, respectively. 
 
\eqn{\qquad}In theory, \eqn{\mathrm{AADP}_{i} = 100\%}, 
i.e., the AADP values of different years are a constant 100\%. However, in practice, there is 
a certain deviation of \eqn{\mathrm{AADP}_{i}} from 100\%. The following approach 
is used to determine the predicted occurrence time. 
When \eqn{\sum_{j=S}^{F}\sum_{w=1}^{24}r_{ijw}/24 = 100\%} (where \eqn{F \geq S}), it follows that \eqn{F} is 
the predicted occurrence time; when \eqn{\sum_{j=S}^{F}\sum_{w=1}^{24}r_{ijw}/24 < 100\%} and 
\eqn{\sum_{j=S}^{F+1}\sum_{w=1}^{24}r_{ijw}/24 > 100\%}, the trapezoid method (Ring and Harris, 1983) 
is used to determine the predicted occurrence time. Let \eqn{\hat{E}_{i}} represent the predicted 
occurrence time of the \eqn{i}th year. Assume that there are \eqn{n}-year phenological records. 
When the starting date \eqn{S} and the temperature-dependent developmental rate model are known, 
the model parameters can be estimated using the Nelder-Mead optiminization method 
(Nelder and Mead, 1965) to minimize the root-mean-square error (RMSE) between the observed and predicted 
occurrence times, i.e.,

\deqn{\mathrm{\mathbf{\hat{P}}} = \mathrm{arg}\,\mathop{\mathrm{min}}\limits_{
\mathrm{\mathbf{P}}}\left\{\mathrm{RMSE}\right\} = \mathrm{arg}\,\mathop{\mathrm{min}}\limits_{
\mathrm{\mathbf{P}}}\sqrt{\frac{\sum_{i=1}^{n}\left(E_{i}-\hat{E}_i\right){}^{2}}{n}}.}

Because \eqn{S} is not determined, a group of candidate \eqn{S} values (in day-of-year) need to 
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
day-of-year, estimated mean daily temperature (= (\code{Tmin} + \code{Tmax})/2) and developmental rate columns}
  \item{MAT}{a matrix consisting of the candidate starting dates and the estimates of 
candidate model parameters with the corresponding RMSEs}
  \item{Dev.accum}{the calculated annual accumulated developmental progresses in different years} 
  \item{Year}{The overlapping years between \code{Year1} and \code{Year2}}
  \item{Time}{The observed occurrence times (day-of-year) in the overlapping years 
between \code{Year1} and \code{Year2}}
  \item{Time.pred}{the predicted occurrence times in different years}
  \item{S}{the determined starting date (day-of-year)}
  \item{par}{the estimates of model parameters}
  \item{RMSE}{the RMSE (in days) between the observed and predicted occurrence times}
  \item{unused.years}{the years that have phenological records but lack climate data}
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
In \code{TDDR}, the first column of \code{Year} saves the years, the second column 
of \code{DOY} saves the day-of-year values, the third column of \code{Temperature} 
saves the estimated mean daily air temperatures (= (\code{Tmin} + \code{Tmax})/2) 
from the starting date to the occurrence times, 
and the fourth column of \code{Rate} saves the calculated developmental rates 
corresponding to the estimated mean daily temperatures. 
}
\seealso{
\code{\link{predADP2}}
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
S.arr0     <- 46

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

# The usage is similar to that of the "ADP" function. There is only a need to 
#   replace "Temp = Temp.val" with "Tmin = Tmin.val, Tmax = Tmax.val" when using 
#   the "ADP2" function.
}