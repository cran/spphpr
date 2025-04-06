
ADP2 <- function( S.arr, expr, ini.val, Year1, Time, Year2, DOY, Tmin, Tmax,
                  DOY.ul=120, fig.opt = TRUE, control = list(), verbose = TRUE ){

  S.arr <- round( S.arr )   
  Time  <- round( Time )
  if( !setequal( unique(Year1), unique(Year2) ) ){
    # Removes the case that 'Year1' isn't in accord with 'Year2'
    Year3        <- intersect(Year1, Year2)
    unused.years <- unique(Year1)[ !(unique(Year1) %in% Year3) ]  
    mat1         <- cbind(Year1, Time)
    mat2         <- cbind(Year2, DOY, Tmin, Tmax)
    new.mat1     <- mat1[Year1 %in% Year3, ]
    new.mat2     <- mat2[Year2 %in% Year3, ]
    Year1        <- new.mat1[,1]
    Time         <- new.mat1[,2]
    Year2        <- new.mat2[,1]
    DOY          <- new.mat2[,2]
    Tmin         <- new.mat2[,3]
    Tmax         <- new.mat2[,4]
  }
  else{
    unused.years <- NULL
  }   
  if( max(S.arr)[1] >= min(Time)[1] ){
    stop("The starting date should be less than the minimal occurrence time.")
  }
  if( min(S.arr)[1] <= 0 ){
    stop("The starting date should be more than 0")
  }
  if(DOY.ul <= max(Time)[1]){
    stop("'DOY.ul' should exceed the maximal occurrence time")
  }
  DOYmax <- as.numeric( tapply(DOY, Year2, max) )
  DOYmin <- as.numeric( tapply(DOY, Year2, min) )
  if( (min(S.arr)[1] < min(DOYmin)[1]) | (max(S.arr)[1] > min(DOYmax)[1]) ){
    stop("The candidates of starting date should be in the range of 'DOY'")
  }
  if( (DOY.ul < min(DOYmin)[1]) | (DOY.ul > min(DOYmax)[1]) ){
    stop("'DOY.ul' should be in the range of 'DOY'")
  }

  expr.accum <- function(P, x){

    if(length(x)==2)    x <- matrix(x, nrow=1, ncol=2)
    T1     <- x[,1]
    T2     <- x[,2]

    temp <- 0
    for(k in 1:nrow(x)){
      Thour <- (T2[k]-T1[k])/2*sin(pi/12*1:24-pi/2)+(T2[k]+T1[k])/2
      Thour <- Thour
      for(q in 1:24){
        temp <- temp + expr(P, Thour[q])/24
      }
    } 
    return(temp)
  }

  ini.val <- as.list(ini.val)
  p       <- length(ini.val)
  r       <- 1
  for (i in 1:p) {
      r <- r * length(ini.val[[i]])
  }
  ini.val <- expand.grid(ini.val)
  len     <- length(S.arr) * r
  #### S (as the first element) and RMSE (as the last element) ####
  ####     will be added to every row of this matrix ##############
  MAT     <- matrix(NA, nrow = len, ncol = (p + 2))
  #################################################################
  q   <- 0
  for(w in 1:length(S.arr)){

    myfun <- function(P){ 
      Time.pred <- c()  
      for(i in 1:length(Year1)){
        temp.Tmin    <- Tmin[ Year2 == Year1[i] ]
        temp.Tmax    <- Tmax[ Year2 == Year1[i] ]
        temp.Tcomb   <- cbind(temp.Tmin, temp.Tmax)
        temp.DOY     <- DOY[ Year2 == Year1[i] ]     
        for(k in which(temp.DOY == S.arr[w]): which(temp.DOY == DOY.ul)){
          temp.progress <- expr.accum(P, temp.Tcomb[which(temp.DOY == S.arr[w]):k, ])
          if(temp.progress >= 1){
            # y1: The upper base length of the trapezoid
            # y2: The length of the line segment parallel to the two bases of the trapezoid 
            # y3: The lower base length of the trapezoid
            # x3: The height of the trapezoid is temp.DOY[k] - temp.DOY[k-1]     
            y1 <- expr.accum(P, temp.Tcomb[which(temp.DOY == S.arr[w]): (k-1), ])
            y2 <- 1
            y3 <- temp.progress   
            x3 <- temp.DOY[k] - temp.DOY[k-1]     
            if(k > 1){     
              Time.pred[i] <- temp.DOY[k-1] + (y2-y1)/(y3-y1) * x3 
            }         
            if(k == 1){
              Time.pred[i] <- 1
            } 
            break
          }
          extreme.temp     <- temp.Tcomb[which(temp.DOY == S.arr[w]): which(temp.DOY == DOY.ul), ]
          extreme.progress <- expr.accum(P, extreme.temp)
          if(extreme.progress < 1){
            Time.pred[i] <- DOY.ul
          } 
        }      
      }
      RMSE <- sqrt(sum((Time.pred-Time)^2)/length(Time))  
      return(RMSE)
    } 

    for (i in 1:nrow(ini.val)) {
      q <- q + 1 
      if(verbose){
        cat("\n")
        print(paste("The current computation progress is ", q, "/", len, sep=""))
        cat("\n")
      }
      Res      <- optim( ini.val[i, ], myfun, control=control )        
      MAT[q, ] <- c(S.arr[w], Res$par, Res$value)  
    }

  }

  M <- MAT[ MAT[,p + 2] == min(MAT[, p + 2]) ]
  M <- matrix(M, ncol=ncol(MAT))
  if (nrow(M)>1){
    index1 <- sample(1:nrow(M), replace=TRUE)[1]
    M      <- M[index1,]
  }

  Dev.accum <- c()
  for(i in 1:length(Year1)){
    temp.Tmin    <- Tmin[ Year2 == Year1[i] ]
    temp.Tmax    <- Tmax[ Year2 == Year1[i] ]
    temp.Tcomb   <- cbind(temp.Tmin, temp.Tmax)
    temp.DOY     <- DOY[ Year2 == Year1[i] ]
    Dev.accum[i] <- expr.accum(M[2:(p+1)], temp.Tcomb[which(temp.DOY == M[1]):which(temp.DOY == Time[i]), ])
  }
  Time.pred <- c()  
  for(i in 1:length(Year1)){
    temp.Tmin    <- Tmin[ Year2 == Year1[i] ]
    temp.Tmax    <- Tmax[ Year2 == Year1[i] ]
    temp.Tcomb   <- cbind(temp.Tmin, temp.Tmax)  
    temp.DOY     <- DOY[ Year2 == Year1[i] ] 
    for(k in which(temp.DOY == M[1]): which(temp.DOY == DOY.ul)){
      temp.progress <- expr.accum(M[2:(p+1)], temp.Tcomb[which(temp.DOY == M[1]): k, ])
      if(temp.progress >= 1){
        # y1: The upper base length of the trapezoid
        # y2: The length of the line segment parallel to the two bases of the trapezoid 
        # y3: The lower base length of the trapezoid
        # x3: The height of the trapezoid is temp.DOY[k] - temp.DOY[k-1]     
        y1           <- expr.accum(M[2:(p+1)], temp.Tcomb[which(temp.DOY == M[1]): (k-1), ])
        y2           <- 1
        y3           <- temp.progress   
        x3           <- temp.DOY[k] - temp.DOY[k-1]     
        if(k > 1){     
          Time.pred[i] <- temp.DOY[k-1] + (y2-y1)/(y3-y1) * x3 
        }      
        if(k == 1){
          Time.pred[i] <- 1
        }  
        break
      }
      extreme.temp     <- temp.Tcomb[which(temp.DOY == M[1]): which(temp.DOY == DOY.ul), ]
      extreme.progress <- expr.accum(M[2:(p+1)], extreme.temp)
      if(extreme.progress < 1){
        Time.pred[i] <- DOY.ul
      } 
    }      
  }
  Tmin.arr <- c()
  Tmax.arr <- c()
  DOY.arr  <- c()
  Year.arr <- c()
  for(i in 1:length(Year1)){
    temp.Tmin  <- Tmin[ Year2 == Year1[i] ]
    temp.Tmax  <- Tmax[ Year2 == Year1[i] ]
    temp.DOY   <- DOY[ Year2 == Year1[i] ]
    var11      <- temp.Tmin[which(temp.DOY == M[1]): which(temp.DOY == Time[i])] 
    var12      <- temp.Tmax[which(temp.DOY == M[1]): which(temp.DOY == Time[i])] 
    var2       <- temp.DOY[which(temp.DOY == M[1]): which(temp.DOY == Time[i])] 
    var3       <- rep(Year1[i], len=length(var11))    
    Tmin.arr   <- c(Tmin.arr, var11)
    Tmax.arr   <- c(Tmax.arr, var12) 
    DOY.arr    <- c(DOY.arr, var2)
    Year.arr   <- c(Year.arr, var3)     
  }

  Tmean.arr <- (Tmin.arr + Tmax.arr)/2
  Rate.arr  <- expr( M[2:(p+1)], Tmean.arr )

  if( fig.opt ){

    dev.new()
    par1 <- par(family="serif")
    par2 <- par(mar=c(5, 5, 2, 2))
    par3 <- par(mgp=c(3, 1, 0))
    on.exit(par(par1))
    on.exit(par(par2))
    on.exit(par(par3))
    plot(Year1, Dev.accum*100, xlab="Year", 
         ylab="Accumulated developmental progress (%)", 
         ylim=c(50, 150), cex.lab=1.5, cex.axis=1.5, type="n")
    abline(h=1*100, lwd=1, col=4, lty=2) 
    points(Year1, Dev.accum*100, pch=1, cex=1.5, col=2)

    dev.new()
    par1 <- par(family="serif")
    par2 <- par(mar=c(5, 5, 2, 2))
    par3 <- par(mgp=c(3, 1, 0))
    on.exit(par(par1))
    on.exit(par(par2))
    on.exit(par(par3))
    plot( Tmean.arr, Rate.arr, cex.lab=1.5, cex.axis=1.5, pch=1, cex=1.5, col=2,
          xlab=expression(paste("Mean daily temperature (", degree, "C)", sep="")), 
          ylab=expression(paste("Calculated developmental rate ( ", {day}^{"-1"}, ")", sep="")) ) 

    dev.new()
    par1 <- par(family="serif")
    par2 <- par(mar=c(5, 5, 2, 2))
    par3 <- par(mgp=c(3, 1, 0))
    on.exit(par(par1))
    on.exit(par(par2))
    on.exit(par(par3))
    ll       <- min(c(Time.pred, Time))[1]
    ul       <- max(c(Time.pred, Time))[1]
    interval <- (ul-ll)/8
    plot( Time, Time.pred, xlim=c(ll-interval, ul+interval), 
          ylim=c(ll-interval, ul+interval), 
          xlab="Observed occurrence time (day-of-year)",
          ylab="Predicted occurrence time (day-of-year)", 
          cex.lab=1.5, cex.axis=1.5, type="n" )
    abline(0, 1, lwd=1, col=4) 
    points(Time, Time.pred, pch=16, cex=1.5, col=2)  

    dev.new()
    par1 <- par(family="serif")
    par2 <- par(mar=c(5, 5, 2, 2))
    par3 <- par(mgp=c(3, 1, 0))
    on.exit(par(par1))
    on.exit(par(par2))
    on.exit(par(par3))
    m.val    <- as.numeric( tapply(Tmean.arr, Year.arr, median) )
    y1.val   <- as.numeric( tapply(Year.arr, Year.arr, mean) )
    ind1     <- sort(m.val, ind=TRUE)$ix
    y2.val   <- y1.val[ind1]
    col1.val <- topo.colors(length(m.val))
    ind2     <- sort(y2.val, ind=TRUE)$ix
    col2.val <- col1.val[ind2]  
    boxplot(Tmean.arr ~ Year.arr, xlab="Year",
         ylab=expression(paste("Mean daily temperature (", degree, "C)", sep="")), 
         cex.lab=1.5, cex.axis=1.5, pch=1, cex=1, col=col2.val) 
  }


  TDDR <- data.frame(Year=Year.arr, DOY=DOY.arr, Temperature=Tmean.arr, Rate=Rate.arr)
  list( TDDR=TDDR, MAT=MAT, Dev.accum = Dev.accum, Year=Year1, Time=Time, Time.pred=Time.pred,    
        S=M[1], par=M[2:(p+1)], RMSE=M[p+2], unused.years=unused.years )
}


