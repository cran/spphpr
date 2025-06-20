
ADD4 <- function( S.arr, T0.arr, Year1, Time, Year2, DOY, Tmin, Tmax, 
                  DOY.ul = 120, fig.opt = TRUE, verbose = TRUE ){

  S.arr  <- round(S.arr)
  DOY.ul <- round(DOY.ul)   
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
  if(min(S.arr)[1] <= 0){ 
    stop("The candidate of starting date should be greater than 0!")         
  }
  if( max(S.arr)[1] >= min(Time)[1] ){ 
    stop("The maximal candidate of starting date should be less than the minimal occurrence time")         
  }
  if(DOY.ul < max(Time)[1]){
    stop("'DOY.ul' should exceed the maximal occurrence time!")
  }
  DOYmax <- as.numeric( tapply(DOY, Year2, max) )
  DOYmin <- as.numeric( tapply(DOY, Year2, min) )
  if( (min(S.arr)[1] < min(DOYmin)[1]) | (max(S.arr)[1] > min(DOYmax)[1]) ){
    stop("The range starting date candidate should be in the range of 'DOY'")
  }
  if( (DOY.ul < min(DOYmin)[1]) | (DOY.ul > min(DOYmax)[1]) ){
    stop("'DOY.ul' should be in the range of 'DOY'!")
  }

  AADD.fun <- function(Tb, T){ 
    # Description:
    # AADD.fun is used to calcualte the annual accumulated degree days.
    # Arguments:    
    # Tb: The base temperature
    # T: Vector saving daily min and max air temperatures (in degrees Celsius)
    if(length(T)==2)    T <- matrix(T, nrow=1, ncol=2)
    T1 <- T[,1]
    T2 <- T[,2]
    AADD0 <- 0 
    for(k in 1:nrow(T)){
      Thour <- (T2[k]-T1[k])/2*sin(pi/12*1:24-pi/2)+(T2[k]+T1[k])/2
      DD    <- 0 # Daily accumulated degree days 
      for(q in 1:24){   
        if(Thour[q] >= Tb){
          DD <- DD + (Thour[q] - Tb)/24
        }
        else{
          DD <- DD + 0
        }
      }
      AADD0 <- AADD0 + DD
    }
    AADD0
  }  

  len       <- length(S.arr)*length(T0.arr)
  counter0  <- 0
  counter1  <- 0
  mAADD.mat <- matrix(NA, nrow=length(S.arr), ncol=length(T0.arr)) 
  RMSE.mat   <- mAADD.mat
  # RMSE.mat:  The RMSEs in days for the different combinations of S and T0 candidates
  # mAADD.mat: The mean AADD values for the different combinations of S and T0 candidates
  for(S in S.arr){
    counter1  <- counter1 + 1
    counter2  <- 0
    Time.pred <- c()
    for(T0 in T0.arr){
      counter2 <- counter2 + 1
      counter0 <- counter0 + 1
      if(verbose){
        Sys.sleep(.005)
        cat(counter0, paste(" of ", len, "\r", sep=""))             
        flush.console()               
        if (counter0 %% len == 0) cat("\n")
      }
      # AADD.arr: To store the temporary AADD values in different years
      AADD.arr <- c()
      for(i in 1:length(Year1)){
        temp.Tmin    <- Tmin[ Year2 == Year1[i] ]
        temp.Tmax    <- Tmax[ Year2 == Year1[i] ]
        temp.Tcomb   <- cbind(temp.Tmin, temp.Tmax)
        temp.DOY     <- DOY[ Year2 == Year1[i] ]
        ind2         <- which(temp.DOY == S): which(temp.DOY == Time[i])
        AADD.arr[i] <- AADD.fun(T0, temp.Tcomb[ind2, ])
      }
      mAADD.mat[counter1, counter2] <- mean(AADD.arr) 
      Time.pred <- c()  
      for(i in 1:length(Year1)){
        temp.Tmin    <- Tmin[ Year2 == Year1[i] ]
        temp.Tmax    <- Tmax[ Year2 == Year1[i] ]
        temp.Tcomb   <- cbind(temp.Tmin, temp.Tmax)
        temp.DOY     <- DOY[ Year2 == Year1[i] ]         
        for(k in which(temp.DOY == S): which(temp.DOY == DOY.ul)){  
          ind2       <- which(temp.DOY == S): k    
          temp.AADD <- AADD.fun(T0, temp.Tcomb[ind2, ])  
          if(temp.AADD >= mean(AADD.arr)){
            # y1: The upper base length of the trapezoid
            # y2: The length of the line segment parallel to the two bases of the trapezoid 
            # y3: The lower base length of the trapezoid
            # x3: The height of the trapezoid is temp.DOY[k] - temp.DOY[k-1]     
            y1           <- AADD.fun(T0, temp.Tcomb[which(temp.DOY == S): (k-1), ])
            y2           <- mean(AADD.arr)
            y3           <- temp.AADD   
            x3           <- temp.DOY[k] - temp.DOY[k-1] 
            if(k > 1){     
              Time.pred[i] <- temp.DOY[k-1] + (y2-y1)/(y3-y1) * x3 
            }          
            if(k == 1){
              Time.pred[i] <- 1
            }            
            break
          }
          extreme.temp <- temp.Tcomb[which(temp.DOY == S): which(temp.DOY == DOY.ul), ]
          extreme.AADD <- AADD.fun( T0, extreme.temp )
          if(extreme.AADD < mean(AADD.arr)){
            Time.pred[i] <- DOY.ul
          } 
        }      
      }
      RMSE.mat[counter1, counter2] <- sqrt(sum((Time.pred-Time)^2)/length(Time))
    }
  }
  location <- which(RMSE.mat == min(RMSE.mat[!is.na(RMSE.mat)]), arr.ind=TRUE)
  if(length(location) > 0){
    if (nrow(location) == 1 ){
      temp1 <- RMSE.mat[location]
      temp2 <- S.arr[location[1]]
      temp3 <- T0.arr[location[2]]
    }
    if (nrow(location) > 1 ){
      temp1 <- RMSE.mat[location[1]]
      temp2 <- S.arr[location[1, 1]]
      temp3 <- T0.arr[location[1, 2]]
    }
  }
  if(length(location) == 0){
    temp1    <- NA
    temp2    <- NA
    temp3    <- NA
    location <- NA
  }
  goal.RMSE <- temp1
  goal.S    <- temp2
  goal.T0   <- temp3
  if(!is.na(goal.RMSE)){
    RMSE.range <- range(RMSE.mat)
    # Below 'AADD.arr' and 'Time.pred' is used on the condition 
    #   that S and T0 are finally determined.
    AADD.arr <- c()
    for(i in 1:length(Year1)){
      temp.Tmin    <- Tmin[ Year2 == Year1[i] ]
      temp.Tmax    <- Tmax[ Year2 == Year1[i] ]
      temp.Tcomb   <- cbind(temp.Tmin, temp.Tmax)
      temp.DOY     <- DOY[ Year2 == Year1[i] ]
      ind2         <- which(temp.DOY == goal.S): which(temp.DOY == Time[i])
      AADD.arr[i] <- AADD.fun(goal.T0, temp.Tcomb[ind2, ])
    }
    Time.pred <- c()
    for(i in 1:length(Year1)){
      temp.Tmin    <- Tmin[ Year2 == Year1[i] ]
      temp.Tmax    <- Tmax[ Year2 == Year1[i] ]
      temp.Tcomb   <- cbind(temp.Tmin, temp.Tmax)
      temp.DOY     <- DOY[ Year2 == Year1[i] ]    
      for(k in which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)){
        ind2       <- which(temp.DOY == goal.S): k    
        temp.AADD <- AADD.fun(goal.T0, temp.Tcomb[ind2, ])    
        if(temp.AADD >= mean(AADD.arr)){
          y1           <- AADD.fun(goal.T0, temp.Tcomb[which(temp.DOY == goal.S): (k-1), ] )
          y2           <- mean(AADD.arr)
          y3           <- temp.AADD   
          x3           <- temp.DOY[k] - temp.DOY[k-1]      
          if(k > 1){     
            Time.pred[i] <- temp.DOY[k-1] + (y2-y1)/(y3-y1) * x3 
          }          
          if(k == 1){
            Time.pred[i] <- 1
          }  
          break
        }
        extreme.temp <- temp.Tcomb[which(temp.DOY == goal.S): which(temp.DOY == DOY.ul), ]
        extreme.AADD <- AADD.fun( goal.T0, extreme.temp )
        if(extreme.AADD < mean(AADD.arr)){
          Time.pred[i] <- DOY.ul
        } 
      }      
    }
    if( fig.opt ){

      dev.new()
      par1 <- par(family="serif")
      par2 <- par(mar=c(5, 5, 2, 2))
      par3 <- par(mgp=c(3, 1, 0))
      on.exit(par(par1))
      on.exit(par(par2))
      on.exit(par(par3))
      if( length(S.arr) > 1 ){
        if( length(T0.arr) > 1 ){
          cols <- terrain.colors(200)
          image(S.arr, T0.arr, RMSE.mat, col = cols, axes = TRUE, cex.axis = 1.5, cex.lab = 1.5, 
            xlab="Starting date (day-of-year)", 
            ylab=expression(paste("Base temperature (", 
            degree, "C)", sep="")))
          contour(S.arr, T0.arr, RMSE.mat, levels = round(seq(RMSE.range[1], 
            RMSE.range[2], len = 20), 4), add = TRUE, cex = 1.5, col = "#696969", labcex=1.5)
          points(goal.S, goal.T0, cex=1.5, pch=16, col=2)
        }   
        if( length(T0.arr) == 1 ){
          plot( S.arr, RMSE.mat, cex.axis = 1.5, cex.lab = 1.5, xlim=range(S.arr),
            ylab="RMSE (days)", 
            xlab="Starting date (day-of-year)", pch=16, cex=1.5, col=4)
          points(goal.S, goal.RMSE, cex=1.5, pch=16, col=2)  
        }
      }
      if(length(S.arr) == 1){
        if(length(T0.arr) > 1){
          plot( T0.arr, RMSE.mat, cex.axis = 1.5, cex.lab = 1.5, xlim=range(T0.arr),
            ylab="RMSE (days)", xlab=expression(paste("Base temperature (", 
            degree, "C)", sep="")), type="l", col=4, lwd=1)
          points(goal.T0, goal.RMSE, cex=1.5, pch=16, col=2)  
        }
        if(length(T0.arr) == 1){
          plot( T0.arr, RMSE.mat, cex.axis = 1.5, cex.lab = 1.5, xlim=range(T0.arr),
            ylab="RMSE (days)", xlab=expression(paste("Base temperature (", 
            degree, "C)", sep="")), pch=16, cex=1.5, col=4)
          points(goal.T0, goal.RMSE, cex=1.5, pch=16, col=2)  
        }
      }

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
      plot(Time, Time.pred, xlim=c(ll-interval, ul+interval), ylim=c(ll-interval, ul+interval), 
           xlab="Observed occurrence time (day-of-year)",
           ylab="Predicted occurrence time (day-of-year)", 
           cex.lab=1.5, cex.axis=1.5, pch=16, cex=1.5, col=2)
      abline(0, 1, lwd=1, col=4) 
      points(Time, Time.pred, pch=16, cex=1.5, col=2)  
    }
  }
  if(is.na(goal.RMSE)){
    AADD.arr <- NA    
    Time.pred <- NA
  }
  list( mAADD.mat=mAADD.mat, RMSE.mat=RMSE.mat, AADD.arr=AADD.arr, 
        Year=Year1, Time=Time, Time.pred=Time.pred, S=goal.S, T0=goal.T0, 
        AADD=mean(AADD.arr), RMSE=goal.RMSE, unused.years=unused.years )
}
