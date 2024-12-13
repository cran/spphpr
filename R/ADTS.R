
ADTS <- function( S.arr, Ea.arr, Year1, Time, Year2, DOY, Temp, DOY.ul = 120, 
                  fig.opt = TRUE, verbose = TRUE ){

  S.arr  <- round(S.arr)
  DOY.ul <- round(DOY.ul)   
  if( !setequal( unique(Year1), unique(Year2) ) ){
    # Removes the case that 'Year1' isn't in accord with 'Year2'
    Year3        <- intersect(Year1, Year2)
    unused.years <- unique(Year1)[ !(unique(Year1) %in% Year3) ]  
    mat1         <- cbind(Year1, Time)
    mat2         <- cbind(Year2, DOY, Temp)
    new.mat1     <- mat1[Year1 %in% Year3, ]
    new.mat2     <- mat2[Year2 %in% Year3, ]
    Year1        <- new.mat1[,1]
    Time         <- new.mat1[,2]
    Year2        <- new.mat2[,1]
    DOY          <- new.mat2[,2]
    Temp         <- new.mat2[,3]
  }   
  else{
    unused.years <- NULL
  }
  T1 <- Temp
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
  AADTS.fun <- function(Ea, T){ 
    # Description:
    # AADTS.fun is used to calcualte the annual accumulative DTS value
    # Arguments:
    # Ea: The activation free energy (in kcal/mol)
    # T:  Daily mean air temperature (in degrees Celsius)
    T    <- T + 273.15
    Ts   <- 298.15 
    # Ts: The reference temperature in K     
    R    <- 1.987   
    # R: The universal gas constant in cal/mol/deg
    AADTS0 <- 0 
    for(k in 1:length(T)){      
        AADTS0 <- AADTS0 + exp((Ea*1000) * (T[k]-Ts)/(R*T[k]*Ts))
    }
    AADTS0
  }  
  len       <- length(S.arr)*length(Ea.arr)
  counter0  <- 0
  counter1  <- 0
  mAADTS.mat <- matrix(NA, nrow=length(S.arr), ncol=length(Ea.arr)) 
  RMSE.mat  <- mAADTS.mat
  # RMSE.mat:  The RMSEs in days for the different combinations of S and Ea candidates
  # mAADTS.mat: The mean AADTS values for the different combinations of S and Ea candidates
  for(S in S.arr){
    counter1  <- counter1 + 1
    counter2  <- 0
    Time.pred <- c()
    for(Ea in Ea.arr){
      counter2 <- counter2 + 1
      counter0 <- counter0 + 1
      if(verbose){
        Sys.sleep(.005)
        cat(counter0, paste(" of ", len, "\r", sep=""))             
        flush.console()               
        if (counter0 %% len == 0) cat("\n")
      }
      # AADTS.arr: To store the temporary AADTS values in different years
      AADTS.arr <- c()
      for(i in 1:length(Year1)){
        temp.T1    <- T1[ Year2 == Year1[i] ]
        temp.DOY   <- DOY[ Year2 == Year1[i] ]
        MeanTemp1  <- temp.T1[which(temp.DOY == S): which(temp.DOY == Time[i])]
        AADTS.arr[i] <- AADTS.fun(Ea, MeanTemp1)
      }
      mAADTS.mat[counter1, counter2] <- mean(AADTS.arr) 
      Time.pred <- c()  
      for(i in 1:length(Year1)){
        temp.T1  <- T1[ Year2 == Year1[i] ]
        temp.DOY <- DOY[ Year2 == Year1[i] ]          
        for(k in which(temp.DOY == S): which(temp.DOY == DOY.ul)){      
          temp.AADTS <- AADTS.fun(Ea, temp.T1[which(temp.DOY == S): k])          
          if(temp.AADTS >= mean(AADTS.arr)){
            # y1: The upper base length of the trapezoid
            # y2: The length of the line segment parallel to the two bases of the trapezoid 
            # y3: The lower base length of the trapezoid
            # x3: The height of the trapezoid is temp.DOY[k] - temp.DOY[k-1]     
            y1           <- AADTS.fun(Ea, temp.T1[which(temp.DOY == S): (k-1)])
            y2           <- mean(AADTS.arr)
            y3           <- temp.AADTS   
            x3           <- temp.DOY[k] - temp.DOY[k-1] 
            if(k > 1){     
              Time.pred[i] <- temp.DOY[k-1] + (y2-y1)/(y3-y1) * x3 
            }          
            if(k == 1){
              Time.pred[i] <- 1
            }            
            break
          }
          extreme.temp <- temp.T1[which(temp.DOY == S): which(temp.DOY == DOY.ul)]
          extreme.AADTS <- AADTS.fun( Ea, extreme.temp )
          if(extreme.AADTS < mean(AADTS.arr)){
            Time.pred[i] <- DOY.ul
          } 
        }      
      }
      RMSE.mat[counter1, counter2] <- sqrt(sum((Time.pred-Time)^2)/length(Time))
    }
  }
  location <- which(RMSE.mat == min(RMSE.mat[!is.na(RMSE.mat)]), arr.ind=T)
  if(length(location) > 0){
    if (nrow(location) == 1 ){
      temp1 <- RMSE.mat[location]
      temp2 <- S.arr[location[1]]
      temp3 <- Ea.arr[location[2]]
    }
    if (nrow(location) > 1 ){
      temp1 <- RMSE.mat[location[1]]
      temp2 <- S.arr[location[1, 1]]
      temp3 <- Ea.arr[location[1, 2]]
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
  goal.Ea   <- temp3
  if(!is.na(goal.RMSE)){
    RMSE.range <- range(RMSE.mat)
    # Below 'AADTS.arr' and 'Time.pred' is used on the condition 
    #   that S and Ea are finally determined.
    AADTS.arr <- c()
    for(i in 1:length(Year1)){
      temp.T1      <- T1[ Year2 == Year1[i] ]
      temp.DOY     <- DOY[ Year2 == Year1[i] ]
      MeanTemp1    <- temp.T1[which(temp.DOY == goal.S): which(temp.DOY == Time[i])]
      AADTS.arr[i] <- AADTS.fun(goal.Ea, MeanTemp1)
    }
    Time.pred <- c()
    for(i in 1:length(Year1)){
      temp.T1    <- T1[ Year2 == Year1[i] ]
      temp.DOY   <- DOY[ Year2 == Year1[i] ]     
      for(k in which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)){
        temp.AADTS   <- AADTS.fun(goal.Ea, temp.T1[which(temp.DOY == goal.S): k])
        if(temp.AADTS > mean(AADTS.arr)){
          y1           <- AADTS.fun(goal.Ea, temp.T1[which(temp.DOY == goal.S): (k-1)])
          y2           <- mean(AADTS.arr)
          y3           <- temp.AADTS   
          x3           <- temp.DOY[k] - temp.DOY[k-1]      
          if(k > 1){     
            Time.pred[i] <- temp.DOY[k-1] + (y2-y1)/(y3-y1) * x3 
          }          
          if(k == 1){
            Time.pred[i] <- 1
          }  
          break
        }
        extreme.temp <- temp.T1[which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)]
        extreme.AADTS <- AADTS.fun( goal.Ea, extreme.temp )
        if(extreme.AADTS < mean(AADTS.arr)){
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
        if( length(Ea.arr) > 1 ){
          cols <- terrain.colors(200)
          image(S.arr, Ea.arr, RMSE.mat, col = cols, axes = TRUE, cex.axis = 1.5, cex.lab = 1.5, 
            xlab="Starting date (day of year)", 
            ylab=expression(paste(italic(E["a"]), 
            " (kcal" %.% "mol"^{"-1"}, ")", sep="")))
          points(goal.S, goal.Ea, cex=1.5, pch=16, col=2)
          contour(S.arr, Ea.arr, RMSE.mat, levels = round(seq(RMSE.range[1], 
            RMSE.range[2], len = 20), 4), add = TRUE, cex = 1.5, col = "#696969", labcex=1.5)
        }   
        if( length(Ea.arr) == 1 ){
          plot( S.arr, RMSE.mat, cex.axis = 1.5, cex.lab = 1.5, xlim=range(S.arr),
            ylab="RMSE (days)", 
            xlab="Starting date (day of year)", pch=16, cex=1.5, col=4)
          points(goal.S, goal.RMSE, cex=1.5, pch=16, col=2)  
        }
      }
      if(length(S.arr) == 1){
        if(length(Ea.arr) > 1){
          plot( Ea.arr, RMSE.mat, cex.axis = 1.5, cex.lab = 1.5, xlim=range(Ea.arr),
            ylab="RMSE (days)", xlab=expression(paste(italic(E["a"]),
            " (kcal" %.% "mol"^{"-1"}, ")", sep="")), type="l", col=4, lwd=1)
          points(goal.Ea, goal.RMSE, cex=1.5, pch=16, col=2)  
        }
        if(length(Ea.arr) == 1){
          plot( Ea.arr, RMSE.mat, cex.axis = 1.5, cex.lab = 1.5, xlim=range(Ea.arr),
            ylab="RMSE (days)", xlab=expression(paste(italic(E["a"]), 
            " (kcal" %.% "mol"^{"-1"}, ")", sep="")), pch=16, cex=1.5, col=4)
          points(goal.Ea, goal.RMSE, cex=1.5, pch=16, col=2)  
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
           xlab="Observed occurrence time (day of year)",
           ylab="Predicted occurrence time (day of year)", 
           cex.lab=1.5, cex.axis=1.5, pch=16, cex=1.5, col=2)
      abline(0, 1, lwd=1, col=4) 
      points(Time, Time.pred, pch=16, cex=1.5, col=2)  
    }
  }
  if(is.na(goal.RMSE)){
    AADTS.arr <- NA    
    Time.pred <- NA
  }
  list( mAADTS.mat=mAADTS.mat, RMSE.mat=RMSE.mat, AADTS.arr=AADTS.arr, 
        Year=Year1, Time=Time, Time.pred=Time.pred, S=goal.S, Ea=goal.Ea, 
        AADTS=mean(AADTS.arr), RMSE=goal.RMSE, unused.years=unused.years )
}
