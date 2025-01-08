
ADD <- function( S.pd = NULL, T0.arr, Year1, Time, Year2, DOY, Temp, DOY.ul = 120, 
                 fig.opt = TRUE, S.def = 54, verbose = TRUE){

  DOY.ul <- round(DOY.ul)   
  if( !setequal( unique(Year1), unique(Year2) ) ){
    # Remove the case that 'Year1' isn't in accord with 'Year2'
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
  if(DOY.ul < max(Time)[1]){
    stop("'DOY.ul' should exceed the maximum occurrence time")
  }
  DOYmax <- as.numeric( tapply(DOY, Year2, max) )
  DOYmin <- as.numeric( tapply(DOY, Year2, min) )
  S.arr  <- min(DOYmin)[1]: min(DOY.ul, min(DOYmax)[1])
  if( ( DOY.ul < min(DOYmin)[1] ) | (DOY.ul > min(DOYmax)[1]) ){
    stop("'DOY.ul' should be in the range of 'DOY'")
  }
  AADD.fun <- function(Tb, T){ 
    # Description:
    # AADD.fun is used to calcualte the annual accumulated degree days.
    # Arguments:    
    # Tb: The base temperature
    # T: Daily mean air temperature (in degrees Celsius)
    AADD0 <- 0 
    for(k in 1:length(T)){   
      if(T[k] >= Tb){
        temp.val <- T[k] - Tb
      }
      else{
        temp.val <- 0
      }
      AADD0 <- AADD0 + temp.val
    }
    AADD0
  }  

  if(!is.null(S.pd)){
    if( !is.numeric(S.pd) )
        stop("'S.pd' should be a number!")
    if( S.pd >= min(Time)[1] )
        stop("'S.pd' should be smaller than the minimum observed occurrence time!")
    S.pd     <- round(S.pd)
    goal.S   <- S.pd
    goal.cor <- NA
    S.arr    <- NA
    cor.arr  <- NA
    search.failure <- NA
  }

  if(is.null(S.pd)){
    cor.arr   <- c()
    for(S in S.arr){ 
      ave.temp.arr <- c()
      for(i in 1:length(Year1)){
        temp.T1  <- T1[ Year2 == Year1[i] ]
        temp.DOY <- DOY[ Year2 == Year1[i] ]
        ave.temp <- mean( temp.T1[which(temp.DOY == S)
                          :which(temp.DOY == Time[i])] )
        ave.temp.arr <- c(ave.temp.arr, ave.temp)
      } 
      cor.arr <- c(cor.arr, cor(Time, ave.temp.arr))
    }
    ind1      <- which( cor.arr == min(cor.arr) )
    goal.S    <- S.arr[ind1]
    goal.cor  <- cor.arr[ind1] 
    search.failure <- 0
    if( goal.S >= min(Time)[1]){
        if( S.def >= min(Time)[1] )
            stop("'S.def' should be smaller than the minimum observed occurrence time")

        goal.S         <- S.def  
        goal.cor       <- NA
        invalid.S      <- S.arr[ind1]
        invalid.cor    <- cor.arr[ind1] 
        search.failure <- 1
        warning("The minimum correlation coefficient method failed to find a suitable starting date!")        
    }
  } 

  RMSE.arr  <- c()
  mAADD.arr <- c() 
  len       <- length(T0.arr)
  counter   <- 0
  Time.pred <- c()
  for(T0 in T0.arr){
    counter <- counter + 1
    if(verbose){
      Sys.sleep(.005)
      cat(counter, paste(" of ", len, "\r", sep=""))             
      flush.console()               
      if (counter %% len == 0) cat("\n")
    } 
    AADD.arr <- c()
    for(i in 1:length(Year1)){
      temp.T1     <- T1[ Year2 == Year1[i] ]
      temp.DOY    <- DOY[ Year2 == Year1[i] ]
      MeanTemp1   <- temp.T1[which(temp.DOY == goal.S): which(temp.DOY == Time[i])]
      AADD.arr[i] <- AADD.fun(T0, MeanTemp1)
    }
    mAADD.arr[counter] <- mean(AADD.arr)
    Time.pred          <- c()  
    for(i in 1:length(Year1)){
      temp.T1  <- T1[ Year2 == Year1[i] ]
      temp.DOY <- DOY[ Year2 == Year1[i] ]           
      for(k in which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)){      
        temp.AADD <- AADD.fun(T0, temp.T1[which(temp.DOY == goal.S): k])          
        if(temp.AADD >= mean(AADD.arr)){
          # y1: The upper base length of the trapezoid
          # y2: The length of the line segment parallel to the two bases of the trapezoid 
          # y3: The lower base length of the trapezoid
          # x3: The height of the trapezoid is temp.DOY[k] - temp.DOY[k-1]     
          y1           <- AADD.fun(T0, temp.T1[which(temp.DOY == goal.S): (k-1)])
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
        extreme.temp    <- temp.T1[which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)]
        extreme.AADD    <- AADD.fun( T0, extreme.temp )
        if(extreme.AADD < mean(AADD.arr)){
          Time.pred[i]  <- DOY.ul
        } 
      }   
    }
    RMSE.arr[counter] <- sqrt(sum((Time.pred-Time)^2)/length(Time))
  }
  ind2       <- which(RMSE.arr == min(RMSE.arr[!is.na(RMSE.arr)]))[1]
  goal.RMSE  <- RMSE.arr[ind2]
  goal.T0    <- T0.arr[ind2]
  RMSE.range <- range(RMSE.arr)
  # In the following script, AADD.arr and Time.pred are used on the condition 
  #     that S and T0 are finally determined.
  AADD.arr <- c()
  for(i in 1:length(Year1)){
    temp.T1     <- T1[ Year2 == Year1[i] ]
    temp.DOY    <- DOY[ Year2 == Year1[i] ]
    MeanTemp1   <- temp.T1[which(temp.DOY == goal.S): which(temp.DOY == Time[i])]
    AADD.arr[i] <- AADD.fun(goal.T0, MeanTemp1)
  }
  Time.pred <- c()  
  for(i in 1:length(Year1)){
    temp.T1  <- T1[ Year2 == Year1[i] ]
    temp.DOY <- DOY[ Year2 == Year1[i] ]           
    for(k in which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)){      
      temp.AADD <- AADD.fun(goal.T0, temp.T1[which(temp.DOY == goal.S): k])          
      if(temp.AADD >= mean(AADD.arr)){
        # y1: The upper base length of the trapezoid
        # y2: The length of the line segment parallel to the two bases of the trapezoid 
        # y3: The lower base length of the trapezoid
        # x3: The height of the trapezoid is temp.DOY[k] - temp.DOY[k-1]     
        y1           <- AADD.fun(goal.T0, temp.T1[which(temp.DOY == goal.S): (k-1)])
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
      extreme.temp    <- temp.T1[which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)]
      extreme.AADD    <- AADD.fun( goal.T0, extreme.temp )
      if(extreme.AADD < mean(AADD.arr)){
        Time.pred[i]  <- DOY.ul
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

    if(length(T0.arr) > 1){
      plot( T0.arr, RMSE.arr, cex.axis = 1.5, cex.lab = 1.5, xlim=range(T0.arr),
        ylab="RMSE (days)", type="l", lwd=1, col=4,
        xlab=expression(paste("Base temperature (", degree, "C)", sep="")) )
      points(goal.T0, goal.RMSE, cex=1.5, pch=16, col=2)  
    }
    if(length(T0.arr) == 1){
      plot( T0.arr, RMSE.arr, cex.axis = 1.5, cex.lab = 1.5, xlim=range(T0.arr),
        ylab="RMSE (days)", pch=16, cex=1.5, col=4, 
        xlab=expression(paste("Base temperature (", degree, "C)", sep="")) )
      points(goal.T0, goal.RMSE, cex=1.5, pch=16, col=2)  
    }  
 
    if( is.null(S.pd) ){
      dev.new()
      par1 <- par(family="serif")
      par2 <- par(mar=c(5, 5, 2, 2))
      par3 <- par(mgp=c(3, 1, 0))
      on.exit(par(par1))
      on.exit(par(par2))
      on.exit(par(par3))
      plot(S.arr, cor.arr, cex.lab=1.5, cex.axis=1.5, xlab="Starting date (day-of-year)",           
         ylab="Correlation coefficient", type="l", lwd=1, col=4) 
      if( search.failure == 0)
        points(goal.S, goal.cor, pch=16, cex=1.5, col=2)
      if( search.failure == 1){
        points(invalid.S, invalid.cor, pch=16, cex=1.5, col=2)
        points(invalid.S, invalid.cor, pch=4, cex=2, col=2) 
        abline(v=min(Time)[1], col=1, lty=2)
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

  list( S.arr=S.arr, cor.coef.arr=cor.arr, cor.coef=goal.cor, 
        search.failure=search.failure, mAADD.arr=mAADD.arr, 
        RMSE.arr=RMSE.arr, AADD.arr=AADD.arr, Year=Year1, Time=Time, Time.pred=Time.pred,  
        S=goal.S, T0=goal.T0, AADD=mean(AADD.arr), RMSE=goal.RMSE, unused.years=unused.years )
}



