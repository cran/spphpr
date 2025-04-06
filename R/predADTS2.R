
predADTS2 <- function(S, Ea, AADTS, Year2, DOY, Tmin, Tmax, DOY.ul=120){

  if( !is.numeric(S) | !is.numeric(Ea) | !is.numeric(DOY.ul))
    stop("'S', 'Ea', and 'DOY.ul' should be a number!")
  S    <- round(S)
  if(S > 365 | S < 1)
    stop("'S' should be a natural number between 1 and 365!")  
  if(DOY.ul > 365 | DOY.ul < 1){
    stop("'DOY.ul' should be a natural number between 1 and 365!") 
  }
  DOY.ul <- round(DOY.ul)

  DOYmax <- as.numeric( tapply(DOY, Year2, max) )
  DOYmin <- as.numeric( tapply(DOY, Year2, min) )
  if( (S < min(DOYmin)[1]) | (S > min(DOYmax)[1]) ){
    stop("The starting date should be in the range of 'DOY'")
  }
  if( (DOY.ul < min(DOYmin)[1]) | (DOY.ul > min(DOYmax)[1]) ){
    stop("'DOY.ul' should be in the range of 'DOY'")
  }
  AADTS.fun <- function(Ea, T){ 
    # Description:
    # AADTS.fun is used to calcualte the annual accumulative DTS value
    # Arguments:
    # Ea: The activation free energy (in kcal/mol)
    # T: Vector saving daily min and max air temperatures (in degrees Celsius)
    if(length(T)==2)    T <- matrix(T, nrow=1, ncol=2)
    T1     <- T[,1]
    T2     <- T[,2]
    AADTS0 <- 0
    for(k in 1:nrow(T)){
      Thour <- (T2[k]-T1[k])/2*sin(pi/12*1:24-pi/2)+(T2[k]+T1[k])/2
      Thour <- Thour + 273.15
      Ts <- 298.15 
      # Ts: The reference temperature in K     
      R  <- 1.987   
      # R: The universal gas constant in cal/mol/deg
      DTS <- 0
      for(q in 1:24){
        DTS <- DTS + exp((Ea*1000) * (Thour[q]-Ts)/(R*Thour[q]*Ts))
      }
      AADTS0 <- AADTS0 + DTS/24
    }
    AADTS0
  }  
  goal.S    <- S
  goal.Ea   <- Ea
  Time.pred <- c()
  uni.Year2 <- sort(unique(Year2))
  for(i in 1:length(uni.Year2)){
    temp.Tmin    <- Tmin[ Year2 == uni.Year2[i] ]
    temp.Tmax    <- Tmax[ Year2 == uni.Year2[i] ]
    temp.Tcomb   <- cbind(temp.Tmin, temp.Tmax)
    temp.DOY     <- DOY[ Year2 == uni.Year2[i] ]      
    for(k in which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)){
      ind2       <- which(temp.DOY == goal.S): k    
      temp.AADTS <- AADTS.fun(goal.Ea, temp.Tcomb[ind2, ])   
      if(temp.AADTS > AADTS){
        y1 <- AADTS.fun(goal.Ea, temp.Tcomb[which(temp.DOY == goal.S): (k-1), ])
        y2 <- AADTS
        y3 <- temp.AADTS   
        x3 <- temp.DOY[k] - temp.DOY[k-1]      
        if(k > 1){     
          Time.pred[i] <- temp.DOY[k-1] + (y2-y1)/(y3-y1) * x3 
        }         
        if(k == 1){
          Time.pred[i] <- 1
        }  
        break
      }
      extreme.temp <- temp.Tcomb[which(temp.DOY == goal.S): which(temp.DOY == DOY.ul), ]
      extreme.AADTS <- AADTS.fun( goal.Ea, extreme.temp )
      if(extreme.AADTS < AADTS){
        Time.pred[i] <- DOY.ul
      } 
    }      
  }

  list( Year=uni.Year2, Time.pred=Time.pred )
}



