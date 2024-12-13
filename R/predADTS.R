
predADTS <- function(S, Ea, AADTS, Year2, DOY, Temp, DOY.ul=120){

  if( !is.numeric(S) | !is.numeric(Ea) | !is.numeric(DOY.ul))
    stop("'S', 'Ea', and 'DOY.ul' should be a number!")
  S    <- round(S)
  if(S > 365 | S < 1)
    stop("'S' should be a natural number between 1 and 365!")  
  if(DOY.ul > 365 | DOY.ul < 1){
    stop("'DOY.ul' should be a natural number between 1 and 365!") 
  }
  DOY.ul <- round(DOY.ul)
  T1 <- Temp

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
    # AADTS.fun is used to calcualte the accumulated annual DTS value
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
  goal.S    <- S
  goal.Ea   <- Ea
  Time.pred <- c()
  uni.Year2 <- sort(unique(Year2))
  for(i in 1:length(uni.Year2)){
    temp.T1  <- T1[ Year2 == uni.Year2[i] ]
    temp.DOY <- DOY[ Year2 == uni.Year2[i] ]     
    for(k in which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)){
      temp.AADTS   <- AADTS.fun(goal.Ea, temp.T1[which(temp.DOY == goal.S): k])
      if(temp.AADTS > AADTS){
        y1 <- AADTS.fun(goal.Ea, temp.T1[which(temp.DOY == goal.S): (k-1)])
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
      extreme.temp <- temp.T1[which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)]
      extreme.AADTS <- AADTS.fun( goal.Ea, extreme.temp )
      if(extreme.AADTS < AADTS){
        Time.pred[i] <- DOY.ul
      } 
    }      
  }

  list( Year=uni.Year2, Time.pred=Time.pred )
}



