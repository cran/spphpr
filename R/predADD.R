
predADD <- function(S, T0, AADD, Year2, DOY, Temp, DOY.ul=120){

  if( !is.numeric(S) | !is.numeric(T0) | !is.numeric(DOY.ul))
    stop("'S', 'T0', and 'DOY.ul' should be a number!")
  S    <- round(S)
  if(S > 365 | S < 1)
    stop("'S' should be a natural number between 1 and 365!")  
  if(DOY.ul > 365 | DOY.ul < 1){
    stop("'DOY.ul' should be a natural number between 1 and 365!") 
  }
  DOY.ul <- round(DOY.ul)
  T1     <- Temp
  DOYmax <- as.numeric( tapply(DOY, Year2, max) )
  DOYmin <- as.numeric( tapply(DOY, Year2, min) )
  if( ( DOY.ul < min(DOYmin)[1] ) | (DOY.ul > min(DOYmax)[1]) ){
    stop("'DOY.ul' should be in the range of 'DOY'")
  }
  AADD.fun <- function(Tb, T){ 
    # Description:
    # AADD.fun is used to calcualte the annual accumulated degree days
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
  goal.T0   <- T0
  goal.S    <- S
  Time.pred <- c()  
  uni.Year2 <- sort(unique(Year2))
  for(i in 1:length(uni.Year2)){
    temp.T1  <- T1[ Year2 == uni.Year2[i] ]
    temp.DOY <- DOY[ Year2 == uni.Year2[i] ]           
    for(k in which(temp.DOY == goal.S): which(temp.DOY == DOY.ul)){      
      temp.AADD <- AADD.fun(goal.T0, temp.T1[which(temp.DOY == goal.S): k])          
      if(temp.AADD >= AADD){
        # y1: The upper base length of the trapezoid
        # y2: The length of the line segment parallel to the two bases of the trapezoid 
        # y3: The lower base length of the trapezoid
        # x3: The height of the trapezoid is temp.DOY[k] - temp.DOY[k-1]     
        y1           <- AADD.fun(goal.T0, temp.T1[which(temp.DOY == goal.S): (k-1)])
        y2           <- AADD
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
      if(extreme.AADD < AADD){
        Time.pred[i]  <- DOY.ul
      } 
    }   
  }

  list( Year=uni.Year2, Time.pred=Time.pred )
}





