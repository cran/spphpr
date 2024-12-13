
predADP <- function(S, expr, theta, Year2, DOY, Temp, DOY.ul=120){

  if( !is.numeric(S) | !is.numeric(DOY.ul))
    stop("'S', and 'DOY.ul' should be a number!")
  S <- round(S)
  if(S > 365 | S < 1)
    stop("'S' should be a natural number between 1 and 365!")  
  if(DOY.ul > 365 | DOY.ul < 1){
    stop("'DOY.ul' should be a natural number between 1 and 365!") 
  }

  DOY.ul <- round(DOY.ul)
  T1 <- Temp

  if( S <= 0 ){
    stop("The starting date should be more than 0")
  }

  DOYmax <- as.numeric( tapply(DOY, Year2, max) )
  DOYmin <- as.numeric( tapply(DOY, Year2, min) )
  if( (S < min(DOYmin)[1]) | (S > min(DOYmax)[1]) ){
    stop("The starting date should be in the range of 'DOY'")
  }
  if( (DOY.ul < min(DOYmin)[1]) | (DOY.ul > min(DOYmax)[1]) ){
    stop("'DOY.ul' should be in the range of 'DOY'")
  }

  expr.accum <- function(P, x){
    temp <- 0
    for(i in 1:length(x)){
      temp <- temp + expr(P, x[i])
    }
    return(temp)
  }
  
  uni.Year2 <- sort( unique(Year2) )
  Time.pred <- c()  
  for(i in 1:length(uni.Year2)){
    temp.T1  <- T1[ Year2 == uni.Year2[i] ]
    temp.DOY <- DOY[ Year2 == uni.Year2[i] ]     
    for(k in which(temp.DOY == S): which(temp.DOY == DOY.ul)){
      MeanTemp      <- temp.T1[which(temp.DOY == S): k]
      temp.progress <- expr.accum(theta, MeanTemp)
      if(temp.progress > 1){
          # y1: The upper base length of the trapezoid
          # y2: The length of the line segment parallel to the two bases of the trapezoid 
          # y3: The lower base length of the trapezoid
          # x3: The height of the trapezoid is temp.DOY[k] - temp.DOY[k-1]     
          y1 <- expr.accum(theta, temp.T1[which(temp.DOY == S): (k-1)])
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
      extreme.temp     <- temp.T1[which(temp.DOY == S): which(temp.DOY == DOY.ul)]
      extreme.progress <- expr.accum(theta, extreme.temp)
      if( extreme.progress < 1 ){
          Time.pred[i] <- DOY.ul
      } 
    }      
  }
  list( Year.pred=uni.Year2, Time.pred=Time.pred )
}









