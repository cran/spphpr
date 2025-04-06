
predADP2 <- function(S, expr, theta, Year2, DOY, Tmin, Tmax, DOY.ul=120){

  if( !is.numeric(S) | !is.numeric(DOY.ul))
    stop("'S', and 'DOY.ul' should be a number!")
  S <- round(S)
  if(S > 365 | S < 1)
    stop("'S' should be a natural number between 1 and 365!")  
  if(DOY.ul > 365 | DOY.ul < 1){
    stop("'DOY.ul' should be a natural number between 1 and 365!") 
  }

  DOY.ul <- round(DOY.ul)

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
  
  uni.Year2 <- sort( unique(Year2) )
  Time.pred <- c()  
  for(i in 1:length(uni.Year2)){
    temp.Tmin    <- Tmin[ Year2 == uni.Year2[i] ]
    temp.Tmax    <- Tmax[ Year2 == uni.Year2[i] ]
    temp.Tcomb   <- cbind(temp.Tmin, temp.Tmax)  
    temp.DOY     <- DOY[ Year2 == uni.Year2[i] ] 
    for(k in which(temp.DOY == S): which(temp.DOY == DOY.ul)){
      temp.progress <- expr.accum(theta, temp.Tcomb[which(temp.DOY == S): k, ])
      if(temp.progress > 1){
          # y1: The upper base length of the trapezoid
          # y2: The length of the line segment parallel to the two bases of the trapezoid 
          # y3: The lower base length of the trapezoid
          # x3: The height of the trapezoid is temp.DOY[k] - temp.DOY[k-1]     
          y1 <- expr.accum(theta, temp.Tcomb[which(temp.DOY == S): (k-1), ])
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
      extreme.temp     <- temp.Tcomb[which(temp.DOY == S): which(temp.DOY == DOY.ul), ]
      extreme.progress <- expr.accum(theta, extreme.temp)
      if( extreme.progress < 1 ){
          Time.pred[i] <- DOY.ul
      } 
    }      
  }
  list( Year.pred=uni.Year2, Time.pred=Time.pred )
}









