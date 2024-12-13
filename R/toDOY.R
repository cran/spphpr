
toDOY <- function(Year, Month, Day){  

  if(length(Year) != length(Month) | length(Year) 
     != length(Day) | length(Month) != length(Day)){
    stop("The lengths of 'Year', 'Month' and 'Day' should be the same!")
  }

  Year  <- as.numeric(Year)
  Month <- as.numeric(Month)
  Day   <- as.numeric(Day)

  M1 <- Month
  D1 <- Day
  for( i in 1:length(Month) ){
    if(round(Month[i]) < 0 | round(Month[i]) > 12)
      stop("'Month' should range between 1 and 12!")
    if(round(Day[i]) < 0 | round(Day[i]) > 31)
      stop("'Day' should range between 1 and 31!")
     if(Month[i] < 10) M1[i] <- paste("0", Month[i], sep="")
     if(Day[i]   < 10) D1[i] <- paste("0", Day[i], sep="")      
  }

  x1 <- paste(Year-1, "-12-31", sep="")
  x2 <- paste(Year, "-", Month, "-", Day, sep="")
  as.numeric(difftime(x2, x1, units="days"))
}


