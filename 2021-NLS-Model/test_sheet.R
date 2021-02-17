for (i in seq(0.00, 1.00, by = 0.01)){
  start.list <- list(p1 = i, p2 = i, p3 = i, p4 = i, p5 = i, p6 = i, p7 = i, p8 = i, p9 = i, p10 = i)
  possibleError <-tryCatch(
    nlsLM(data = data, 
          
          sample ~ (MR*p1) + (OF*p2) + (LG*p3) + (CC*p4) + (SC*p5) + (A*p6) + (W*p7) + (FC*p8) + (RC*p9) + (TE * p10),
          
          start = start.list, lower = rep(0,10), upper = rep(1, 10),
          control = nls.lm.control(maxiter = 1000, maxfev = 10000),
          trace = FALSE),
    error=function(e) e
  )
  if(inherits(possibleError, "error")) {next} 
  else {
      print(i)
      
    } 

}

for (i in seq(0.00, 1.00, by = 0.001)){
  start.list <- list(p1 = i, p2 = i, p3 = i, p4 = i, p5 = i, p6 = i, p7 = i, p8 = i, p9 = i, p10 = i)
  possibleError <-tryCatch(
  nlsLM(data = data, 
      
      sample ~ (MR*p1) + (OF*p2) + (LG*p3) + (CC*p4) + (SC*p5) + (A*p6) + (W*p7) + (FC*p8) + (RC*p9) + (TE * p10),
      
      start = start.list, lower = rep(0,10), upper = rep(100, 10),
      control = nls.lm.control(maxiter = 1000, maxfev = 10000),
      trace = FALSE),
  error=function(e) e
  )
  if(inherits(possibleError, "error")) {next}
  else{
    print(i)
  }
}