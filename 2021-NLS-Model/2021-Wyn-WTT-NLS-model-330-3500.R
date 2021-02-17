# NLS is Non-linear Least Squares####
#'A regression technique
#'A signal built from Gaussian curves (can be described by three parameters height, st. dev., and mean)
#'NLS figures out those parameters. The fit is determined by the sum of products between the 
#'KDEs for each on-shore sample times a constant from 0-1. The constants are a percentage out of 100% that each source 
#'sample contributes to a single target sample. We then figure out the fit for each offshore sample based 
#'on the KDEs for the onshore samples. For the pie chart, they're normalized to 100%.

# Load Libraries####
library(MASS)
library(IsoplotR)
library(ggplot2)
library(readxl)
library(minpack.lm)
library(zoo)
library(xlsx)
library(RColorBrewer)

############################################ Load Data

file_name = '2021-Wyn-WTT-NLS-330-3500.xlsx'

# Wynyard Fm (Target) data####
Wynyard <- read_excel(file_name,
                      sheet = "Wynyard Formation")

# Source Data data####
Granites <- read_excel(file_name, 
                         sheet = "Granites")
Mount_Read <- read_excel(file_name, 
                          sheet = "Mount Read Volcanics")
Oonah <- read_excel(file_name, 
                       sheet = "Oonah Formation")
Luina <- read_excel(file_name,  
                    sheet = "Luina Group")
Crimson <- read_excel(file_name, 
                    sheet = "Crimson Creek Formation")
Success <- read_excel(file_name, 
                    sheet = "Success Creek Group")
Arthur <- read_excel(file_name,  
                      sheet = "Arthur Lineament")
Wings <- read_excel(file_name, 
                      sheet = "Wings Sandstone")
Forest <- read_excel(file_name, 
                      sheet = "Forest Conglomerate")
Rocky <- read_excel(file_name, 
                      sheet = "Rocky Cape Group")
Tyennan <- read_excel(file_name,  
                      sheet = "Tyennan Element")

##################### Target sample density objects####
mina = 430
maxa = 1900

target <- data.frame(Wynyard[,c(1,2)])
tmp <- kde(target[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Wynyard.dens <- tmp$y/max(tmp$y)

##################### Source density objects####
source0 <- data.frame(Granites[,c(1,2)])
tmp <- kde(source0[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Granites.dens <- tmp$y/max(tmp$y)

source1 <- data.frame(Mount_Read[,c(1,2)])
tmp <- kde(source1[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Mount_Read.dens <- tmp$y/max(tmp$y)

source2 <- data.frame(Oonah[,c(1,2)])
tmp <- kde(source2[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Oonah.dens <- tmp$y/max(tmp$y)

source3 <- data.frame(Luina[,c(1,2)])
tmp <- kde(source3[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Luina.dens <- tmp$y/max(tmp$y)

source4 <- data.frame(Crimson[,c(1,2)])
tmp <- kde(source4[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Crimson.dens <- tmp$y/max(tmp$y)

source5 <- data.frame(Success[,c(1,2)])
tmp <- kde(source5[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Success.dens <- tmp$y/max(tmp$y)

source6 <- data.frame(Arthur[,c(1,2)])
tmp <- kde(source6[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Arthur.dens <- tmp$y/max(tmp$y)

source7 <- data.frame(Wings[,c(1,2)])
tmp <- kde(source7[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Wings.dens <- tmp$y/max(tmp$y)

source8 <- data.frame(Forest[,c(1,2)])
tmp <- kde(source8[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Forest.dens <- tmp$y/max(tmp$y)

source9 <- data.frame(Rocky[,c(1,2)])
tmp <- kde(source9[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Rocky.dens <- tmp$y/max(tmp$y)

source10 <- data.frame(Tyennan[,c(1,2)])
tmp <- kde(source10[,1], adaptive = TRUE, plot = FALSE, from = mina, to = maxa)
Tyennan.dens <- tmp$y/max(tmp$y)

########### Mixing model for Wynyard ####
#Wynyard ####
data <- data.frame(sample = Wynyard.dens, MR = Mount_Read.dens, OF = Oonah.dens, LG = Luina.dens, CC = Crimson.dens,SC = Success.dens,
                   A = Arthur.dens, W = Wings.dens, FC = Forest.dens, RC = Rocky.dens, TE = Tyennan.dens)

#Find Valid Starting values for NLS
all_runs_data = data.frame(matrix(ncol=22,nrow=0))
xc <- c("start_value","p1","p2","p3","p4","p5","p6","p7","p8","p9","p10",
        "std_p1","std_p2","std_p3","std_p4","std_p5","std_p6",
        "std_p7","std_p8","std_p9","std_p10","rss")



for (i in seq(0.00, 1.00, by = 0.01)){
  start.list <- list(p1 = i, p2 = i, p3 = i, p4 = i, 
                     p5 = i, p6 = i, p7 = i, p8 = i, p9 = i, p10 = i)
  possibleError <-tryCatch(
    nlsLM(data = data, 
          
          sample ~ (MR*p1) + (OF*p2) + (LG*p3) + (CC*p4) + (SC*p5) + (A*p6) + (W*p7) + (FC*p8) + (RC*p9) + (TE*p10),
          
          start = start.list, lower = rep(0,10), upper = rep(1, 10),
          control = nls.lm.control(maxiter = 1000, maxfev = 10000),
          trace = FALSE),
    error=function(e) e
  )
  if(inherits(possibleError, "error")) {
    next} else {
      fit <-nlsLM(data = data, 
                  
                  sample ~ (MR*p1) + (OF*p2) + (LG*p3) + (CC*p4) + (SC*p5) + (A*p6) + (W*p7) + (FC*p8) + (RC*p9) + (TE * p10),
                  
                  start = start.list, lower = rep(0,10), upper = rep(1, 10),
                  control = nls.lm.control(maxiter = 1000, maxfev = 10000),
                  trace = FALSE)
      
      pars <- summary(fit)$par[,1]
      std_error <- summary(fit)$par[,2]
      pred <- predict(fit)
      rss <- deviance(fit)
      new_row <- c(i,pars,std_error,rss)
      all_runs_data <- rbind(all_runs_data, new_row)
      
    } 
  colnames(all_runs_data) <- xc
}

# Select best runs and give one more iteration
all_runs_data = all_runs_data[order(all_runs_data$rss),]
best_runs <-all_runs_data[1:100,]

all_runs_optimized = data.frame(matrix(ncol=22,nrow=0))
best_runs_kdes = data.frame(matrix(ncol=0,nrow=514))

for (i in seq(1,100)){
  p_values <-best_runs[,2:11]
  p_values_it <- p_values[i,]
  v1 = p_values_it[,1]
  v2 = p_values_it[,2]
  v3 = p_values_it[,3]
  v4 = p_values_it[,4]
  v5 = p_values_it[,5]
  v6 = p_values_it[,6]
  v7 = p_values_it[,7]
  v8 = p_values_it[,8]
  v9 = p_values_it[,9]
  v10 = p_values_it[,10]
  
  start.list <- list(p1 = v1, p2 = v2, p3 = v3, p4 = v4, 
                     p5 = v5, p6 = v6, p7 = v7, p8 = v8, p9 = v9, p10 = v10)
  fit <-nlsLM(data = data, 
              
              sample ~ (MR*p1) + (OF*p2) + (LG*p3) + (CC*p4) + (SC*p5) + (A*p6) + (W*p7) + (FC*p8) + (RC*p9) + (TE * p10),
              
              start = start.list, lower = rep(0,10), upper = rep(1, 10),
              control = nls.lm.control(maxiter = 1000, maxfev = 10000),
              trace = FALSE)
  
  pars <- summary(fit)$par[,1]
  std_error <- summary(fit)$par[,2]
  rss <- deviance(fit)
  new_row <- c(i,pars,std_error,rss)
  all_runs_optimized <- rbind(all_runs_optimized, new_row)
  
  
  pred <- predict(fit)
  best_runs_kdes <-cbind.data.frame(best_runs_kdes,pred,deparse.level = 0)
}
colnames(all_runs_optimized) <- xc


plot(x = seq(from = mina, to = maxa,length.out = length(data$sample)),
     y = data$sample, type = 'l',
     main = "Wynyard",
     ylab = "Normalized Density",
     xlab = "Age (Ma)", lwd=2)

no_cols =ncol(best_runs_kdes)
cols <- brewer.pal(4,'Set2')
for (k in (1:no_cols)){
  y_run <-best_runs_kdes[,k]
  lines(x = seq(from = mina, to = maxa,length.out = length(data$sample)),
        y = y_run, col = "purple", lwd=1)
}
legend("topright", legend = c("Data", "100 Best NLS fits"), lty = c(1,1), 
       lwd = c(4,2), col = c("black", "purple"))


###Output####
age_df <- data.frame("age_Ma" = tmp$x, "best_fit_kde" = best_runs_kdes, "target_kde"=Wynyard.dens)
write.csv (age_df,file="model_results_adaptive_kde.csv")

write.csv (all_runs_optimized,file="model_results_adaptive_kde_percents.csv")

################################################################