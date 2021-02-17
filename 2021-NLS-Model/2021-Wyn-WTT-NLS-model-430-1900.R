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
library("xlsx")

############################################ Load Data

# Wynyard Fm (Target) data####
Wynyard <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx',
                      sheet = "Wynyard Formation")

# Source Data data####
Mount_Read <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx', 
                          sheet = "Mount Read Volcanics")
Oonah <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx', 
                       sheet = "Oonah Formation")
Luina <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx', 
                    sheet = "Luina Group")
Crimson <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx', 
                    sheet = "Crimson Creek Formation")
Success <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx', 
                    sheet = "Success Creek Group")
Arthur <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx', 
                      sheet = "Arthur Lineament")
Wings <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx', 
                      sheet = "Wings Sandstone")
Forest <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx', 
                      sheet = "Forest Conglomerate")
Rocky <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx', 
                      sheet = "Rocky Cape Group")
Tyennan <- read_excel('2021-Wyn-WTT-NLS-430-1900.xlsx', 
                      sheet = "Tyennan Element")

##################### Target sample density objects####
target <- data.frame(Wynyard[,c(1,2)])
tmp <- kde(target[,1], adaptive = TRUE, plot = FALSE, from = 430, to = 1900)
Wynyard.dens <- tmp$y/max(tmp$y)

##################### Source density objects####
source1 <- data.frame(Mount_Read[,c(1,2)])
tmp <- kde(source1[,1], adaptive = TRUE, plot = FALSE, from = 430, to = 1900)
Mount_Read.dens <- tmp$y/max(tmp$y)

source2 <- data.frame(Luina[,c(1,2)])
tmp <- kde(source1[,1], adaptive = TRUE, plot = FALSE, from = 430, to = 1900)
Luina.dens <- tmp$y/max(tmp$y)

source3 <- data.frame(Crimson[,c(1,2)])
tmp <- kde(source1[,1], adaptive = TRUE, plot = FALSE, from = 430, to = 1900)
Crimson.dens <- tmp$y/max(tmp$y)

source4 <- data.frame(Success[,c(1,2)])
tmp <- kde(source1[,1], adaptive = TRUE, plot = FALSE, from = 430, to = 1900)
Success.dens <- tmp$y/max(tmp$y)

source5 <- data.frame(Arthur[,c(1,2)])
tmp <- kde(source1[,1], adaptive = TRUE, plot = FALSE, from = 430, to = 1900)
Arthur.dens <- tmp$y/max(tmp$y)

source6 <- data.frame(Wings[,c(1,2)])
tmp <- kde(source1[,1], adaptive = TRUE, plot = FALSE, from = 430, to = 1900)
Wings.dens <- tmp$y/max(tmp$y)

source7 <- data.frame(Forest[,c(1,2)])
tmp <- kde(source1[,1], adaptive = TRUE, plot = FALSE, from = 430, to = 1900)
Forest.dens <- tmp$y/max(tmp$y)

source8 <- data.frame(Rocky[,c(1,2)])
tmp <- kde(source1[,1], adaptive = TRUE, plot = FALSE, from = 430, to = 1900)
Rocky.dens <- tmp$y/max(tmp$y)

source9 <- data.frame(Tyennan[,c(1,2)])
tmp <- kde(source1[,1], adaptive = TRUE, plot = FALSE, from = 430, to = 1900)
Tyennan.dens <- tmp$y/max(tmp$y)

########### Mixing model for Wynyard ####
#Wynyard ####
data <- data.frame(sample = Wynyard.dens,
                   MR = Mount_Read.dens, LG = Luina.dens, CC = Crimson.dens,SC = Success.dens,
                   A = Arthur.dens, W = Wings.dens, FC = Forest.dens, RC = Rocky.dens, TE = Tyennan.dens)

start.list <- list(p1 = i, p2 = i, p3 = i, p4 = i, 
                     p5 = i, p6 = i, p7 = i, p8 =i, p9 =i)
fit <- nlsLM(data = data, 
               
               sample ~ (MR*p1) + (LG*p2) + (CC*p3) + (SC*p4) + (A*p5) + (W*p6) + (FC*p7) + (RC*p8) + (RC*p9),
               
               start = start.list, lower = rep(0,9), upper = rep(1, 9),
               control = nls.lm.control(maxiter = 1000, maxfev = 10000),
               trace = TRUE)


pars <- summary(fit)$par[,1]
pred <- predict(fit) 

plot(x = seq(from = 430, to = 4000,length.out = length(data$sample)),
     y = data$sample, type = 'l',
     main = "Wynyard",
     ylab = "Normalized Density",
     xlab = "Age (Ma)", lwd=2)
lines(x = seq(from = 430, to = 4000,length.out = length(data$sample)),
      y = pred, col = 'blue', lwd=2)
legend("topright", legend = c("Data", "NLS fit"), lty = c(1,1), 
       lwd = c(4,4), col = c("black", "blue"))

unknown <- (sum(data$sample) - sum(pred))/sum(data$sample) #unknowns are just the subtracted integrals
known <- sum(pred)/sum(data$sample)
unknown + known
df <- data.frame(
  group = c("?", "Rocky_Cape", "Arthur",
            "Mount_Read", "Luina", "Crimson", "Oonah","Forest", "Wurawina"),
  value = c(unknown, known*(pars/sum(pars)))
)
ggplot(data = df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9",
                             "#009E73", "#F0E442", "#0072B2", "pink","red","green")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank())

###Unknown ID###

difference <- (data$sample - pred)
normalized_difference <- (difference/data$sample)
roll_diff <- rollmean(difference,20)
#abs_diff <- abs(difference)
#roll_abs_diff <- rollmean(abs_diff,50)

plot(x = seq(from = 430, to = 4000,length.out = length(roll_diff)),
     y = (roll_diff), type = 'h',
     main = "roll_diff",
     ylab = "Normalized Density",
     xlab = "Age (Ma)", lwd=2,
     ylim = c(-0.5,1))
lines(x = seq(from = 430, to = 4000,length.out = length(difference)),
      y = data$sample, col = 'red', lwd=2)
lines(x = seq(from = 430, to = 4000,length.out = length(difference)),
      y = pred, col = 'blue', lwd=2)

###Output####
age_df <- data.frame("age_Ma" = tmp$x, "best_fit_kde" = pred, "target_kde"=Wynyard.dens)
write.csv (age_df,file="original_finalmodel_adaptive_kde.csv")
write.xlsx (age_df,file="original_finalmodel_adaptive_kde.xlsx")

write.csv (df,file="original_finalmodel_adaptive_kde_percents.csv")
write.xlsx (df,file="original_finalmodel_adaptive_kde_percents.xlsx")

################################################################