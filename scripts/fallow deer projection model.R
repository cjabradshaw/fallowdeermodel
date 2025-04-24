############################################
## fallow deer (Dama dama) projection model
## CJA Bradshaw
## April 2025
############################################

# load libraries
library(readr)
library(dplyr)
library(plotly)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(viridis)
library(patchwork)
library(deSolve)

# load source functions
source("matrixOperators.r")

# functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# estimate species fecundity allometrically
fecund.pred.func <- function(mass) { # mass in grams
  fecund <- exp(2.719 - 0.211*log(mass))
  return(fecund)
}

# estimate species adult survival allometrically
surv.pred.func <- function(mass) { # mass in grams
  surv <- exp(-exp(-0.5 - 0.25*log(mass)))
  return(surv)
}


#########################################################################
## allometric weighting of vital rates from several cervids to estimate
## fallow deer vital rates

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Longevity parameter setup
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# defines parameters related to deer life span & age distribution for population model
# 'longev' defines maximum lifespan in years
longev <- 18
age.vec <- seq(0, longev, 1)  
lage <- length(age.vec)
sex.ratio <- 0.5 # female only projection but assumed equal male-to-female ratio
stages <- lage  # number of life stages corresponding to the age vector length

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## fertility data
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# deer (standard errors assumed 10% mean fertility)

## calculate correction factors for fallow deer using allometrically predicted fertilities (see main report)
# allometric prediction of fecundity (number of offspring/female/year)
mass.vec <- c(34.1278,52.99999,63.612305,66.1728,82.5,177.5229,240,390)
spp.lab <- c("hog","sika","fallow","chital","rusa","sambar","red","elk")
mass.dat <- data.frame(spp=spp.lab, mass=mass.vec)
mass.dat

fecund.vec <- fecund.pred.func(mass.vec*1000)
fecund.dat <- cbind(mass.dat, fecund=fecund.vec)
fecund.dat

# Cervus elaphus (red deer)
# Mayle, B.A. 1996 Progress in predictive management of deer populations in British woodlands. 
# Forest Ecology and Management 88, 187-198. (doi:10.1016/S0378-1127(96)03824-8)
f1.med <- rep(NA, stages)
f1.med[1:12] <- c(0, 0.425, 0.425, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76)
f1.med <- f1.med * fecund.dat[3,3]/fecund.dat[7,3]  # convert to fallow deer equivalent
f1.se <- f1.med*0.1
f1.med

# Capreolus capreolus (roe deer) (mass = 22.5 kg)
# Mayle, B.A. 1996 Progress in predictive management of deer populations in British woodlands. 
# Forest Ecology and Management 88, 187-198. (doi:10.1016/S0378-1127(96)03824-8)
fecund.roe <- fecund.pred.func(22.5*1000)
f2.med <- rep(NA, stages)
f2.med[1:9] <- c(0, 0.45, 0.57, 0.57, 0.57, 0.57, 0.57, 0.2, 0.2)
f2.med <- f2.med * fecund.dat[2,3]/fecund.roe  # convert to fallow deer equivalent
f2.se <- f2.med*0.1
f2.med

# Cervus elaphus (red deer)
# Benton, T.G., Grant, A. & Clutton-Brock, T.H. 1995 Does environmental stochasticity matter?
# Analysis of red deer life-histories on Rum. Evolutionary Ecology 9, 559-574. (doi:10.1007/BF01237655)
f3.med <- rep(NA, stages)
f3.med[1:19] <- c(0, 0, 0.086, 0.202, 0.241, 0.209, 0.212, 0.215, 0.219, 0.225, 0.206, 0.231, 0.249, 0.19, 0.142,
                  0.082, 0.062, 0.078, 0.025)
f3.med <- f3.med * fecund.dat[3,3]/fecund.dat[7,3]  # convert to fallow deer equivalent
f3.se <- f3.med* 0.1
f3.med

# Dama dama mesopotamica (Persion fallow deer)
# Saltz, D. 1998 A long-term systematic approach to planning reintroductions: the Persian fallowdeer and
# the Arabian oryx in Israel. Animal Conservation 1, 245-252. (doi:10.1111/j.1469-1795.1998.tb00035.x)
f4.med <- rep(NA, stages)
f4.med[1:14] <- c(0, 0, 0.44,	0.29,	0.19,	0.26,	0.31,	0.32,	0.34,	0.34,	0.32,	0.2, 0.18, 0.08)
f4.se <- f4.med * 0.1
f4.med

# Cervus elaphus
# Lowe, V.P.W. 1969 Population dynamics of the red deer (Cervus elaphus L.) on Rhum.
# Journal of Animal Ecology 38, 425-457. (doi:10.2307/2782)
f5.med <- rep(NA, stages)
f5.med[1:16] <- c(0, 0, 0.311, 0.278, 0.302, 0.4, 0.476, 0.358, 0.447, 0.289, 0.283, 0.285, 0.283, 0.282, 0.285, 0.284)
f5.med <- f5.med * fecund.dat[3,3]/fecund.dat[7,3]  # convert to fallow deer equivalent
f5.se <- f5.med * 0.1
f5.med

# Cervus elaphus
# Caughley, G. 1977 Analysis of Vertebrate Populations. Chichester, United Kingdom, John Wiley and Sons
f6.med <- rep(NA, stages)
f6.med[1:12] <- c(0, 0, 0.063, 0.415,	0.4, 0.455, 0.414, 0.486, 0.476, 0.455, 0.5, 0.375)
f6.med <- f6.med * fecund.dat[3,3]/fecund.dat[7,3]  # convert to fallow deer equivalent
f6.se <- f6.med * 0.1
f6.med

# plot fertility data for each of the above
par(mfrow=c(3,2))
plot(age.vec, f1.med, type="l", xlab="age", ylab="fertility", ylim=c(0,1))
lines(age.vec, f1.med - 1.96*f1.se, lty=2, col="red")
lines(age.vec, f1.med + 1.96*f1.se, lty=2, col="red")

plot(age.vec, f2.med, type="l", xlab="age", ylab="fertility", ylim=c(0,1))
lines(age.vec, f2.med - 1.96*f2.se, lty=2, col="red")
lines(age.vec, f2.med + 1.96*f2.se, lty=2, col="red")

plot(age.vec, f3.med, type="l", xlab="age", ylab="fertility", ylim=c(0,1))
lines(age.vec, f3.med - 1.96*f3.se, lty=2, col="red")
lines(age.vec, f3.med + 1.96*f3.se, lty=2, col="red")

plot(age.vec, f4.med, type="l", xlab="age", ylab="fertility", ylim=c(0,1))
lines(age.vec, f4.med - 1.96*f4.se, lty=2, col="red")
lines(age.vec, f4.med + 1.96*f4.se, lty=2, col="red")

plot(age.vec, f5.med, type="l", xlab="age", ylab="fertility", ylim=c(0,1))
lines(age.vec, f5.med - 1.96*f5.se, lty=2, col="red")
lines(age.vec, f5.med + 1.96*f5.se, lty=2, col="red")

plot(age.vec, f6.med, type="l", xlab="age", ylab="fertility", ylim=c(0,1))
lines(age.vec, f6.med - 1.96*f6.se, lty=2, col="red")
lines(age.vec, f6.med + 1.96*f6.se, lty=2, col="red")

par(mfrow=c(1,1))

# data.frame of age-specific mean survival probabilities
fertmn.dat <- data.frame(f1.med, f2.med, f3.med, f4.med, f5.med, f6.med)
fert.mn <- apply(fertmn.dat, 1, median, na.rm=T)
fertsepr.mn <- rep(0.1,length(age.vec))

# data.frame of approximated lower CI age specific survival values
fertlo.dat <- data.frame(f1.med - 1.96*f1.se, f2.med - 1.96*f2.se, f3.med - 1.96*f3.se, 
                         f4.med - 1.96*f4.se, f5.med - 1.96*f5.se, f6.med - 1.96*f6.se)
colnames(fertlo.dat) <- c("f1lo", "f2lo", "f3lo", "f4lo", "f5lo", "f6lo")

# data.frame of approximated upper CI age specific survival values
fertup.dat <- data.frame(f1.med + 1.96*f1.se, f2.med + 1.96*f2.se, f3.med + 1.96*f3.se, 
                         f4.med + 1.96*f4.se, f5.med + 1.96*f5.se, f6.med + 1.96*f6.se)
colnames(fertup.dat) <- c("f1up", "f2up", "f3up", "f4up", "f5up", "f6up")

plot(age.vec, fert.mn, type="l", ylim=c(0,1))
lines(age.vec, fert.mn - 1.96*fertsepr.mn, lty=2, col="red")
lines(age.vec, ifelse((fert.mn + 1.96*fertsepr.mn) > 1, 1, fert.mn + 1.96*fertsepr.mn), lty=2, col="red")

# fit function to mean fertility to smooth
# power series: Y=a1*X^b1 + c1*X^d1
fert.dat <- data.frame(age.vec, fert.mn)
colnames(fert.dat) <- c("agemod", "fertmn")
fert.dat.trunc <- fert.dat[-c(1:2),]
fert.dat.trunc

param.init <- c(-0.003852,
                1.798,
                0.3622,
                0.2421) # parameter initialisation
fit.ps <- nls(fertmn ~ a1*agemod^b1 + c1*agemod^d1, 
                data = fert.dat.trunc,
                algorithm = "default",
                start = c(a1 = param.init[1], b1 = param.init[2], c1 = param.init[3], d1 = param.init[4]),
                trace = TRUE,      
                nls.control(maxiter = 100000, tol = 1e-05, minFactor = 1/1024))
coef(fit.ps) # fitted coefficients

# plot smoothed fertility curve
plot(age.vec, fert.dat[,2], pch=19, xlab="age", type="l", ylab="fertility", ylim=c(0, 0.5),col="black")
age.pred.vec <- seq(0,longev,1)
pred.fert <- as.numeric(coef(fit.ps)[1]) * age.pred.vec^as.numeric(coef(fit.ps)[2]) +
  as.numeric(coef(fit.ps)[3]) * age.pred.vec^as.numeric(coef(fit.ps)[4])
pred.fert[1:2] <- 0
lines(age.pred.vec, pred.fert, lty=2, lwd=1, col="red")



##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## survival data
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
surv.vec <- surv.pred.func(mass.vec*1000)
surv.dat <- cbind(mass.dat, surv=surv.vec)
surv.dat

# Dama dama (fallow deer)
# Focardi, S., Toso, S. & Pecchioli, E. 1996 The population modelling of fallow deer and wild boar in
# a Mediterranean ecosystem. Forest Ecology and Management 88, 7-14. (doi:10.1016/S0378-1127(96)03804-2)
s1.med <- rep(NA, stages)
s1.med[1:19] <- c(0.72,	rep(0.93, 18))
s1.se <- s1.med*0.1
s1.med

# Cervus elaphus (red deer)
# Mayle, B.A. 1996 Progress in predictive management of deer populations in British woodlands. 
# Forest Ecology and Management 88, 187-198. (doi:10.1016/S0378-1127(96)03824-8).
s2.med <- rep(NA, stages)
s2.med[1:13] <- c(0.48,	0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.3, 0.3, 0.3, 0.3,	0.3, 0.3)
s2.med <- s2.med * surv.dat[3,3]/surv.dat[7,3]  # convert to fallow deer equivalent
s2.se <- s2.med*0.1
s2.med

# Capreolus capreolus (roe deer); mass = 22.5 kg
# Mayle, B.A. 1996 Progress in predictive management of deer populations in British woodlands. 
# Forest Ecology and Management 88, 187-198. (doi:10.1016/S0378-1127(96)03824-8).
surv.roe <- surv.pred.func(22.5*1000)
s3.med <- rep(NA, stages)
s3.med[1:9] <- c(0.88,	0.9, 0.95, 0.95, 0.95, 0.95, 0.95, 0.3,	0.3)
s3.med <- s3.med * surv.dat[3,3]/surv.roe  # convert to fallow deer equivalent
s3.se <- s3.med* 0.1
s3.med

# Cervus elaphus (red deer)
# Benton, T.G., Grant, A. & Clutton-Brock, T.H. 1995 Does environmental stochasticity matter?
# Analysis of red deer life-histories on Rum. Evolutionary Ecology 9, 559-574. (doi:10.1007/BF01237655)
s4.med <- rep(NA, stages)
s4.med[1:19] <- c(0.902, 0.962,	0.956, 0.956, 0.949, 0.939,	0.971, 0.982, 0.942, 0.941,	0.881, 0.82, 0.837, 0.755, 
            0.639, 0.498,	0.405, 0.238,	0)
s4.med <- s4.med * surv.dat[3,3]/surv.dat[7,3]  # convert to fallow deer equivalent
s4.se <- s4.med * 0.1
s4.med

# Dama dama mesopotamica (Persian fallow)
# Saltz, D. 1998 A long-term systematic approach to planning reintroductions: the Persian fallowdeer and
# the Arabian oryx in Israel. Animal Conservation 1, 245-252. (doi:10.1111/j.1469-1795.1998.tb00035.x)
s5.med <- rep(NA, stages)
s5.med[1:14] <- c(0.66, 0.93,	0.93,	0.93,	0.93,	0.93,	0.93,	0.93,	0.93,	0.91, 0.87, 0.79, 0.63, 0.31)
s5.se <- s5.med * 0.1
s5.med

# Cervus elaphus (red deer)
# Lowe, V.P.W. 1969 Population dynamics of the red deer (Cervus elaphus L.) on Rhum.
# Journal of Animal Ecology 38, 425-457. (doi:10.2307/2782)
s6.med <- rep(NA, stages)
s6.med[1:16] <- c(0.863,	0.9027,	0.8922,	0.8792,	0.8626,	0.8407,	0.8105,	0.4984,	0.3273,	0.8588,	0.8354,	0.8025,
            0.7532,	0.6712,	0.5076,	0)
s6.med <- s6.med * surv.dat[3,3]/surv.dat[7,3]  # convert to fallow deer equivalent
s6.se <- s6.med * 0.1
s6.med

# Odocoileus virginianus (white tailed deer); mass = 65.32 kg
# Wiskirchen, K.H., Jacobsen, T.C., Ditchkoff, S.S., Demarais, S. & Grand, J.B. 2023 Adult white-tailed
# deer survival in hunted populations on public and private lands. 
# Wildlife Society Bulletin 47, e1391. (doi:10.1002/wsb.1391).
surv.whitetail <- surv.pred.func(65.32*1000)
s7.med <- rep(NA, stages)
s7.med[1:16] <- c(NA,	0.733, 0.733,	0.733, 0.83, 0.83,	0.83,	0.83,	0.83,	0.83,	0.83,	0.83,	0.83,	0.83,	0.83,	0.83)
s7.med <- s7.med * surv.dat[3,3]/surv.whitetail  # convert to fallow deer equivalent
s7.se <- s7.med * 0.1
s7.med

# Cervus canadensis (elk)
# Unsworth, J.W., Kuck, L., Scott, M.D. & Garton, E.O. 1993 Elk mortality in the Clearwater Drainage
# of northcentral Idaho. Journal of Wildlife Management 57, 495-502. (doi:10.2307/3809273)
s8.med <- rep(NA, stages)
s8.med[2:9] <- rep(0.886, 8)
s8.med <- s8.med * surv.dat[3,3]/surv.dat[8,3]  # convert to fallow deer equivalent
s8.se <- s8.med * 0.1
s8.med

# Odocoileus hemionus (mule deer - site PHU1); mass = 120 kg
# Pac, D.F., Mackie, R.J. & Jorgensen, H.E. 1991 Mule Deer Population Organization,
# Behavior and Dynamics in a Northern Rocky Mountain Environment. (Helena, Montana, USA)
surv.mule <- surv.pred.func(120*1000)
s9.med <- rep(NA, stages)
s9.med[2:16] <- c(0.8843,	1, 0.9057, 0.9486,	0.8971,	0.9229,	0.9143,	0.9357,	0.8586,	0.9314,	0.82,
                0.6486,	0.6443,	0.3871,	0.5071)
s9.med <- s9.med * surv.dat[3,3]/surv.mule  # convert to fallow deer equivalent
s9.se <- s9.med * 0.1
s9.med

# Odocoileus hemionus (mule deer - site PHU4); mass = 120 kg
# Pac, D.F., Mackie, R.J. & Jorgensen, H.E. 1991 Mule Deer Population Organization,
# Behavior and Dynamics in a Northern Rocky Mountain Environment. (Helena, Montana, USA)
s10.med <- rep(NA, stages)
s10.med[2:14] <- c(0.8457, 0.82, 0.9057, 0.79, 0.8329, 0.79, 0.8629, 0.7514, 0.7129, 0.6271, 0.6271,
                 0.5029, 0.3529)
s10.med <- s10.med * surv.dat[3,3]/surv.mule  # convert to fallow deer equivalent
s10.se <- s10.med * 0.1
s10.med

# plot survival data for each of the above, 
par(mfrow=c(2,5))
plot(age.vec, s1.med, type="l", xlab="age", ylab="survival", ylim=c(0,1), main = "s1")
lines(age.vec, s1.med - 1.96*s1.se, lty=2, col="red")
lines(age.vec, s1.med + 1.96*s1.se, lty=2, col="red")

plot(age.vec, s2.med, type="l", xlab="age", ylab="survival", ylim=c(0,1), main = "s2")
lines(age.vec, s2.med - 1.96*s2.se, lty=2, col="red")
lines(age.vec, s2.med + 1.96*s2.se, lty=2, col="red")

plot(age.vec, s3.med, type="l", xlab="age", ylab="survival", ylim=c(0,1), main = "s3")
lines(age.vec, s3.med - 1.96*s3.se, lty=2, col="red")
lines(age.vec, s3.med + 1.96*s3.se, lty=2, col="red")

plot(age.vec, s4.med, type="l", xlab="age", ylab="survival", ylim=c(0,1), main = "s4")
lines(age.vec, s4.med - 1.96*s4.se, lty=2, col="red")
lines(age.vec, s4.med + 1.96*s4.se, lty=2, col="red")

plot(age.vec, s5.med, type="l", xlab="age", ylab="survival", ylim=c(0,1), main = "s5")
lines(age.vec, s5.med - 1.96*s5.se, lty=2, col="red")
lines(age.vec, s5.med + 1.96*s5.se, lty=2, col="red")

plot(age.vec, s6.med, type="l", xlab="age", ylab="survival", ylim=c(0,1), main = "s6")
lines(age.vec, s6.med - 1.96*s6.se, lty=2, col="red")
lines(age.vec, s6.med + 1.96*s6.se, lty=2, col="red")

plot(age.vec, s7.med, type="l", xlab="age", ylab="survival", ylim=c(0,1), main = "s7")
lines(age.vec, s7.med - 1.96*s7.se, lty=2, col="red")
lines(age.vec, s7.med + 1.96*s7.se, lty=2, col="red")

plot(age.vec, s8.med, type="l", xlab="age", ylab="survival", ylim=c(0,1), main = "s8")
lines(age.vec, s8.med - 1.96*s8.se, lty=2, col="red")
lines(age.vec, s8.med + 1.96*s8.se, lty=2, col="red")

plot(age.vec, s9.med, type="l", xlab="age", ylab="survival", ylim=c(0,1), main = "s9")
lines(age.vec, s9.med - 1.96*s9.se, lty=2, col="red")
lines(age.vec, s9.med + 1.96*s9.se, lty=2, col="red")

plot(age.vec, s10.med, type="l", xlab="age", ylab="survival", ylim=c(0,1), main = "s10")
lines(age.vec, s10.med - 1.96*s10.se, lty=2, col="red")
lines(age.vec, s10.med + 1.96*s10.se, lty=2, col="red")

par(mfrow = c(1, 1))  # Reset to single plot layout


# data.frame of age-specific mean survival probabilities
survmn.dat <- data.frame(s1.med, s2.med, s3.med, s4.med, s5.med, s6.med, s7.med, s8.med, s9.med, s10.med)
surv.mn <- apply(survmn.dat, 1, median, na.rm=T)
survsepr.mn <- rep (0.1,length(age.vec))

# data.frame of approximated lower CI age specific survival values
survlo.dat <- data.frame(s1.med - 1.96*s1.se, s2.med - 1.96*s2.se, s3.med - 1.96*s3.se, 
                         s4.med - 1.96*s4.se, s5.med - 1.96*s5.se, s6.med - 1.96*s6.se, 
                         s7.med - 1.96*s7.se, s8.med - 1.96*s8.se, s9.med - 1.96*s9.se, s10.med - 1.96*s10.se)
colnames(survlo.dat) <- c("s1lo", "s2lo", "s3lo", "s4lo", "s5lo", "s6lo", "s7lo", "s8lo","s9lo", "s10lo")

# data.frame of approximated upper CI age specific survival values
survup.dat <- data.frame(s1.med + 1.96*s1.se, s2.med + 1.96*s2.se, s3.med + 1.96*s3.se, 
                         s4.med + 1.96*s4.se, s5.med + 1.96*s5.se, s6.med + 1.96*s6.se, 
                         s7.med + 1.96*s7.se, s8.med + 1.96*s8.se, s9.med + 1.96*s9.se, s10.med + 1.96*s10.se)
colnames(survup.dat) <- c("s1up", "s2up", "s3up", "s4up", "s5up", "s6up", "s7up", "s8up","s9up", "s10up")

plot(age.vec, surv.mn, type="l", ylim=c(0,1))
lines(age.vec, surv.mn - 1.96*survsepr.mn, lty=2, col="red")
lines(age.vec, ifelse((surv.mn + 1.96*survsepr.mn) > 1, 1, surv.mn + 1.96*survsepr.mn), lty=2, col="red")

# Siler hazard model
mort.dat <- data.frame(age.vec, 1-surv.mn)
colnames(mort.dat) <- c("agemod", "mortmn")
mort.dat

a1 <- 0.2 # initial fawn mortality rate (also known as αt)
b1 <- 1.131 # rate of mortality decline (also known as bt)
a2 <- 0.013 # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1 - 0.95 # initial adult mortality rate (also known as βt)
b3 <- 0.086 # rate of mortality increase
param.init <- c(a1, b1, a2, a3, b3) # parameter initialisation

fit.siler <- nls(mortmn ~ a1 * exp(-b1*agemod) + a2 + a3 * exp(b3 * agemod), 
                 data = mort.dat,
                 algorithm = "default",
                 start = c(a1 = param.init[1], b1 = param.init[2], a2 = param.init[3],
                           a3 = param.init[4], b3 = param.init[5]),
                 trace = TRUE,      
                 nls.control(maxiter = 100000, tol = 1e-05, minFactor = 1/1024))
coef(fit.siler) # fitted coefficients

# plot smoothed mortality curve
plot(age.vec, mort.dat[,2], pch=19, xlab="age", type="l", ylab="mortality", ylim=c(0, 1),col="black")
age.pred.vec <- seq(0,longev,1)
pred.mort <- as.numeric(coef(fit.siler)[1]) * exp(-as.numeric(coef(fit.siler)[2]) * age.pred.vec) +
                                                   as.numeric(coef(fit.siler)[3]) + as.numeric(coef(fit.siler)[4]) *
                                                   exp(as.numeric(coef(fit.siler)[5]) * age.pred.vec)
lines(age.pred.vec, pred.mort, lty=2, lwd=1, col="red")


# plot smoothed survival curve
plot(age.vec, surv.mn, pch=19, xlab="age", type="l", ylab="survival", ylim=c(0, 1),col="black")
age.pred.vec <- seq(0,longev,1)
pred.surv <- 1 - pred.mort
pred.surv.se <- pred.surv * 0.1
lines(age.pred.vec, pred.surv, lty=2, lwd=1, col="red")
points(age.vec, survmn.dat[,1],pch=19,cex=0.5,col="light grey")
points(age.vec, survmn.dat[,2],pch=19,cex=0.5,col="light grey")
points(age.vec, survmn.dat[,3],pch=19,cex=0.5,col="light grey")
points(age.vec, survmn.dat[,4],pch=19,cex=0.5,col="light grey")
points(age.vec, survmn.dat[,5],pch=19,cex=0.5,col="light grey")
points(age.vec, survmn.dat[,6],pch=19,cex=0.5,col="light grey")
points(age.vec, survmn.dat[,7],pch=19,cex=0.5,col="light grey")
points(age.vec, survmn.dat[,8],pch=19,cex=0.5,col="light grey")
points(age.vec, survmn.dat[,9],pch=19,cex=0.5,col="light grey")
points(age.vec, survmn.dat[,10],pch=19,cex=0.5,col="light grey")
points(age.vec, survlo.dat[,1],pch=19,cex=0.5,col="pink")
points(age.vec, survlo.dat[,2],pch=19,cex=0.5,col="pink")
points(age.vec, survlo.dat[,3],pch=19,cex=0.5,col="pink")
points(age.vec, survlo.dat[,4],pch=19,cex=0.5,col="pink")
points(age.vec, survlo.dat[,5],pch=19,cex=0.5,col="pink")
points(age.vec, survlo.dat[,6],pch=19,cex=0.5,col="pink")
points(age.vec, survlo.dat[,7],pch=19,cex=0.5,col="pink")
points(age.vec, survlo.dat[,8],pch=19,cex=0.5,col="pink")
points(age.vec, survlo.dat[,9],pch=19,cex=0.5,col="pink")
points(age.vec, survlo.dat[,10],pch=19,cex=0.5,col="pink")
points(age.vec, survup.dat[,1],pch=19,cex=0.5,col="light green")
points(age.vec, survup.dat[,2],pch=19,cex=0.5,col="light green")
points(age.vec, survup.dat[,3],pch=19,cex=0.5,col="light green")
points(age.vec, survup.dat[,4],pch=19,cex=0.5,col="light green")
points(age.vec, survup.dat[,5],pch=19,cex=0.5,col="light green")
points(age.vec, survup.dat[,6],pch=19,cex=0.5,col="light green")
points(age.vec, survup.dat[,7],pch=19,cex=0.5,col="light green")
points(age.vec, survup.dat[,8],pch=19,cex=0.5,col="light green")
points(age.vec, survup.dat[,9],pch=19,cex=0.5,col="light green")
points(age.vec, survup.dat[,10],pch=19,cex=0.5,col="light green")


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## deterministic projection without density feedback
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## create matrix
popmat <- matrix(data = 0, nrow=stages, ncol=stages)
diag(popmat[2:stages,]) <- pred.surv[-length(pred.surv)]
popmat[stages,stages] <- pred.surv[stages]
popmat[1,] <- pred.fert
colnames(popmat) <- age.vec
rownames(popmat) <- age.vec
popmat.orig <- popmat ## save original matrix

## matrix properties
max.lambda(popmat.orig) ## 1-yr lambda
max.r(popmat.orig) # rate of population change, 1-yr
# from allometric prediction: r = 0.2703408 (should be substantially higher than measured values)
ssd <- stable.stage.dist(popmat.orig) ## stable stage distribution
plot(age.vec, ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
gen.l <- G.val(popmat.orig, longev) # mean generation length
gen.l

## allometric prediction of equilibrium densities
## for mammal herbivores (indviduals per km^2)
## D = 10^(4.196 - 0.74*log10(M)) (Damuth 1981 Population density and body size in mammals Nature 290:699–700. doi:10.1038/290699a0 Google Scholar)
Deq <- 10^(4.196 - 0.74*log10(1000*mass.dat$mass))
Deq
Deq.dat <- data.frame(spp=spp.lab, mass=mass.dat$mass, Deq=Deq)
Deq.dat

## high habitat suitability (> 0.5) cells in SDM prediction
nhighsuit.cells <- 40180
cell.dimlon <- 0.008333333333333328014
cell.size <- (110.574*cell.dimlon)^2
cell.size

# convert to population size
curr.pop.Deq <- round(nhighsuit.cells * cell.size * Deq.dat[3,3], 0)
curr.pop.Deq

## initial population vector
curr.pop <- curr.pop.Deq * sex.ratio  # multiply population estimate by sex ratio to get female population
# create the initial population vector using the stable stage distribution
init.vec <- curr.pop * stable.stage.dist(popmat)
# Plot the initial population distribution across ages
plot(age.vec, init.vec, xlab = "age (years)", ylab = "Nf", type = "l")  # Line plot of the female population across ages


#################################
## population projection model ##
#################################
# set time parameters for projection
yr.now <- 2025  # current year (starting point of the projection)
yr.end <- 2050  # end year for the projection
t <- (yr.end - yr.now)  # Total number of years to project (e.g., 24 years)

# use original population matrix (before any modifications)
popmat <- popmat.orig

## storage
# create a matrix to store population counts for each stage over time
# Rows represent age stages, and columns represent time (years)
n.mat <- matrix(0, nrow = stages, ncol = (t + 1))  # Initialize with zeros

# set initial population vector (from stable stage distribution) as first column in n.mat
n.mat[, 1] <- init.vec  # represents population at year 1

## population projection loop
# loop through each year and project population using matrix multiplication
for (i in 1:t) {
  n.mat[, i + 1] <- popmat %*% n.mat[, i]  # multiply population matrix by current population vector
}

## plot results
## shows how total population of females evolves over projection period, 
## based on fertility, survival, and population dynamics
# create vector of years for x-axis of the plot
yrs <- seq(yr.now, yr.end, 1)

# plot total population (sum across all age groups) over time
plot(yrs, as.vector(colSums(n.mat)), type = "l", 
     xlab = "year", ylab = "Nf", 
     xlim = c(yr.now, yr.end))  # set x-axis limits to projection period

# export for external plotting
det.proj.out <- data.frame(year=yrs, N=2*round(as.vector(colSums(n.mat)),0))
write.csv(det.proj.out, file="detprojout.csv", row.names=FALSE)


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## stochastic projection without density feedback
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
s.arr <- m.arr <- array(data=NA, dim=c(t+1, longev+1, iter))

for (e in 1:iter) {
  popmat <- popmat.orig
  
  n.mat <- s.mat <- m.mat <- matrix(0, nrow=longev+1,ncol=(t+1))
  n.mat[,1] <- init.vec
  s.mat[,1] <- pred.surv
  m.mat[,1] <- pred.fert
  
  for (i in 1:t) {
    # stochastic survivals following beta distribution
    s.alpha <- estBetaParams(pred.surv, pred.surv.se^2)$alpha
    s.beta <- estBetaParams(pred.surv, pred.surv.se^2)$beta
    s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
    
    # stochastic fertilty sampler (Gaussian)
    fert.stch <- rnorm(length(popmat[,1]), pred.fert, 0.1*pred.fert)
    m.arr[i,,e] <- ifelse(fert.stch < 0, 0, fert.stch)
    
    totN.i <- sum(n.mat[,i], na.rm=T)

    diag(popmat[2:(stages),]) <- (s.stoch[-(stages)])
    popmat[stages,stages] <- (s.stoch[stages])
    popmat[1,] <- m.arr[i,,e]
    n.mat[,i+1] <- popmat %*% n.mat[,i]
    
    s.arr[i,,e] <- s.stoch
    
  } # end i loop
  
  n.sums.mat[e,] <- as.vector(colSums(n.mat))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,n.md,type="l", main = "", xlab="year", ylab="N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,n.up,lty=2,col="red",lwd=1.5)

s.add <- m.add  <- rep(0, stages)
for (m in 1:iter) {
  s.add <- rbind(s.add, s.arr[ceiling(gen.l):(t+1),,m])
  m.add <- rbind(m.add, m.arr[ceiling(gen.l):(t+1),,m])
  if (m %% itdiv==0) print(m)
}
s.add <- s.add[-1,]
m.add <- m.add[-1,]

s.md <- apply(s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
s.up <- apply(s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
s.lo <- apply(s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(age.vec,s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(s.lo),1.05*max(s.up)))
lines(age.vec,s.lo,lty=2,col="red",lwd=1.5)
lines(age.vec,s.up,lty=2,col="red",lwd=1.5)

m.md <- apply(m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
m.up <- apply(m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
m.lo <- apply(m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(age.vec,m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(m.lo),1.05*max(m.up)))
lines(age.vec,m.lo,lty=2,col="red",lwd=1.5)
lines(age.vec,m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))

# export for external plotting
stoch.proj.out <- data.frame(year=yrs, N=2*round(n.md,0), N.up=2*round(n.up,0), N.lo=2*round(n.lo,0))
write.csv(stoch.proj.out, file="stochprojout.csv", row.names=FALSE)



##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## stochastic projection with catastrophes and without density feedback
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10
cat.pr <- 0.14 # probability of a catastrophe occurring per generation

n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
s.arr <- m.arr <- array(data=NA, dim=c(t+1, longev+1, iter))

for (e in 1:iter) {
  popmat <- popmat.orig
  
  n.mat <- s.mat <- m.mat <- matrix(0, nrow=longev+1,ncol=(t+1))
  n.mat[,1] <- init.vec
  s.mat[,1] <- pred.surv
  m.mat[,1] <- pred.fert
  
  for (i in 1:t) {
    # stochastic survivals following beta distribution
    s.alpha <- estBetaParams(pred.surv, pred.surv.se^2)$alpha
    s.beta <- estBetaParams(pred.surv, pred.surv.se^2)$beta
    s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
    
    # catastrophe
    if (rbinom(1, 1, cat.pr/gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      s.stoch <- s.stoch * (rbeta(1, cat.alpha, cat.beta))
    }
    
    # stochastic fertilty sampler (Gaussian)
    fert.stch <- rnorm(length(popmat[,1]), pred.fert, 0.1*pred.fert)
    m.arr[i,,e] <- ifelse(fert.stch < 0, 0, fert.stch)
    
    totN.i <- sum(n.mat[,i], na.rm=T)
 
    diag(popmat[2:(stages),]) <- (s.stoch[-(stages)])
    popmat[stages,stages] <- (s.stoch[stages])
    popmat[1,] <- m.arr[i,,e]
    n.mat[,i+1] <- popmat %*% n.mat[,i]
    
    s.arr[i,,e] <- s.stoch
    
  } # end i loop
  
  n.sums.mat[e,] <- as.vector(colSums(n.mat))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,n.md,type="l", main = "", xlab="year", ylab="N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,n.up,lty=2,col="red",lwd=1.5)

s.add <- m.add  <- rep(0, stages)
for (m in 1:iter) {
  s.add <- rbind(s.add, s.arr[ceiling(gen.l):(t+1),,m])
  m.add <- rbind(m.add, m.arr[ceiling(gen.l):(t+1),,m])
  if (m %% itdiv==0) print(m) 
}
s.add <- s.add[-1,]
m.add <- m.add[-1,]

s.md <- apply(s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
s.up <- apply(s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
s.lo <- apply(s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(age.vec,s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(s.lo),1.05*max(s.up)))
lines(age.vec,s.lo,lty=2,col="red",lwd=1.5)
lines(age.vec,s.up,lty=2,col="red",lwd=1.5)

m.md <- apply(m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
m.up <- apply(m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
m.lo <- apply(m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(age.vec,m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(m.lo),1.05*max(m.up)))
lines(age.vec,m.lo,lty=2,col="red",lwd=1.5)
lines(age.vec,m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))

# export for external plotting
stochcat.proj.out <- data.frame(year=yrs, N=2*round(n.md,0), N.up=2*round(n.up,0), N.lo=2*round(n.lo,0))
write.csv(stochcat.proj.out, file="stochcatprojout.csv", row.names=FALSE)



##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## stochastic projection with catastrophes and density feedback
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# density feedback on survival
K.pop <- round(1.5 * curr.pop, 0)  # carrying capacity (e.g., 1.5 times current population)
K.pop

# survival multipliers for different population densities
surv.mult.start <- 1.0
surv.mult.lo <- 0.997
surv.mult.mid <- 0.985
surv.mult.end <- 0.875

K.lo <- 1.1 * curr.pop
K.mid <- 1.25 * curr.pop

# Survival multipliers based on population sizes
K.vec <- c(curr.pop,K.lo,K.mid,K.pop)
surv.mult.vec <- c(surv.mult.start, surv.mult.lo, surv.mult.mid, surv.mult.end)
plot(K.vec, surv.mult.vec, pch=19)# Plot the relationship between population size and survival multipliers

## fit logistic function to survival multipliers based on population size
# logistic function: y = a / (exp((b - X) / c) + 1)
DD.dat <- data.frame(K.vec, surv.mult.vec)

# initial parameter estimation for logistic model
param.init.est <- as.numeric(getInitial(surv.mult.vec ~ SSlogis(K.vec, a, b, c), dat=DD.dat))
param.init.est

# Fit a logistic function to model survival multipliers based on population size
fit.lgst <- nls(
  surv.mult.vec ~ a / (exp((b - K.vec) / c) + 1), 
  data = DD.dat,
  algorithm = "port",
  start = c(a = param.init.est[1], b = param.init.est[2], c = param.init.est[3]),
  trace = TRUE,      
  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/2048)
)

a.lp <- coef(fit.lgst)[1]
b.lp <- coef(fit.lgst)[2]
c.lp <- coef(fit.lgst)[3]

# plot fitted logistic function for survival multipliers
plot(K.vec, surv.mult.vec, pch = 19, xlab = "Nf", ylab = "reduction in survival")
K.pred.vec <- seq(curr.pop, K.pop, 10)  # generate sequence of population sizes
pred.surv.mult <- coef(fit.lgst)[1] / (exp((coef(fit.lgst)[2] - K.pred.vec) / coef(fit.lgst)[3]) + 1)
lines(K.pred.vec, pred.surv.mult, lty = 2, lwd = 1, col = "red")

## project
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10
cat.pr <- 0.14 # probability of a catastrophe occurring per generation
yr.now <- 2025  # current year (starting point of the projection)
yr.end <- 2050  # end year for the projection
yrs <- seq(yr.now, yr.end, 1)

t <- (yr.end - yr.now)  # Total number of years to project (e.g., 24 years)

n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
s.arr <- m.arr <- array(data=NA, dim=c(t+1, longev+1, iter))
for (e in 1:iter) {
  popmat <- popmat.orig
  
  n.mat <- s.mat <- m.mat <- matrix(0, nrow=longev+1,ncol=(t+1))
  n.mat[,1] <- init.vec
  s.mat[,1] <- pred.surv
  m.mat[,1] <- pred.fert
  
  for (i in 1:t) {
    # stochastic survivals following beta distribution
    s.alpha <- estBetaParams(pred.surv, pred.surv.se^2)$alpha
    s.beta <- estBetaParams(pred.surv, pred.surv.se^2)$beta
    s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
    
    # catastrophe
    if (rbinom(1, 1, cat.pr/gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      s.stoch <- s.stoch * (rbeta(1, cat.alpha, cat.beta))
    }
    
    # stochastic fertilty sampler (Gaussian)
    fert.stch <- rnorm(length(popmat[,1]), pred.fert, 0.1*pred.fert)
    m.arr[i,,e] <- ifelse(fert.stch < 0, 0, fert.stch)
    
    totN.i <- sum(n.mat[,i], na.rm=T)
    pred.red <- a.lp/(exp((b.lp-totN.i)/c.lp)+1) # density reduction factor
    
    diag(popmat[2:(stages),]) <- (s.stoch[-(stages)]) * pred.red
    popmat[stages,stages] <- (s.stoch[stages]) * pred.red
    popmat[1,] <- m.arr[i,,e]
    n.mat[,i+1] <- popmat %*% n.mat[,i]
    
    s.arr[i,,e] <- s.stoch * pred.red
    
  } # end i loop
  
  n.sums.mat[e,] <- as.vector(colSums(n.mat))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,n.md,type="l", main = "", xlab="year", ylab="N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,n.up,lty=2,col="red",lwd=1.5)
abline(h=K.pop, col="blue", lty=2, lwd=0.9)

s.add <- m.add  <- rep(0, stages)
for (m in 1:iter) {
  s.add <- rbind(s.add, s.arr[ceiling(gen.l):(t+1),,m])
  m.add <- rbind(m.add, m.arr[ceiling(gen.l):(t+1),,m])
  if (m %% itdiv==0) print(m) 
}
s.add <- s.add[-1,]
m.add <- m.add[-1,]

s.md <- apply(s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
s.up <- apply(s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
s.lo <- apply(s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(age.vec,s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(s.lo),1.05*max(s.up)))
lines(age.vec,s.lo,lty=2,col="red",lwd=1.5)
lines(age.vec,s.up,lty=2,col="red",lwd=1.5)

m.md <- apply(m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
m.up <- apply(m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
m.lo <- apply(m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(age.vec,m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(m.lo),1.05*max(m.up)))
lines(age.vec,m.lo,lty=2,col="red",lwd=1.5)
lines(age.vec,m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))

# export for external plotting
stochcatDF.proj.out <- data.frame(year=yrs, N=2*round(n.md,0), N.up=2*round(n.up,0), N.lo=2*round(n.lo,0))
write.csv(stochcatDF.proj.out, file="stochcatDFprojout.csv", row.names=FALSE)


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## stochastic projection with catastrophes and density feedback,
## proportional harvests over 10 years
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## project
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10
cat.pr <- 0.14 # probability of a catastrophe occurring per generation
yr.now <- 2025  # current year (starting point of the projection)
yr.end <- 2035  # end year for the projection
yrs <- seq(yr.now, yr.end, 1)
t <- (yr.end - yr.now)  # Total number of years to project (e.g., 24 years)

## create proportional harvest vector
harv.vec <- seq(0, 0.5, 0.05)
length(harv.vec)

# storage vectors for harvest proportion increment
harv.tot.sum.md <- harv.tot.sum.up <- harv.tot.sum.lo <- rep(NA, length(harv.vec))
n.min.md <- n.min.up <- n.min.lo <- rep(NA, length(harv.vec))
n.min.prop.md <- n.min.prop.up <- n.min.prop.lo <- rep(NA, length(harv.vec))
harv.tot.sum.prop.md <- harv.tot.sum.prop.up <- harv.tot.sum.prop.lo <- rep(NA, length(harv.vec))

for (h in 1:length(harv.vec)) {
  
  n.sums.mat <- harv.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    popmat <- popmat.orig
    
    n.mat <- harv.mat <- matrix(0, nrow=longev+1,ncol=(t+1))
    n.mat[,1] <- init.vec
    
    for (i in 1:t) {
      # stochastic survivals following beta distribution
      s.alpha <- estBetaParams(pred.surv, pred.surv.se^2)$alpha
      s.beta <- estBetaParams(pred.surv, pred.surv.se^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
      
      # catastrophe
      if (rbinom(1, 1, cat.pr/gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        s.stoch <- s.stoch * (rbeta(1, cat.alpha, cat.beta))
      }
      
      # stochastic fertilty sampler (Gaussian)
      fert.stch <- rnorm(length(popmat[,1]), pred.fert, 0.1*pred.fert)
      
      totN.i <- sum(n.mat[,i], na.rm=T)
      pred.red <- a.lp/(exp((b.lp-totN.i)/c.lp)+1) # density reduction factor
      
      diag(popmat[2:(stages),]) <- (s.stoch[-(stages)]) * pred.red
      popmat[stages,stages] <- (s.stoch[stages]) * pred.red
      popmat[1,] <- m.arr[i,,e]
      n.mat[,i+1] <- popmat %*% n.mat[,i]
      
      # remove harvested individuals
      harvested.tot <- round(sum(n.mat[,i+1]) * harv.vec[h], 0)
      stage.dist.i <- n.mat[,i+1]/sum(n.mat[,i+1])
      harvested.i <- round(stage.dist.i * harvested.tot, 0)
      harv.mat[,i+1] <- harvested.i # record harvested individuals
      n.mat[,i+1] <- n.mat[,i+1] - harvested.i
      
    } # end i loop
    
    n.sums.mat[e,] <- as.vector(colSums(n.mat))
    harv.sums.mat[e,] <- as.vector(colSums(harv.mat))
    
    if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  harv.md <- apply(harv.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  harv.up <- apply(harv.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  harv.lo <- apply(harv.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  par(mfrow=c(1,2))
  plot(yrs,n.md,type="l", main = "", xlab="year", ylab="N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
  lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
  lines(yrs,n.up,lty=2,col="red",lwd=1.5)
  abline(h=K.pop, col="blue", lty=2, lwd=0.9)

  plot(yrs,harv.md,type="l", main = "", xlab="year", ylab="N harvested", lwd=2, ylim=c(0.95*min(harv.lo),1.05*max(harv.up)))
  lines(yrs,harv.lo,lty=2,col="red",lwd=1.5)
  lines(yrs,harv.up,lty=2,col="red",lwd=1.5)
  par(mfrow=c(1,2))
  
  ## total harvested over projection period
  harv.tot.sum.md[h] <- sum(harv.md, na.rm=T)
  harv.tot.sum.up[h] <- sum(harv.up, na.rm=T)
  harv.tot.sum.lo[h] <- sum(harv.lo, na.rm=T)
  
  # minimum population size achieved
  n.min.md[h] <- min(n.md, na.rm=T)
  n.min.up[h] <- min(n.up, na.rm=T)
  n.min.lo[h] <- min(n.lo, na.rm=T)
  
  # minimum proportion of N0 achieved
  n.min.prop.md[h] <- min(n.md, na.rm=T)/curr.pop
  n.min.prop.up[h] <- min(n.up, na.rm=T)/curr.pop
  n.min.prop.lo[h] <- min(n.lo, na.rm=T)/curr.pop
  
  # proportion total harvested relative to initial population size
  harv.tot.sum.prop.md[h] <- sum(harv.md, na.rm=T)/curr.pop
  harv.tot.sum.prop.up[h] <- sum(harv.up, na.rm=T)/curr.pop
  harv.tot.sum.prop.lo[h] <- sum(harv.lo, na.rm=T)/curr.pop
  
  Sys.sleep(2) # pause for two seconds
  print("********************************")
  print(paste("harvest proportion = ", harv.vec[h], sep=""))
  print("********************************")
  Sys.sleep(2) # pause for two seconds
  
} # end h loop

## plot relationship between harvest proportion and minimum n (including upper and lower estimates)
harv.min.df <- data.frame(harv.vec, n.min.md, n.min.up, n.min.lo)
plot1 <- ggplot(harv.min.df, aes(x=harv.vec, y=n.min.md)) +
  geom_line() +
  geom_ribbon(aes(ymin=n.min.lo, ymax=n.min.up), alpha=0.2) +
  labs(x="annual harvest proportion", y="minimum population size over 10 years") +
  theme_minimal()

## plot relationship between harvest proportion and total harvested (including upper and lower estimates)
harv.tot.df <- data.frame(harv.vec, harv.tot.sum.md, harv.tot.sum.up, harv.tot.sum.lo)
plot2 <- ggplot(harv.tot.df, aes(x=harv.vec, y=harv.tot.sum.md)) +
  geom_line() +
  geom_ribbon(aes(ymin=harv.tot.sum.lo, ymax=harv.tot.sum.up), alpha=0.2) +
  labs(x="annual harvest proportion", y="total individuals harvested over 10 years") +
  theme_minimal()

## plot relationship between harvest proportion and minimum proportion of N0 (including upper and lower estimates)
harv.min.prop.df <- data.frame(harv.vec, n.min.prop.md, n.min.prop.up, n.min.prop.lo)
plot3 <- ggplot(harv.min.prop.df, aes(x=harv.vec, y=n.min.prop.md)) +
  geom_line() +
  geom_ribbon(aes(ymin=n.min.prop.lo, ymax=n.min.prop.up), alpha=0.2) +
  labs(x="annual harvest proportion", y="minimum population size as proportion of N0") +
  theme_minimal()

## plot relationship between harvest proportion and total harvested as proportion of N0 (including upper and lower estimates)
harv.tot.sum.prop.df <- data.frame(harv.vec, harv.tot.sum.prop.md, harv.tot.sum.prop.up, harv.tot.sum.prop.lo)
plot4 <- ggplot(harv.tot.sum.prop.df, aes(x=harv.vec, y=harv.tot.sum.prop.md)) +
  geom_line() +
  geom_ribbon(aes(ymin=harv.tot.sum.prop.lo, ymax=harv.tot.sum.prop.up), alpha=0.2) +
  labs(x="annual harvest proportion", y="total individuals harvested over 10 years as proportion of N0") +
  theme_minimal()

gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol=2)

## export for external plotting
harvNnProp.proj.out <- data.frame(harv.vec, 2*n.min.md, 2*n.min.up, 2*n.min.lo, n.min.prop.md, n.min.prop.up, n.min.prop.lo)
write.csv(harvNnProp.proj.out, file="harvNnPropprojout.csv", row.names=FALSE)

harvTharv.prop.df <- data.frame(harv.vec, 2*harv.tot.sum.md, 2*harv.tot.sum.up, 2*harv.tot.sum.lo,
                                 harv.tot.sum.prop.md, harv.tot.sum.prop.up, harv.tot.sum.prop.lo)
write.csv(harvTharv.prop.df, file="harvTharvpropout.csv", row.names=FALSE)

