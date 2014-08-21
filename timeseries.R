# R code used in initial time series analysis of data.

setEPS()

makeEPS = T # Boolean determining whether to save EPS files with plots.
obsYear = 365 # Number of observations per year; adjust if using weekly or monthly data.
confLevel = 0.90 # Confidence level for computing confidence intervals.

library(quantreg)

temp <- read.csv("data/ohare/ohare_1960_2010_temp.csv",head=TRUE,sep=",")
temp=ts(temp[,4],start=1,frequency=1)
# For weekly or monthly data, use:
# temp=ts(temp[,5],start=1,frequency=1)

# To remove trend and seasonality, we first try the method
# of section 1.5.2 in Brockwell and Davis (2002).
if (obsYear %% 2 == 1) {
  ma1 = filter(temp, sides = 2, rep(1, obsYear) / obsYear)
} else {
  ma1 = filter(temp, sides = 2, c(0.5, rep(1, obsYear - 1), 0.5) / 52)
}

w <- rep(0, obsYear)
for (k in 1:floor(obsYear / 2)) {
  n = 0
  nx = 0
  j = 1
  while ((k + j * obsYear > floor(obsYear / 2)) && (k + j * obsYear <= length(temp) - floor(obsYear / 2))) {
    n = n + 1
    nx = nx + temp[k+j*obsYear] - ma1[k+j*obsYear]
    j = j + 1
  }
  w[k] = nx / n
}
for (k in (floor(obsYear / 2) + 1):obsYear) {
  n = 0
  nx = 0
  j = 0
  while ((k + j * obsYear > floor(obsYear / 2)) && (k + j * obsYear <= length(temp) - floor(obsYear / 2))) {
    n = n + 1
    nx = nx + temp[k + j * obsYear] - ma1[k + j * obsYear]
    j = j + 1
  }
  w[k] = nx / n
}

s_h <- rep(0,obsYear)
for (k in 1:obsYear) {
  s_h[k] = w[k] - (1 / obsYear)*sum(w)
}
for (k in (obsYear+1):length(temp)) {
  s_h[k] = s_h[k - obsYear]
}

d_t = temp - s_h

fit = rq(d_t~time(d_t), na.action=NULL, tau=0.5)
trend = fitted(fit)

# Finally, we arrive at the estimated noise series.
met1 = temp - trend - s_h


# Alternatively, we can use differencing to
# eliminate trend and seasonality.
lagged = diff(temp, lag=obsYear)
met2 = diff(lagged)

# We also try combining differencing and kernel
# smoothing.
met3 = ksmooth(time(temp), temp, "normal", bandwidth=obsYear*2)[2]
met3 = unlist(met3,use.names=FALSE)
met3 = diff(met3,differences=3)

# Finally, we try periodic regression, as described in Shumway and Stoffer (2011).
wk = time(temp)/365
cs = cos(2*pi*wk)
sn = sin(2*pi*wk)
fit = lm(temp~wk+cs+sn,na.action=NULL)
met4 = resid(fit)

# We compare the results for the first five years of data of
# removing trend and seasonality using each of
# the four methods side by side.
if (makeEPS) postscript("Detrending comparison.eps") else dev.new()
par(mfrow=c(2,1))
plot(temp[1:(obsYear*5)], main="Original data", ylab="Air Temperature at 1200 UTC", xlab="Observation number", xaxp=c(0,obsYear*5,5), cex=0.05)
#plot(met1[1:(obsYear*5)], main="Detrending with regression", ylab="Air Temperature at 1200 UTC", xlab="Observation number", xaxp=c(0,obsYear*5,5), cex=0.05)
#plot(met2[1:(obsYear*5)], main="Detrending with differencing", ylab="Air Temperature at 1200 UTC", xlab="Observation number", xaxp=c(0,obsYear*5,5), cex=0.05)
#plot(met3[1:(obsYear*5)], main="Detrending with differencing and kernel smoothing", ylab="Air Temperature at 1200 UTC", xlab="Observation number", xaxp=c(0,obsYear*5,5), cex=0.05)
plot(met4[1:(obsYear*5)], main="Detrending with periodic regression", ylab="Air Temperature at 1200 UTC", xlab="Observation number", xaxp=c(0,obsYear*5,5), cex=0.05)
if (makeEPS) dev.off()

# Plot ACF comparisons.
if (makeEPS) postscript("acfPacf.eps") else dev.new()
par(mfrow=c(3,1))
acf(met1, obsYear, main="Detrending with regression", ci=confLevel)
pacf(met2, obsYear*3, main="Detrending with differencing", xaxp=c(0, obsYear*3, 3), ci=confLevel)
pacf(met3, obsYear*3, main="Detrending with differencing and kernel smoothing", xaxp=c(0, obsYear*3, 3), ci=confLevel)
if (makeEPS) dev.off()
if (makeEPS) postscript("ACF comparison 1.eps") else dev.new()
par(mfrow=c(2,1))
acf(temp, obsYear*3, main="Original data", xaxp=c(0, obsYear*5, 5), ci=0.90)
acf(met1, obsYear*3, main="Detrending with regression", xaxp=c(0, obsYear*5, 5), ci=confLevel)
if (makeEPS) dev.off()
if (makeEPS) postscript("ACF comparison 2.eps") else dev.new()
par(mfrow=c(2,1))
acf(temp, obsYear*3, main="Original data", xaxp=c(0, obsYear*5, 5), ci=confLevel)
acf(met2, obsYear*3, main="Detrending with differencing", xaxp=c(0, obsYear*5, 5), ci=confLevel)
if (makeEPS) dev.off()
if (makeEPS) postscript("ACF comparison 3.eps") else dev.new()
par(mfrow=c(2,1))
acf(temp, obsYear*3, main="Original data", xaxp=c(0, obsYear*5, 5), ci=confLevel)
acf(met3, obsYear*5, main="Detrending with differencing and kernel smoothing", xaxp=c(0, obsYear*5, 5), ci=confLevel)
if (makeEPS) dev.off()
if (makeEPS) postscript("ACF comparison 4.eps") else dev.new()
par(mfrow=c(2,1))
acf(temp, obsYear*3, main="Original data", xaxp=c(0, obsYear*5, 5), ci=confLevel)
acf(met4, obsYear*5, main="Detrending with periodic regression", xaxp=c(0, obsYear*5, 5), ci=confLevel)
if (makeEPS) dev.off()
# Plot PACF comparisons.
if (makeEPS) postscript("PACF comparison 1.eps") else dev.new()
par(mfrow=c(2,1))
pacf(temp, obsYear*5, main="Original data", xaxp=c(0, obsYear*5, 5), ci=confLevel)
pacf(met1, obsYear*5, main="Detrending with regression", xaxp=c(0, obsYear*5, 5), ci=confLevel)
if (makeEPS) dev.off()
if (makeEPS) postscript("PACF comparison 2.eps") else dev.new()
par(mfrow=c(2,1))
pacf(temp, obsYear*5, main="Original data", xaxp=c(0, obsYear*5, 5), ci=confLevel)
pacf(met2, obsYear*5, main="Detrending with differencing", xaxp=c(0, obsYear*5, 5), ci=confLevel)
if (makeEPS) dev.off()
if (makeEPS) postscript("PACF comparison 3.eps") else dev.new()
par(mfrow=c(2,1))
pacf(temp, obsYear*5, main="Original data", xaxp=c(0, obsYear*5, 5), ci=confLevel)
pacf(met3, obsYear*5, main="Detrending with differencing and kernel smoothing", xaxp=c(0, obsYear*5, 5), ci=confLevel)
if (makeEPS) dev.off()
if (makeEPS) postscript("PACF comparison 4.eps") else dev.new()
par(mfrow=c(2,1))
pacf(temp, obsYear*5, main="Original data", xaxp=c(0, obsYear*5, 5), ci=confLevel)
pacf(met4, obsYear*5, main="Detrending with periodic regression", xaxp=c(0, obsYear*5, 5), ci=confLevel)
if (makeEPS) dev.off()
