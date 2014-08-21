# R code used in final analysis of data.
# lprq2 function is based on (and nearly identical to) lprq function from quantreg package.

# We define a function for representing an ARMA
# process as an AR process. We assume the ARMA process
# is invertible and use the result in Brockwell and Davis (2002), p. 86.
# In general, the resulting AR process will have infinite
# order; this function returns at least lag.min coefficients.
"ARMAtoAR" <- function(ar, ma, lag.min) {
  p <- 1:1 # this will represent pi_j terms
  p[1] = 1 # p is shifted one to  the right since we can't have p[0]
  for (j in 1:lag.min) {
    if (j > length(ar)) ar[j] = 0
    s = 0
    for (k in 1:length(ma)) {
      if (j - k >= 0) s = s + ma[k]*p[j-k+1]
    }
    p[j+1] = -ar[j] - s
  }
  for (k in 1:length(ma)) {
    for (j in 1:lag.min) {
      ar[j+k-1] = ar[j+k-1] + ma[k]*p[j]
    }
  }
  ar
}

# We define a new function for performing locally linear quantile
# regression that also returns confidence interval bounds.
"lprq2" <- function (x, y, h, tau = 0.5, m = 50) 
{
    xx <- seq(min(x), max(x), length = m)
    fv <- xx
    dv <- xx
    lb <- xx
    ub <- xx
    for (i in 1:length(xx)) {
        z <- x - xx[i]
        wx <- dnorm(z/h)
        r <- rq(y ~ z, weights = wx, tau = tau, ci = FALSE)
        fv[i] <- r$coef[1]
        dv[i] <- r$coef[2]
        ci <- summary.rq(r, se = "rank", alpha = 1 - confLevel)$coefficients
        lb[i] <- ci[3]
        ub[i] <- ci[5]
        ## mcmb <- rqmcmb(z, y, tau = tau)
        ## ci <- rqmcmb.ci(mcmb, alpha = 1 - confLevel)
        ## lb[i] <- ci[1,3]
        ## ub[i] <- ci[1,4]
    }
    list(xx = xx, fv = fv, dv = dv, lb = lb, ub = ub)
}

setEPS()
makeEPS = T # Boolean determining whether to save plots to EPS files.
confLevel = 0.90 # Confidence level to use throughout.

library(quantreg)

# We read the data in from the file.
temp <- read.csv("../data/ohare/ohare_1960_2010_temp_weekly.csv",head=TRUE,sep=",")
timeTemp <- list(vector(),vector(),vector(),vector(),vector())
dec <- list(vector(),vector(),vector(),vector(),vector())
startWeek = -1 # index of first week of year
startYear = -1 # index of first year
curYear = -1 # year of last observation
idx = 1 # decade with which we are dealing
for (i in 1:nrow(temp)) {
  if (startWeek == -1) { # first observation
    startWeek = temp[i,3]
    curYear = startYear = temp[i,2]
  }
  if (temp[i,2] > curYear) { # new year
    startWeek = temp[i,3]
    curYear = temp[i,2]
  }
  if ((temp[i,2] > startYear+9) && !(idx==5)) { # new decade
    idx = idx + 1
    startYear = temp[i,2]
  }
  timeTemp[[idx]] = append(timeTemp[[idx]], temp[i,3] - startWeek + 1)
  dec[[idx]] = append(dec[[idx]], temp[i,5])
}

# We compare the results of using different
# bandwidths.
for (i in 1:1) {
  if (makeEPS) postscript(paste("bandwidths-", i, ".eps")) else dev.new()
  fit1 <- lprq(timeTemp[[i]], dec[[i]],
                   (max(timeTemp[[i]]) - min(timeTemp[[i]])) / 20,
              m = max(timeTemp[[i]]), tau=0.5) # (X_max - X_min) / 20
  fit2 <- lprq(timeTemp[[i]], dec[[i]],
                   (max(timeTemp[[i]]) - min(timeTemp[[i]])) / 10,
              m = max(timeTemp[[i]]), tau=0.5) # (X_max - X_min) / 10
  fit3 <- lprq(timeTemp[[i]], dec[[i]],
                   (max(timeTemp[[i]]) - min(timeTemp[[i]])) / 5,
              m = max(timeTemp[[i]]), tau=0.5)
  fit4 <- lprq(timeTemp[[i]], dec[[i]],
                   1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
              m = max(timeTemp[[i]]), tau=0.5)
  plot(timeTemp[[i]], dec[[i]], main=paste("Bandwidth comparison, decade ", i), ylab="Air Temperature at 1200 UTC", xlab="Week number",
       cex=0.05, col="red")
  lines(fit1$fv[1:max(timeTemp[[i]])], col="blue")
  lines(fit2$fv[1:max(timeTemp[[i]])], col="black")
  lines(fit3$fv[1:max(timeTemp[[i]])], col="green")
  lines(fit4$fv[1:max(timeTemp[[i]])], col="yellow")
  legend(0, 25, c("(X_max - X_min)/20", "(X_max - X_min) / 10",
                  "(X_max - X_min)/5", "Xiao et. al (2003) choice"),
                  lty="solid", col=c("blue", "black", "green", "yellow"))
  if (makeEPS) dev.off()
}
# We look at the residuals obtained using different
# bandwidths.
i = 1 # Similar results are obtained for other decades.
fit <- list(vector(), vector(), vector(), vector(), vector())
resid <- list(vector(), vector(), vector(), vector(), vector())
fit[[i]] <- lprq(timeTemp[[i]], dec[[i]],
                 (max(timeTemp[[i]]) - min(timeTemp[[i]])) / 5,
                 m = max(timeTemp[[i]]), tau=0.5) # (X_max - X_min) / 5
for (j in 1:length(dec[[i]])) {
  resid[[i]] = append(resid[[i]], dec[[i]][j] -
         fit[[i]]$fv[timeTemp[[i]][j]])
}
if (makeEPS) postscript(paste("band-resid-", i, "-r5.eps")) else dev.new()
par(mfrow=c(3,1))
plot(timeTemp[[i]], dec[[i]],main=paste("Residuals for Decade", i, "using bandwidth (X_max - X_min) / 5"), cex=0.05,col="red")
points(resid[[i]], cex=0.05, col="blue")
lines(fit[[i]]$fv[1:max(timeTemp[[i]])], col="black")
acf(resid[[i]], ci=confLevel, main=paste("Residuals for Decade ", i))
pacf(resid[[i]], ci=confLevel, main=paste("Residuals for Decade ", i))
if (makeEPS) dev.off()
fit[[i]] <- lprq(timeTemp[[i]], dec[[i]],
                 1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
                 m = max(timeTemp[[i]]), tau=0.5)
for (j in 1:length(dec[[i]])) {
  resid[[i]] = append(resid[[i]], dec[[i]][j] -
         fit[[i]]$fv[timeTemp[[i]][j]])
}
if (makeEPS) postscript(paste("band-resid-", i, "-xiao.eps")) else dev.new()
par(mfrow=c(3,1))
plot(timeTemp[[i]], dec[[i]],main=paste("Residuals for Decade", i, "using Xiao et al. (2003) bandwidth"), cex=0.05,col="red")
points(resid[[i]], cex=0.05, col="blue")
lines(fit[[i]]$fv[1:max(timeTemp[[i]])], col="black")
acf(resid[[i]], ci=confLevel, main=paste("Residuals for Decade ", i))
pacf(resid[[i]], ci=confLevel, main=paste("Residuals for Decade ", i))
if (makeEPS) dev.off()

# Step 1 from section 2.2 of Xiao et al. (2003): calculate estimated residuals of
# local polynomial regression.
resid <- list(vector(),vector(),vector(),vector(),vector())
fit<-list(1)
for (i in 1:5) {
  fit[[i]] <- lprq(timeTemp[[i]], dec[[i]],
                   1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
              m = max(timeTemp[[i]]), tau=0.5)
  for (j in 1:length(dec[[i]])) {
    resid[[i]] = append(resid[[i]], dec[[i]][j] -
           fit[[i]]$fv[timeTemp[[i]][j]])
  }
}
# We plot the estimated residuals as well as their ACF's and PACF's.
for (i in 1:5) {
  if (makeEPS && i == 1)
    postscript("qr-arma-dec1.eps")
  else
    dev.new()
  par(mfrow=c(3,1))
  plot(timeTemp[[i]], dec[[i]], main=paste("Residuals for Decade ", i),
       cex=0.05,col="red")
  points(resid[[i]], cex=0.05, col="blue")
  lines(fit[[i]]$fv[1:max(timeTemp[[i]])], col="black")
  acf(resid[[i]], ci=confLevel, main=paste("Residuals for Decade ", i))
  pacf(resid[[i]], ci=confLevel, main=paste("Residuals for Decade ", i))
  if (makeEPS && i == 1) dev.off()
}

# Step 2: We fit AR/ARMA models to the residuals and perform
# transforms. For each decade, an AR model seems to fit best.
# In practice, though, it turns out that the transforms do
# not work well and the data seem to be closer to IID
# without them, except in the case of decade 1.
# We only keep the transformed data for decade 1.
for (i in 1:5) {
  a_h <- ar(resid[[i]])
  Y = dec[[i]]
  m_h = rep(fit[[i]]$fv, length.out=length(Y))
  Y_tr <- 1:length(Y)
  for (i2 in 1:length(Y_tr)) {
    s = 0
    for (j in 1:a_h$order) {
      if (!(i2 - j > 0)) break
      s = s + a_h$ar[j]*(Y[i2-j] - m_h[i2-j])
    }
    Y_tr[i2] = Y[i2] - s
  }
  # We look at transformed data.
  if (makeEPS) postscript(paste("data-transformed-", i, ".eps")) else dev.new()
  plot(1:length(Y), rep(c(-25,25),length.out=length(Y)), type="n",
  		    ylab="Air Temperature at 1200 UTC", xlab="Observation number")
  points(1:length(Y), Y, col="black", cex=0.05)
  points(1:length(Y_tr), Y_tr, col="red", cex=0.05)
  legend(0,-20, c("Original data", "Transformed data"),
         col=c("black", "red"), pch=20)
  if (makeEPS) dev.off()
  fit[[i]] <- lprq(timeTemp[[i]], Y_tr,
                   1.06*sd(timeTemp[[i]])*(length(Y_tr)^(-1/5)),
                   m = max(timeTemp[[i]]), tau=0.5)
  resid[[i]] <- vector()
  for (j in 1:length(Y_tr)) {
    resid[[i]] = append(resid[[i]], Y_tr[j] -
           fit[[i]]$fv[timeTemp[[i]][j]])
  }
  if (makeEPS) postscript(paste("data-tr-resid-", i, ".eps")) else dev.new()
  par(mfrow=c(3,1))
  plot(1:length(Y_tr), Y_tr,
       main=paste("Residuals for Decade", i, "after transform"), cex=0.05,
       col="red")
  points(1:length(Y_tr), resid[[i]], cex=0.05, col="blue")
  lines(1:length(Y), rep(fit[[i]]$fv,length.out=length(Y)), col="black")
  acf(resid[[i]], ci=confLevel, main=paste("Residuals for Decade ", i))
  pacf(resid[[i]], ci=confLevel, main=paste("Residuals for Decade ", i))
  if (makeEPS) dev.off()
  if (i == 1) dec[[i]] = Y_tr
}

# Now we compare the results of nonparametric regression
# across decades.
colors <- c("red", "blue", "black", "green", "yellow")
par(mfrow=c(1,1))
# First, the 95th percentile.
fit <- list(1)
for (i in 1:5) {
  fit[[i]] <- lprq2(timeTemp[[i]], dec[[i]],
                    1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
                   m = max(timeTemp[[i]]), tau=0.95)
}
if (makeEPS) postscript("decades-95th.eps") else dev.new()
plot(1:54, rep(c(-10, 25), length.out=54), type="n", main="95th percentile",
     ylab="Air Temperature at 1200 UTC", xlab="Week number")
legend(2, 25, c("Decade 1", "Decade 2", "Decade 3", "Decade 4",
                "Decade 5"), lty="solid", col=colors)
for (i in 1:5) {
  lines(1:max(timeTemp[[i]]), fit[[i]]$fv, col=colors[i])
  lines(1:max(timeTemp[[i]]), fit[[i]]$lb, col=colors[i], lty="dotted")
  lines(1:max(timeTemp[[i]]), fit[[i]]$ub, col=colors[i], lty="dotted")
}
if (makeEPS) dev.off()

# Next, the median.
fit <- list(1)
for (i in 1:5) {
  fit[[i]] <- lprq2(timeTemp[[i]], dec[[i]],
                    1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
                   m = max(timeTemp[[i]]), tau=0.50)
}
if (makeEPS) postscript("decades-median.eps") else dev.new()
plot(1:54, rep(c(-15, 25),length.out=54),
     type="n", main="Median", ylab="Air Temperature at 1200 UTC", xlab="Week number")
legend(2, 25, c("Decade 1", "Decade 2", "Decade 3", "Decade 4",
                "Decade 5"), lty="solid", col=colors)
for (i in 1:5) {
  lines(1:max(timeTemp[[i]]), fit[[i]]$fv, col=colors[i])
  lines(1:max(timeTemp[[i]]), fit[[i]]$lb, col=colors[i], lty="dotted")
  lines(1:max(timeTemp[[i]]), fit[[i]]$ub, col=colors[i], lty="dotted")

}
if (makeEPS) dev.off()

# Finally, the 5th percentile.
fit <- list(1)
for (i in 1:5) {
  fit[[i]] <- lprq2(timeTemp[[i]], dec[[i]],
                    1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
                   m = max(timeTemp[[i]]), tau=0.05)
}
if (makeEPS) postscript("decades-5th.eps") else dev.new()
plot(1:54, rep(c(-25, 25), length.out=54), type="n", main="5th percentile",
     ylab="Air Temperature at 1200 UTC", xlab="Week number")
legend(2, 25, c("Decade 1", "Decade 2", "Decade 3", "Decade 4",
                "Decade 5"), lty="solid", col=colors)
for (i in 1:5) {
  lines(1:max(timeTemp[[i]]), fit[[i]]$fv, col=colors[i])
  lines(1:max(timeTemp[[i]]), fit[[i]]$lb, col=colors[i], lty="dotted")
  lines(1:max(timeTemp[[i]]), fit[[i]]$ub, col=colors[i], lty="dotted")
}
if (makeEPS) dev.off()

# In addition, for each decade, we look at the results for
# the different quantiles.
# Note: There is a downward spike in the lower bound of the
# confidence interval for the 5th percentile for decade 3.
# This spike is not accounted for by outliers in the data.
for (i in 1:5) {
  if (makeEPS) postscript(paste("decade-",i,"detailed.eps")) else dev.new()
  plot(1:max(timeTemp[[i]]), rep(c(-20, 25), length.out=max(timeTemp[[i]])),
                                 type="n", main=paste("Decade",i),
       yaxp=c(-20, 25, 9))
  legend(2, 25, c("95th percentile", "Median", "5th percentile"),
         lty="solid", col=c("red", "black", "blue"))
  fit[[i]] <- lprq2(timeTemp[[i]], dec[[i]],
                    1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
                   m = max(timeTemp[[i]]), tau=0.05)
  lines(1:max(timeTemp[[i]]), fit[[i]]$fv, col="blue")
  lines(1:max(timeTemp[[i]]), fit[[i]]$lb, col="blue", lty="dotted")
  lines(1:max(timeTemp[[i]]), fit[[i]]$ub, col="blue", lty="dotted")
  fit[[i]] <- lprq2(timeTemp[[i]], dec[[i]],
                    1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
                   m = max(timeTemp[[i]]), tau=0.50)
  lines(1:max(timeTemp[[i]]), fit[[i]]$fv, col="black")
  lines(1:max(timeTemp[[i]]), fit[[i]]$lb, col="black", lty="dotted")
  lines(1:max(timeTemp[[i]]), fit[[i]]$ub, col="black", lty="dotted")
  fit[[i]] <- lprq2(timeTemp[[i]], dec[[i]],
                    1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
                   m = max(timeTemp[[i]]), tau=0.95)
  lines(1:max(timeTemp[[i]]), fit[[i]]$fv, col="red")
  lines(1:max(timeTemp[[i]]), fit[[i]]$lb, col="red", lty="dotted")
  lines(1:max(timeTemp[[i]]), fit[[i]]$ub, col="red", lty="dotted")
  if (makeEPS) dev.off()
}

# Finally, we make pairwise comparisons across decades,
# comparing decades 1 and 4, 1 and 5, 2 and 4, and
# 2 and 5. We look at the 95th and 5th percentiles
# and the median.
"pairwiseCmp" <- function(i1, i2, fit_95, fit_med, fit_5) {
  if (makeEPS) postscript(paste("Decade comparison-", i1, "and", i2, ".eps")) else dev.new()
  plot(1:max(c(timeTemp[[i1]], timeTemp[[i2]])),
       rep(c(-20, 25), length.out=max(c(timeTemp[[i1]], timeTemp[[i2]]))),
       type="n", main=paste("Decade",i1,"versus decade",i2,"at the 95th, 50th, and 5th percentiles"), ylab="Air Temperature at 1200 UTC", xlab="Week number",
       yaxp=c(-20, 25, 9))
  legend(2, 25, c(paste("Decade",i1),paste("Decade",i2)),
         lty="solid", col=c("black", "red"))
  lines(1:max(timeTemp[[i1]]), fit_95[[i1]]$fv, col="black")
  lines(1:max(timeTemp[[i1]]), fit_95[[i1]]$lb, col="black", lty="dotted")
  lines(1:max(timeTemp[[i1]]), fit_95[[i1]]$ub, col="black", lty="dotted")
  lines(1:max(timeTemp[[i1]]), fit_med[[i1]]$fv, col="black")
  lines(1:max(timeTemp[[i1]]), fit_med[[i1]]$lb, col="black", lty="dotted")
  lines(1:max(timeTemp[[i1]]), fit_med[[i1]]$ub, col="black", lty="dotted")
  lines(1:max(timeTemp[[i1]]), fit_5[[i1]]$fv, col="black")
  lines(1:max(timeTemp[[i1]]), fit_5[[i1]]$lb, col="black", lty="dotted")
  lines(1:max(timeTemp[[i1]]), fit_5[[i1]]$ub, col="black", lty="dotted")
  lines(1:max(timeTemp[[i2]]), fit_95[[i2]]$fv, col="red")
  lines(1:max(timeTemp[[i2]]), fit_95[[i2]]$lb, col="red", lty="dotted")
  lines(1:max(timeTemp[[i2]]), fit_95[[i2]]$ub, col="red", lty="dotted")
  lines(1:max(timeTemp[[i2]]), fit_med[[i2]]$fv, col="red")
  lines(1:max(timeTemp[[i2]]), fit_med[[i2]]$lb, col="red", lty="dotted")
  lines(1:max(timeTemp[[i2]]), fit_med[[i2]]$ub, col="red", lty="dotted")
  lines(1:max(timeTemp[[i2]]), fit_5[[i2]]$fv, col="red")
  lines(1:max(timeTemp[[i2]]), fit_5[[i2]]$lb, col="red", lty="dotted")
  lines(1:max(timeTemp[[i2]]), fit_5[[i2]]$ub, col="red", lty="dotted")
}
fit_95 <- list(vector(),vector(),vector(),vector(),vector())
fit_med <- list(vector(),vector(),vector(),vector(),vector())
fit_5 <- list(vector(),vector(),vector(),vector(),vector())
for (i in 1:5) {
  fit_95[[i]] <- lprq2(timeTemp[[i]], dec[[i]],
                    1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
                    m = max(timeTemp[[i]]), tau=0.95)
  fit_med[[i]] <- lprq2(timeTemp[[i]], dec[[i]],
                    1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
                    m = max(timeTemp[[i]]), tau=0.50)
  fit_5[[i]] <- lprq2(timeTemp[[i]], dec[[i]],
                    1.06*sd(timeTemp[[i]])*(length(dec[[i]])^(-1/5)),
                    m = max(timeTemp[[i]]), tau=0.05)
}
pairwiseCmp(1, 4, fit_95, fit_med, fit_5)
pairwiseCmp(1, 5, fit_95, fit_med, fit_5)
pairwiseCmp(2, 4, fit_95, fit_med, fit_5)
pairwiseCmp(2, 5, fit_95, fit_med, fit_5)
