# Page 9, Figure 4
#******************************
pmean <- get_posterior_mean(fit2)
slope <- pmean[7:94,5]
dates <- date::as.date(c("12Feb90","11Feb95","10Feb2000","11May2001",
                         "12Jul2001","18Oct2001","2Apr2002","20Jan2003",
                         "1Jan2004","1Jan2006","1Jan2007","1Jan2008",
                         "1Jan2009","1Jan2010","1Jan2011","1Jan2012",
                         "1Jan2013"),order="dmy")
days <- julian(as.Date(dates,"%d%b%Y"))

thick <- matrix(nrow = 88,ncol = 17)
for(i in 1:88){
  thick[i,] <- 0.25 - slope[i] * (days - days[1])/365
}

x.max <- max(days) + 100
x.min <- min(days) - 600

display.days <- days[-c(4,5,7)]
display.dates <- dates[-c(4,5,7)]

par(las = 1,mgp = c(4,1,0), mar = c(4,7,4,2), mfrow = c(2,2))
plot(days,thick[1,],type = "n", xlim = c(x.min,x.max),
     ylim = c(0, 0.26), xlab = "Date", ylab = "Fitted Wall Thickness",
     xaxt = "n")
axis(1,display.days,display.dates,las = 2, cex.axis = 0.75)

for(i in 1:22){
  lines(days,thick[4*i-3,], type = "l", lwd = 1, lty = 1, col = "red")
}
abline(h = 0.05, lwd = 1, col = "blue")
abline(v = 12072, lwd = 0.5, col = "grey",lty = 2)
abline(v = 9172, lwd = 0.5, col = "grey",lty = 2)
mtext("Quadrant 1 (top)", side = 3, line = 1, cex = 1.2)
mtext("Date", side = 1, line = 3.8, cex = 0.8)

plot(days,thick[1,],type = "n", xlim = c(x.min,x.max),
     ylim = c(0, 0.26), xlab = "Date", ylab = "Fitted Wall Thickness",
     xaxt = "n")
axis(1,display.days,display.dates,las = 2, cex.axis = 0.75)

for(i in 1:22){
  lines(days,thick[4*i-2,], type = "l", lwd = 1, lty = 1, col = "blue")
}
abline(h = 0.05, lwd = 1, col = "blue")
abline(v = 12072, lwd = 0.5, col = "grey",lty = 2)
abline(v = 9172, lwd = 0.5, col = "grey",lty = 2)
mtext("Quadrant 2 (right)", side = 3, line = 1, cex = 1.2)
mtext("Date", side = 1, line = 3.8, cex = 0.8)

plot(days,thick[1,],type = "n", xlim = c(x.min,x.max),
     ylim = c(0, 0.26), xlab = "Date", ylab = "Fitted Wall Thickness",
     xaxt = "n")
axis(1,display.days,display.dates,las = 2, cex.axis = 0.75)

for(i in 1:22){
  lines(days,thick[4*i-1,], type = "l", lwd = 1, lty = 1, col = "green")
}
abline(h = 0.05, lwd = 1, col = "blue")
abline(v = 12072, lwd = 0.5, col = "grey",lty = 2)
abline(v = 9172, lwd = 0.5, col = "grey",lty = 2)
mtext("Quadrant 3 (bottom)", side = 3, line = 1, cex = 1.2)
mtext("Date", side = 1, line = 3.8, cex = 0.8)

plot(days,thick[1,],type = "n", xlim = c(x.min,x.max),
     ylim = c(0, 0.26), xlab = "Date", ylab = "Fitted Wall Thickness",
     xaxt = "n")
axis(1,display.days,display.dates,las = 2, cex.axis = 0.75)

for(i in 1:22){
  lines(days,thick[4*i,], type = "l", lwd = 1, lty = 1, col = "deeppink")
}
abline(h = 0.05, lwd = 1, col = "blue")
abline(v = 12072, lwd = 0.5, col = "grey",lty = 2)
abline(v = 9172, lwd = 0.5, col = "grey",lty = 2)
mtext("Quadrant 4 (left)", side = 3, line = 1, cex = 1.2)
mtext("Date", side = 1, line = 3.8, cex = 0.8)

# Page 10, Figure 5
#******************************
mubeta1 <- pmean[1,5]
sigmabeta <- pmean[5,5]
sigma <- pmean[6,5]

day0 <- 7347
day.seq <- seq(from = day0, to = day0 + 365*30, by = 10)

post2 <- as.array(fit2)
iter <- dim(post2)[1]
post2 <- rbind(post2[,1,],post2[,2,],post2[,3,],post2[,4,]) # Merge 4 chains

samplesize = 1000
index <- sample(1:(4*iter), samplesize)
post2.sub <- post2[index,]

failure.prob <- matrix(nrow = length(index), ncol = length(day.seq))

for(i in 1:length(index)){
  mcmc.mu.beta <- post2.sub[i, 1]
  mcmc.sigma.beta <- post2.sub[i, 5]
  failure.prob[i,] <- pnorm((log((day.seq - 7347)/365) - log(0.2) + mcmc.mu.beta)/mcmc.sigma.beta)
}

failure.prob.mean <- colMeans(failure.prob)
failure.prob.lower <- apply(failure.prob, 2, quantile, probs = 0.025)
failure.prob.upper <- apply(failure.prob, 2, quantile, probs = 0.975)
failure.prob.lower10 <- apply(failure.prob, 2, quantile, probs = 0.1)
failure.prob.upper10 <- apply(failure.prob, 2, quantile, probs = 0.9)

failure.prob.mean.lognormal <- qnorm(failure.prob.mean)
failure.prob.lower.lognormal <- qnorm(failure.prob.lower)
failure.prob.upper.lognormal <- qnorm(failure.prob.upper)
failure.prob.lower10.lognormal <- qnorm(failure.prob.lower10)
failure.prob.upper10.lognormal <- qnorm(failure.prob.upper10)

SplidaOptions(SPLIDA.DateOnPlot)
probpaper("Lognormal", x.range = c(4,30), title.option = "blank",
          grid = FALSE, cex = 1.4, cex.labs = 1.4, cex.tic.lab = 1.4,
          y.range = c(1/10000, 0.2),
          xlab = "Years Since Pipeline Installation (Feb 1990)",
          ylab = "Fraction Failing")
lines(log(x), failure.prob.mean.lognormal, col = "blue",
      lwd = 2.5, lty = 1)
lines(log(x), failure.prob.lower.lognormal, col = "deeppink",
      lwd = 2.5, lty = 2)
lines(log(x), failure.prob.upper.lognormal, col = "deeppink",
      lwd = 2.5, lty = 2)
lines(log(x), failure.prob.lower10.lognormal, col = "green",
      lwd = 2.5, lty = 3)
lines(log(x), failure.prob.upper10.lognormal, col = "green",
      lwd = 2.5, lty = 3)

legend("topleft",
       legend = c("Cdf Estimate","95% Credible Intervals","80% Credible Intervals"),
       col = c("blue","deeppink","green"), lty = c(1,2,3),
       merge = TRUE, border = c(0,0,1,0), lwd = c(2.5, 2.5, 2.5))
