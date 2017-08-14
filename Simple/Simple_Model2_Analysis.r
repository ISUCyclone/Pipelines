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

# Page 10, Figure 5(a)
#******************************
mubeta1 <- pmean[1,5]
sigmabeta <- pmean[5,5]
sigma <- pmean[6,5]

day0 <- 7347
day.seq <- seq(from = day0, to = day0 + 365*30, by = 10)

post2 <- as.array(fit2)
iter2 <- dim(post2)[1]
post2 <- rbind(post2[,1,],post2[,2,],post2[,3,],post2[,4,]) # Merge 4 chains

samplesize2 = 1000
index <- sample(1:(4*iter2), samplesize2)
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

# Page 10, Figure 5(b)
#******************************
tc <- 12072
t.seq.new <- seq(from = tc, to = tc + 365*20, by = 10)

remaining.life.time <- matrix(nrow = length(index), ncol = length(t.seq.new))
failure.prob.tc <- NULL
failure.prob <- matrix(nrow = length(index), ncol = length(t.seq.new))

for(i in 1:length(index)){
  mcmc.mu.beta <- post2.sub[i, 1]
  mcmc.sigma.beta <- post2.sub[i, 5]
  failure.prob.tc <- pnorm((log((tc-7347)/365) - log(0.2) + mcmc.mu.beta)/mcmc.sigma.beta)
  failure.prob[i,] <- pnorm((log((t.seq.new-7347)/365)-log(0.2)+mcmc.mu.beta)/mcmc.sigma.beta)
  remaining.life.time[i,] <- (failure.prob[i,] - failure.prob.tc)/(1-failure.prob.tc)
}

remaining.life.time.mean <- colMeans(remaining.life.time, na.rm = FALSE, dims = 1)
remaining.life.time.lower <- apply(remaining.life.time, 2, quantile, probs = 0.025)
remaining.life.time.upper <- apply(remaining.life.time, 2, quantile, probs = 0.975)
remaining.life.time.lower10 <- apply(remaining.life.time, 2, quantile, probs = 0.1)
remaining.life.time.upper10 <- apply(remaining.life.time, 2, quantile, probs = 0.9)

xc <- (t.seq.new -tc)/365
remaining.life.time.mean.lognormal <- qnorm(remaining.life.time.mean)
remaining.life.time.lower.lognormal <- qnorm(remaining.life.time.lower)
remaining.life.time.upper.lognormal <- qnorm(remaining.life.time.upper)
remaining.life.time.lower10.lognormal <- qnorm(remaining.life.time.lower10)
remaining.life.time.upper10.lognormal <- qnorm(remaining.life.time.upper10)

SplidaOptions(SPLIDA.DateOnPlot = FALSE)
probpaper("Lognormal", x.range = c(0.05, 20), title.option = "blank", grid = FALSE,
          cex = 1.4, cex.labs = 1.4, cex.tic.lab = 1.4, 
          y.range = c(1/100000, 0.2), xlab = "Years Since January 2003",
          ylab = "Fraction Failing")

lines(log(xc),remaining.life.time.mean.lognormal,col="blue",lwd=2.5,lty=1)
lines(log(xc), remaining.life.time.lower.lognormal, lty=2, col="deeppink", lwd=2.5)
lines(log(xc), remaining.life.time.upper.lognormal, lty=2, col="deeppink", lwd=2.5)
lines(log(xc), remaining.life.time.lower10.lognormal, lty=3, col="green", lwd=2.5)
lines(log(xc), remaining.life.time.upper10.lognormal, lty=3, col="green", lwd=2.5)

legend("topleft", legend = c("cdf Estimate","95% Credible Intervals","80% Credible Intervals" ), col=c("blue","deeppink","green"), 
       lty=c(1,2,3), merge=TRUE , border=c(0,0,1,0), lwd=c(2.5, 2.5, 2.5))

# Page 12, Figure 6
tc <- 12072
t.seq.new <- seq(from=tc, to=tc+365*10,by=5)

p <- c(0.1,0.2,0.3,0.4)
M <- 100
p.star <- 1-(1-p)^(1/M)

Ftc <- numeric(length = length(index))
tt <- matrix(nrow = 4,ncol = length(index))
for(j in 1:4){
  for(i in 1:length(index)){
    mcmc.mu <- post2.sub[i, 1]
    mcmc.sigma <- post2.sub[i, 5]
    Ftc[i] <- pnorm((log((tc - t0)/365) - log(0.2) + mcmc.mu)/mcmc.sigma)
    tt[j,i] <- 365 * exp(mcmc.sigma * qnorm(p.star[j]*(1-Ftc[i]) + Ftc[i]) - mcmc.mu + log(0.2)) + t0
  }
}

tt.year <- (tt - tc)/365

density.model2 <- apply(tt.year,1,density)

post1 <- as.array(fit1)
iter1 <- dim(post1)[1]
post1 <- rbind(post1[,1,],post1[,2,],post1[,3,],post1[,4,])

samplesize1 <- 1000
post1.sub <- post1[sample(4*iter1, samplesize1),]

Ftc <- numeric(length = samplesize1)
tt1 <- matrix(nrow = 4,ncol = samplesize1)
for(j in 1:4){
  for(i in 1:samplesize1){
    mcmc.mu <- post1.sub[i, 1]
    mcmc.sigma <- post1.sub[i, 2]
    Ftc[i] <- pnorm((log((tc - t0)/365) - log(0.2) + mcmc.mu)/mcmc.sigma)
    tt1[j,i] <- 365 * exp(mcmc.sigma * qnorm(p.star[j]*(1-Ftc[i]) + Ftc[i]) - mcmc.mu + log(0.2)) + t0
  }
}



tt.year <- (tt - tc)/365
tt1.year <- (tt1 - tc)/365

density.model2 <- apply(tt.year,1,density)
density.model1 <- apply(tt1.year,1,density)

par(las=1)
plot(density.model2[[1]]$x,density.model2[[1]]$y*1,ylim=c(-4, 10),xlim=c(0, 6),type='n', 
     ylab="Quantiles of the Posterior Minimum Remaining Lifetime Distribution", 
     xlab="Years Since January 2003", yaxt="n")

axis(2,c(-2,2,5,8), c(0.1, 0.2, 0.3, 0.4),las=1,cex.axis=0.75)

lines(density.model2[[1]]$x,density.model2[[1]]$y-2, lty=2, col="blue", lwd=2.5)
lines(density.model2[[2]]$x,density.model2[[2]]$y+2, lty=2, col="blue", lwd=2.5)
lines(density.model2[[3]]$x,density.model2[[3]]$y+5, lty=2, col="blue", lwd=2.5)
lines(density.model2[[4]]$x,density.model2[[4]]$y+8, lty=2, col="blue", lwd=2.5)

lines(density.model1[[1]]$x,density.model1[[1]]$y-2, lty=4, col="red", lwd=2.5)
lines(density.model1[[2]]$x,density.model1[[2]]$y+2, lty=4, col="red", lwd=2.5)
lines(density.model1[[3]]$x,density.model1[[3]]$y+5, lty=4, col="red", lwd=2.5)
lines(density.model1[[4]]$x,density.model1[[4]]$y+8, lty=4, col="red", lwd=2.5)

abline(h=-2.1, col="grey")
abline(h=1.9, col="grey")
abline(h=4.9, col="grey")
abline(h=7.9, col="grey")

legend("bottomright", legend = c("Model (1)","Model (2)"), 
       col=c("red","blue"), 
       lty=c(4,2), merge=TRUE ,lwd=c(2.5, 2.5),cex=1.)

