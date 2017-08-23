post3 <- as.array(fit3)
post3 <- rbind(post3[,1,],post3[,2,],post3[,3,],post3[,4,])

beta <- post3[,43:75]
init <- post3[,76:108]
thickness.init <- post3[,10:42]

pmean <- get_posterior_mean(fit3)[,5]
beta.mean <- pmean[43:75]
init.mean <- pmean[76:108]
thickness.init.mean <- pmean[10:42]

data <- post3[,c(1:9,109)]

index <- sample(32000, 1000)
data <- data[index,]

Y=structure(.Data=c(
  0.45,0.45,0.45,0.37,0.44,0.43,0.43,0.37,0.45,0.45,0.45,0.4,
  0.43,0.42,0.43,0.43,0.42,0.41,0.41,0.35,0.41,0.42,0.42,0.42,
  0.46,0.46,0.46,0.43,0.43,0.42,0.43,0.42,0.45,0.45,0.45,0.44,
  0.45,0.44,0.45,0.44,0.44,0.43,0.44,0.44,0.44,0.44,0.44,0.44,
  0.43,0.43,0.43,0.43,0.42,0.43,0.42,0.42,0.43,0.44,0.44,0.43,
  0.45,0.44,0.45,0.44,0.45,0.44,0.45,0.41,0.42,0.42,0.42,0.39,
  0.44,0.44,0.44,0.41,0.42,0.42,0.42,0.42,0.41,0.41,0.41,0.37,
  0.42,0.42,0.42,0.39,0.43,0.44,0.43,0.43,0.44,0.44,0.44,0.42,
  0.43,0.44,0.45,0.45,0.43,0.43,0.43,0.42,0.44,0.43,0.43,0.42,
  0.44,0.43,0.43,0.4,0.51,0.51,0.52,0.47,0.51,0.51,0.52,0.49,
  0.55,0.55,0.56,0.55,0.51,0.5,0.51,0.51,0.52,0.52,0.52,0.52),
  .Dim=c(33, 4))
thickness.fitted.measure <- matrix(Y, 33,4, byrow=T)

# Soft failure
Df.elbow <- mean(Y[1:12,1])*.2; Df.elbow 
Df.pipe <- mean(Y[13:28,1])*.2; Df.pipe
Df.tee <- mean(Y[29:33,1])*.2; Df.tee

library(date)
date.tmpt <- as.date(c("1apr2000","1nov2000","1jul2001","1jan2004", "1jan2025"),order="dmy")
JulianDate <- julian(as.Date(date.tmpt, "%d%b%Y"))

t.lower <- min(JulianDate)
t.upper <- max(JulianDate)

#t.seq <- signif(seq(t.lower, t.upper, by=365),5)
t.seq <- signif(seq(t.lower, t.upper, length=100),5)

failure.prob.elbow <- matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq))
failure.prob.pipe <-  matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq))
failure.prob.tee <- matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq))

N <-1000
set.seed(198411)

for (i in 1: dim(data)[1]) {
  
  mcmc.run <- data[i,]
  
  ##############
  # Step 1
  ##############
  
  # For each type of the TML (elbow, pipe and tee), generate N=100,000
  # y0.elbows, y0.pipes, y0.tees,  beta_1s and T_Is.
  y0.elbow <- rnorm(N, mean = mcmc.run[1], sd =  mcmc.run[4])
  y0.pipe <- rnorm(N, mean = mcmc.run[2], sd =  mcmc.run[4])
  y0.tee <- rnorm(N, mean =  mcmc.run[3],  sd =  mcmc.run[4])
  beta1 <- exp(rnorm(N, mean = mcmc.run[10], sd = mcmc.run[6]))
  TI <- exp(rnorm(N, mean = mcmc.run[7], sd = mcmc.run[8]))
  
  for (j in 1: length(t.seq)){
    t <- rep(t.seq[j],N)
    failure.prob.elbow[i,j] <- sum(y0.elbow-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.elbow )/N
    failure.prob.pipe[i,j] <- sum(y0.pipe-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.pipe )/N
    failure.prob.tee[i,j] <- sum(y0.tee-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.tee )/N
  }
}

write.table(failure.prob.elbow, file = "failure.prob.elbow.lognormal.csv", sep=",")
write.table(failure.prob.pipe, file = "failure.prob.pipe.lognormal.csv", sep=",")
write.table(failure.prob.tee, file = "failure.prob.tee.lognormal.csv", sep=",")

failure.prob.elbow <- read.table(file="failure.prob.elbow.lognormal.csv",sep=",")
failure.prob.pipe <- read.table(file="failure.prob.pipe.lognormal.csv",sep=",")
failure.prob.tee <- read.table(file="failure.prob.tee.lognormal.csv",sep=",")

failure.prob.elbow.mean <- colMeans(failure.prob.elbow, na.rm = FALSE, dims = 1)
failure.prob.pipe.mean <- colMeans(failure.prob.pipe, na.rm = FALSE, dims = 1)
failure.prob.tee.mean <- colMeans(failure.prob.tee, na.rm = FALSE, dims = 1)

# Figure 14, (a)
#****************************
failure.prob.elbow.lower<- apply(failure.prob.elbow, 2, quantile, probs=0.025)
failure.prob.elbow.upper<- apply(failure.prob.elbow, 2, quantile, probs=0.975)
failure.prob.elbow.lower10<- apply(failure.prob.elbow, 2, quantile, probs=0.1)
failure.prob.elbow.upper10<- apply(failure.prob.elbow, 2, quantile, probs=0.9)

par(las=1)
plot((t.seq-min(t.seq))/365, failure.prob.elbow.mean, type="n", xlab="Years Since Installation", 
     ylab="Probability", xlim=c(min((t.seq-min(t.seq))/365),25),
     ylim=c(0.0, 0.2))

lines((t.seq-min(t.seq))/365, failure.prob.elbow.mean, lty=1, col="blue", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.elbow.lower, lty=2, col="deeppink", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.elbow.upper, lty=2, col="deeppink", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.elbow.lower10, lty=4, col="green", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.elbow.upper10, lty=4, col="green", lwd=3.5)

legend("bottomright", legend = c("cdf Estimate","95% Credible Intervals","80% Credible Intervals" ), col=c("blue","deeppink","green"), 
lty=c(1,2,4), merge=TRUE , border=c(0,0,1,0), lwd=c(2.5, 2.5, 2.5),cex=0.7)

mtext("TML Component: Elbow", line=1)

# Figure 14, (b)
#****************************
failure.prob.pipe.lower<- apply(failure.prob.pipe, 2, quantile, probs=0.025)
failure.prob.pipe.upper<- apply(failure.prob.pipe, 2, quantile, probs=0.975)
failure.prob.pipe.lower10<- apply(failure.prob.pipe, 2, quantile, probs=0.1)
failure.prob.pipe.upper10<- apply(failure.prob.pipe, 2, quantile, probs=0.9)

plot((t.seq-min(t.seq))/365, failure.prob.pipe.mean, type="n", xlab="Years Since Installation", 
     ylab="Probability", xlim=c(min((t.seq-min(t.seq))/365),25),
     ylim=c(0.0, 0.2))

lines((t.seq-min(t.seq))/365, failure.prob.pipe.mean, lty=1, col="blue", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.pipe.lower, lty=2, col="deeppink", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.pipe.upper, lty=2, col="deeppink", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.pipe.lower10, lty=4, col="green", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.pipe.upper10, lty=4, col="green", lwd=3.5)

legend("bottomright", legend = c("cdf Estimate","95% Credible Intervals","80% Credible Intervals" ), col=c("blue","deeppink","green"), 
lty=c(1,2,4), merge=TRUE , border=c(0,0,1,0), lwd=c(2.5, 2.5, 2.5),cex=0.7)

mtext("TML Component: Pipe", line=1)

# Figure 14, (c)
#****************************
failure.prob.tee.lower<- apply(failure.prob.tee, 2, quantile, probs=0.025)
failure.prob.tee.upper<- apply(failure.prob.tee, 2, quantile, probs=0.975)
failure.prob.tee.lower10<- apply(failure.prob.tee, 2, quantile, probs=0.1)
failure.prob.tee.upper10<- apply(failure.prob.tee, 2, quantile, probs=0.9)

plot((t.seq-min(t.seq))/365, failure.prob.tee.mean, type="n", xlab="Years Since Installation", 
     ylab="Probability", xlim=c(min((t.seq-min(t.seq))/365),25),
     ylim=c(0.0, 0.2))

lines((t.seq-min(t.seq))/365, failure.prob.tee.mean, lty=1, col="blue", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.tee.lower, lty=2, col="deeppink", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.tee.upper, lty=2, col="deeppink", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.tee.lower10, lty=4, col="green", lwd=3.5)
lines((t.seq-min(t.seq))/365, failure.prob.tee.upper10, lty=4, col="green", lwd=3.5)

legend("bottomright", legend = c("cdf Estimate","95% Credible Intervals","80% Credible Intervals" ), col=c("blue","deeppink","green"), 
lty=c(1,2,4), merge=TRUE , border=c(0,0,1,0), lwd=c(2.5, 2.5, 2.5),cex=0.7)

mtext("TML Component: Tee", line=1)

# Page 23, Figure 13
#****************************
markyear <- seq(0,25,by = 5)
markday <- markyear*365+t.seq[1]

par(mfrow = c(2,2), mgp = c(4,1,0), mar = c(5,7,4,2))
plot(t.seq, failure.prob.elbow.mean, xlab = "Years After Installation",ylab = "Probability",
     type = "n", xaxt = "n",xlim = c(min(markday), max(markday)), ylim = c(0,0.6))
lines(t.seq, failure.prob.elbow.mean, col = "blue", lty = 2)
axis(1,at =  markday,labels =  markyear)

fprob.elbow.weibull <- read.table("failure.prob.elbow.weibull.csv", sep = ",")
fprob.elbow.mean.weibull <- colMeans(fprob.elbow.weibull, na.rm = FALSE, dims = 1)
lines(t.seq, fprob.elbow.mean.weibull, col = "red", lty = 3)
mtext("TML Component Elbow",cex = 0.7)
legend("topleft",cex = 0.6,bty = "n",lty = c(2,3), col = c("blue", "red"), c("Lognormal Corrosion Rate","Weibull Corrosion Rate"),box.lwd = 0.5)

plot(t.seq, failure.prob.pipe.mean, xlab = "Years After Installation",ylab = "Probability",
     type = "n", xaxt = "n",xlim = c(min(markday), max(markday)), ylim = c(0,0.6))
lines(t.seq, failure.prob.pipe.mean, col = "blue", lty = 2)
axis(1,at =  markday,labels =  markyear)

fprob.pipe.weibull <- read.table("failure.prob.pipe.weibull.csv", sep = ",")
fprob.pipe.mean.weibull <- colMeans(fprob.pipe.weibull, na.rm = FALSE, dims = 1)
lines(t.seq, fprob.pipe.mean.weibull, col = "red", lty = 3)
mtext("TML Component Pipe",cex = 0.7)
legend("topleft",cex = 0.6,bty = "n",lty = c(2,3), col = c("blue", "red"), c("Lognormal Corrosion Rate","Weibull Corrosion Rate"),box.lwd = 0.5)

plot(t.seq, failure.prob.tee.mean, xlab = "Years After Installation",ylab = "Probability",
     type = "n", xaxt = "n",xlim = c(min(markday), max(markday)), ylim = c(0,0.6))
lines(t.seq, failure.prob.tee.mean, col = "blue", lty = 2)
axis(1,at =  markday,labels =  markyear)

fprob.tee.weibull <- read.table("failure.prob.tee.weibull.csv", sep = ",")
fprob.tee.mean.weibull <- colMeans(fprob.tee.weibull, na.rm = FALSE, dims = 1)
lines(t.seq, fprob.tee.mean.weibull, col = "red", lty = 3)
mtext("TML Component Pipe",cex = 0.7)
legend("topleft",cex = 0.6,bty = "n",lty = c(2,3), col = c("blue", "red"), c("Lognormal Corrosion Rate","Weibull Corrosion Rate"),box.lwd = 0.5)

# Page 26, Figure 16
#****************************
# Remaining Life Time of the Current Circuit
tc <- 12418
t.seq.new <- seq(from=tc, to=tc+365*20,by=10)

remaining.life.time.elbow <- matrix(NA, nrow=dim(data)[1], ncol=length(t.seq.new))
remaining.life.time.pipe <- matrix(NA, nrow=dim(data)[1], ncol=length(t.seq.new))
remaining.life.time.tee <- matrix(NA, nrow=dim(data)[1], ncol=length(t.seq.new))

failure.prob.elbow <- matrix(NA, nrow=dim(data)[1], ncol=length(t.seq.new))
failure.prob.pipe <- matrix(NA, nrow=dim(data)[1], ncol=length(t.seq.new))
failure.prob.tee <- matrix(NA, nrow=dim(data)[1], ncol=length(t.seq.new))

failure.prob.tc.elbow <- NULL
failure.prob.tc.pipe <- NULL
failure.prob.tc.tee <- NULL

N <-1000
set.seed(198411)

for (i in 1: dim(data)[1]) {
  mcmc.run <- data[i,]
  y0.elbow <- rnorm(N, mean = mcmc.run[1], sd =  mcmc.run[4])
  y0.pipe <- rnorm(N, mean = mcmc.run[2], sd =  mcmc.run[4])
  y0.tee <- rnorm(N, mean =  mcmc.run[3],  sd =  mcmc.run[4])
  beta1 <- exp(rnorm(N, mean = mcmc.run[10], sd = mcmc.run[6]))
  TI <- exp(rnorm(N, mean = mcmc.run[7], sd = mcmc.run[8]))
  
  for (j in 1: length(t.seq.new)){
    t <- rep(t.seq.new[j],N)
    failure.prob.elbow[i,j] <- sum(y0.elbow-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.elbow )/N
    failure.prob.pipe[i,j] <- sum(y0.pipe-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.pipe )/N
    failure.prob.tee[i,j] <- sum(y0.tee-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.tee )/N
  }
  failure.prob.tc.elbow <-  sum(y0.elbow-beta1*(tc-TI)/365*ifelse(t>TI, 1, 0)< Df.elbow )/N
  failure.prob.tc.pipe <-  sum(y0.elbow-beta1*(tc-TI)/365*ifelse(t>TI, 1, 0)< Df.elbow )/N
  failure.prob.tc.tee <-  sum(y0.elbow-beta1*(tc-TI)/365*ifelse(t>TI, 1, 0)< Df.elbow )/N
  
  remaining.life.time.elbow[i,] <- (failure.prob.elbow[i,]-failure.prob.tc.elbow)/(1-failure.prob.tc.elbow)
  remaining.life.time.pipe[i,] <- (failure.prob.pipe[i,]-failure.prob.tc.pipe)/(1-failure.prob.tc.pipe)
  remaining.life.time.tee[i,] <- (failure.prob.tee[i,]-failure.prob.tc.tee)/(1-failure.prob.tc.tee)
}

t.lower <- min(t.seq.new)
t.upper <- max(t.seq.new)

write.table(remaining.life.time.elbow, file = "remaining.life.time.elbow.lognormal.csv", sep=",")
write.table(remaining.life.time.pipe, file = "remaining.life.time.pipe.lognormal.csv", sep=",")
write.table(remaining.life.time.tee, file = "remaining.life.time.tee.lognormal.csv", sep=",")

remaining.life.time.elbow <- read.table(file="remaining.life.time.elbow.lognormal.csv",sep=",")
remaining.life.time.pipe <- read.table(file="remaining.life.time.pipe.lognormal.csv",sep=",")
remaining.life.time.tee <- read.table(file="remaining.life.time.tee.lognormal.csv",sep=",")

remaining.life.time.mean.elbow <- colMeans(remaining.life.time.elbow , na.rm = FALSE, dims = 1)
remaining.life.time.lower.elbow <- apply(remaining.life.time.elbow , 2, quantile, probs=0.025)
remaining.life.time.upper.elbow <- apply(remaining.life.time.elbow , 2, quantile, probs=0.975)
remaining.life.time.lower10.elbow <- apply(remaining.life.time.elbow , 2, quantile, probs=0.1)
remaining.life.time.upper10.elbow <- apply(remaining.life.time.elbow , 2, quantile, probs=0.9)

plot((t.seq.new-min(t.seq.new))/365, remaining.life.time.mean.elbow, type="n", 
     xlab="Years Since the Last Inspection", 
     ylab="Probability", xlim=c(min((t.seq-min(t.seq))/365),24),
     ylim=c(0.0, 0.2))

lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.mean.elbow, lty=1, col="blue", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.lower.elbow, lty=2, col="deeppink", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.upper.elbow, lty=2, col="deeppink", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.lower10.elbow, lty=4, col="green", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.upper10.elbow, lty=4, col="green", lwd=3.5)

legend("bottomright", legend = c("cdf Estimate","95% Credible Intervals","80% Credible Intervals" ), col=c("blue","deeppink","green"), 
       lty=c(1,2,3), merge=TRUE , border=c(0,0,1,0), lwd=c(2.5, 2.5, 2.5),cex=0.7)

mtext("TML Component: Elbow", line=1)

remaining.life.time.mean.pipe <- colMeans(remaining.life.time.pipe, na.rm = FALSE, dims = 1)
remaining.life.time.lower.pipe <- apply(remaining.life.time.pipe, 2, quantile, probs=0.025)
remaining.life.time.upper.pipe  <- apply(remaining.life.time.pipe, 2, quantile, probs=0.975)
remaining.life.time.lower10.pipe  <- apply(remaining.life.time.pipe, 2, quantile, probs=0.1)
remaining.life.time.upper10.pipe  <- apply(remaining.life.time.pipe, 2, quantile, probs=0.9)

plot((t.seq.new-min(t.seq.new))/365, remaining.life.time.mean.pipe, type="n", 
     xlab="Years Since the Last Inspection", 
     ylab="Probability", xlim=c(min((t.seq-min(t.seq))/365),24),
     ylim=c(0.0, 0.2))

lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.mean.pipe, lty=1, col="blue", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.lower.pipe, lty=2, col="deeppink", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.upper.pipe, lty=2, col="deeppink", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.lower10.pipe, lty=4, col="green", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.upper10.pipe, lty=4, col="green", lwd=3.5)

legend("bottomright", legend = c("cdf Estimate","95% Credible Intervals","80% Credible Intervals" ), col=c("blue","deeppink","green"), 
lty=c(1,2,3), merge=TRUE , border=c(0,0,1,0), lwd=c(2.5, 2.5, 2.5),cex=0.7)

mtext("TML Component: Pipe", line=1)

remaining.life.time.mean.tee <- colMeans(remaining.life.time.tee, na.rm = FALSE, dims = 1)
remaining.life.time.lower.tee <- apply(remaining.life.time.tee, 2, quantile, probs=0.025)
remaining.life.time.upper.tee  <- apply(remaining.life.time.tee, 2, quantile, probs=0.975)
remaining.life.time.lower10.tee  <- apply(remaining.life.time.tee, 2, quantile, probs=0.1)
remaining.life.time.upper10.tee  <- apply(remaining.life.time.tee, 2, quantile, probs=0.9)

plot((t.seq.new-min(t.seq.new))/365, remaining.life.time.mean.tee, type="n", xlab="Years Since the Last Inspection",
     ylab="Probability", xlim=c(min((t.seq-min(t.seq))/365),24),
     ylim=c(0.0, 0.2))

lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.mean.tee, lty=1, col="blue", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.lower.tee, lty=2, col="deeppink", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.upper.tee, lty=2, col="deeppink", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.lower10.tee, lty=4, col="green", lwd=3.5)
lines((t.seq.new-min(t.seq.new))/365, remaining.life.time.upper10.tee, lty=4, col="green", lwd=3.5)

legend("bottomright", legend = c("cdf Estimate","95% Credible Intervals","80% Credible Intervals" ), col=c("blue","deeppink","green"), 
lty=c(1,2,3), merge=TRUE , border=c(0,0,1,0), lwd=c(2.5, 2.5, 2.5),cex=0.7)

mtext("TML Component: Tee", line=1)

# Page 25, Figure 15
#****************************
markyear <- seq(0,25,by = 5)
markday <- markyear*365+t.seq.new[1]

par(mfrow = c(2,2), mgp = c(4,1,0), mar = c(5,7,4,2))
plot(t.seq.new, remaining.life.time.mean.elbow, xlab = "Years After the Last Inspection",ylab = "Probability",
     type = "n", xaxt = "n",xlim = c(min(markday), max(markday)), ylim = c(0,0.6))
lines(t.seq.new, remaining.life.time.mean.elbow, col = "blue", lty = 2)
axis(1,at =  markday[-6],labels =  markyear[-6])

rlt.elbow.weibull <- read.table("remaining.life.time.elbow.weibull.csv", sep = ",")
rlt.elbow.mean.weibull <- colMeans(rlt.elbow.weibull, na.rm = FALSE, dims = 1)
lines(t.seq.new, rlt.elbow.mean.weibull, col = "red", lty = 3)
mtext("TML Component Elbow",cex = 0.7)
legend("bottomright",cex = 0.6,bty = "n",lty = c(2,3), col = c("blue", "red"), c("Lognormal Corrosion Rate","Weibull Corrosion Rate"),box.lwd = 0.5)

plot(t.seq.new, remaining.life.time.mean.pipe, xlab = "Years After the Last Inspection",ylab = "Probability",
     type = "n", xaxt = "n",xlim = c(min(markday), max(markday)), ylim = c(0,0.6))
lines(t.seq.new, remaining.life.time.mean.pipe, col = "blue", lty = 2)
axis(1,at =  markday[-6],labels =  markyear[-6])

rlt.pipe.weibull <- read.table("remaining.life.time.pipe.weibull.csv", sep = ",")
rlt.pipe.mean.weibull <- colMeans(rlt.pipe.weibull, na.rm = FALSE, dims = 1)
lines(t.seq.new, rlt.pipe.mean.weibull, col = "red", lty = 3)
mtext("TML Component Pipe",cex = 0.7)
legend("bottomright",cex = 0.6,bty = "n",lty = c(2,3), col = c("blue", "red"), c("Lognormal Corrosion Rate","Weibull Corrosion Rate"),box.lwd = 0.5)


plot(t.seq.new, remaining.life.time.mean.tee, xlab = "Years After the Last Inspection",ylab = "Probability",
     type = "n", xaxt = "n",xlim = c(min(markday), max(markday)), ylim = c(0,0.6))
lines(t.seq.new, remaining.life.time.mean.tee, col = "blue", lty = 2)
axis(1,at =  markday[-6],labels =  markyear[-6])

rlt.tee.weibull <- read.table("remaining.life.time.tee.weibull.csv", sep = ",")
rlt.tee.mean.weibull <- colMeans(rlt.tee.weibull, na.rm = FALSE, dims = 1)
lines(t.seq.new, rlt.tee.mean.weibull, col = "red", lty = 3)
mtext("TML Component Tee",cex = 0.7)
legend("bottomright",cex = 0.6,bty = "n",lty = c(2,3), col = c("blue", "red"), c("Lognormal Corrosion Rate","Weibull Corrosion Rate"),box.lwd = 0.5)

# Page 27, Figure 17
#****************************
# Elbow
remaining.life.time <- remaining.life.time.elbow
p <- c(0.1, 0.2, 0.3, 0.4)
M <- 100
p.star <- 1-(1-p)^(1/M)

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[1])
  min.diff0 <- min(diff0)
  index <- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar1.elbow <- t.out

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[2])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar2.elbow <- t.out

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[3])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar3.elbow <- t.out

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[4])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <-  t.seq.new[index]
}

t.out.pstar4.elbow <- t.out

# Pipe
remaining.life.time <- remaining.life.time.pipe

p <- c(0.1, 0.2, 0.3, 0.4)
M <- 100

p.star <- 1-(1-p)^(1/M)


t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[1])
  min.diff0 <- min(diff0)
  index <- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar1.pipe <- t.out

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[2])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar2.pipe <- t.out

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[3])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar3.pipe <- t.out

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[4])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <-  t.seq.new[index]
}

t.out.pstar4.pipe <- t.out

# Tee
remaining.life.time <- remaining.life.time.tee

p <- c(0.1, 0.2, 0.3, 0.4)
M <- 100

p.star <- 1-(1-p)^(1/M)

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[1])
  min.diff0 <- min(diff0)
  index <- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar1.tee <- t.out

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[2])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar2.tee <- t.out

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[3])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar3.tee <- t.out

t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[4])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <-  t.seq.new[index]
}

t.out.pstar4.tee <- t.out

t.out.pstar1.elbow <- (t.out.pstar1.elbow-min(t.seq))/365
t.out.pstar2.elbow <- (t.out.pstar2.elbow-min(t.seq))/365
t.out.pstar3.elbow <- (t.out.pstar3.elbow-min(t.seq))/365
t.out.pstar4.elbow <- (t.out.pstar4.elbow-min(t.seq))/365

t.out.pstar1.pipe <- (t.out.pstar1.pipe-min(t.seq))/365
t.out.pstar2.pipe <- (t.out.pstar2.pipe-min(t.seq))/365
t.out.pstar3.pipe <- (t.out.pstar3.pipe-min(t.seq))/365
t.out.pstar4.pipe <- (t.out.pstar4.pipe-min(t.seq))/365

t.out.pstar1.tee <- (t.out.pstar1.tee-min(t.seq))/365
t.out.pstar2.tee <- (t.out.pstar2.tee-min(t.seq))/365
t.out.pstar3.tee <- (t.out.pstar3.tee-min(t.seq))/365
t.out.pstar4.tee <- (t.out.pstar4.tee-min(t.seq))/365

tout1.elbow <- density(t.out.pstar1.elbow,  bw = 1.3)
tout2.elbow <- density(t.out.pstar2.elbow,  bw = 1.3)
tout3.elbow <- density(t.out.pstar3.elbow,  bw = 1.3)
tout4.elbow <- density(t.out.pstar4.elbow,  bw = 1.3)

tout1.pipe <- density(t.out.pstar1.pipe,  bw = 1.3)
tout2.pipe <- density(t.out.pstar2.pipe,  bw = 1.3)
tout3.pipe <- density(t.out.pstar3.pipe,  bw = 1.3)
tout4.pipe <- density(t.out.pstar4.pipe,  bw = 1.3)

tout1.tee <- density(t.out.pstar1.tee,  bw = 1.3)
tout2.tee <- density(t.out.pstar2.tee,  bw = 1.3)
tout3.tee <- density(t.out.pstar3.tee,  bw = 1.3)
tout4.tee <- density(t.out.pstar4.tee,  bw = 1.3)


par(las=1)
plot(tout1.elbow$x, -tout1.elbow$y*4, ylim=c(-0.5, 4),xlim=c(0, 18),type='n', 
     ylab="Quantiles of the Posterior Minimum Remaining Lifetime Distribution", 
     xlab="Years Since January 2004", yaxt="n")

axis(2,c(0,1,2,3), c(0.1, 0.2, 0.3, 0.4),las=1,cex.axis=0.75)

lines(tout1.elbow$x,tout1.elbow$y*4, lty=4, col="red", lwd=2.5)
lines(tout2.elbow$x,tout2.elbow$y*4+1, lty=4, col="red", lwd=2.5)
lines(tout3.elbow$x,tout3.elbow$y*4+2.0, lty=4, col="red", lwd=2.5)
lines(tout4.elbow$x,tout4.elbow$y*4+3, lty=4, col="red", lwd=2.5)

lines(tout1.pipe$x,tout1.pipe$y*4, lty=2, col="blue", lwd=2.5)
lines(tout2.pipe$x,tout2.pipe$y*4+1, lty=2, col="blue", lwd=2.5)
lines(tout3.pipe$x,tout3.pipe$y*4+2.0, lty=2, col="blue", lwd=2.5)
lines(tout4.pipe$x,tout4.pipe$y*4+3, lty=2, col="blue", lwd=2.5)

lines(tout1.tee$x,tout1.tee$y*4, lty=1, col="green", lwd=2.5)
lines(tout2.tee$x,tout2.tee$y*4+1, lty=1, col="green", lwd=2.5)
lines(tout3.tee$x,tout3.tee$y*4+2.0, lty=1, col="green", lwd=2.5)
lines(tout4.tee$x,tout4.tee$y*4+3, lty=1, col="green", lwd=2.5)

legend("topright", legend = c("Elbow","Pipe","Tee"), 
       col=c("red","blue","green"), 
       lty=c(4,2,1), merge=TRUE ,lwd=c(2.5, 2.5, 2.5),cex=0.6)

abline(h=0, col="grey")
abline(h=1, col="grey")
abline(h=2, col="grey")
abline(h=3, col="grey")

# Page 20, Figure 10
#****************************
index1 <- sample(32000, 150)
beta.sub <- as.data.frame(beta[index1,])
colnames(beta.sub) <- c(1:33)
boxplot(beta.sub, ylab = expression(paste("Corrosion Rate: ","\u03b2"[1[i]])), xlab = "TML Number",ylim = c(0,0.25), col = "yellow")

# Page 21, Figure 11
#****************************
ydate <- as.date(c("1Apr2000","1Nov2000","1Jul2001","1Jan2004","1Jan2006","1Jan2008","1Jan2010","1Jan2012","1Jan2014"))
yday <- julian(as.Date(ydate,"%d%b%Y"))
init.sub <- as.data.frame(init[index1,])
colnames(init.sub) <- c(1:33)
boxplot(init.sub, yaxt = "n",ylab = expression(paste("Initiation Time: ","T"[I[i]])), xlab = "TML Number", ylim = c(10500,16500), col = "yellow")
axis(side = 2, yday, ydate,cex.axis = 0.6)

# Page 13, Figure 7
#****************************
Date = NULL;Thickness = NULL;TML = NULL

for(i in 1:33){
  for(j in 1:4){
    Date <- c(Date, JulianDate[j])
    TML <- c(TML, i)
    Thickness <- c(Thickness, thickness.fitted.measure[i,j])
  }
}

thick.data.frame <- data.frame(Date = Date, Thickness = Thickness, TML = as.factor(TML))

mypanel1 <- function(x,y){
  panel.lines(x,y)
  panel.points(x,y,col = as.factor(x),pch =16)
}

xyplot(data = thick.data.frame, Thickness~ Date|TML,ylab = "Wall Thickness",panel = mypanel1,
       scales = list(x = list(at = Date[1:4], labels = date.tmpt[-5], rot = 90, cex = 0.6)),
       key=list(space="right",lines=list(col=c("black","red","green","blue"),type = "p", pch = 16),
                text=list(c("Apr2000","Nov2000","Jul2001","Jan2004"))))

# Page 18, Figure 8
#*************************
pmean <- get_posterior_mean(fit3)[,5]
#initial.time <- pmean[76:108]
initial.time <- apply(init,2,median)

# Assume that all the corrosions(if started) are between the 3rd and 4th inspections.
# If not, changes should be made for the "TMLs with Initiation"

wi <- which((initial.time - JulianDate[4])<0)
woi <- c(1:33)[-wi]
# Without Initiation
thick.no.init <- thickness.fitted.measure[woi,]

Date = NULL;Thickness = NULL;TML = NULL;InitTime = NULL; 
CorrRate = NULL; InitThick = NULL;

#initial.time.woi <- initial.time[woi]
initial.time.woi <- initial.time[woi]
corrosion.rate.woi <- beta.mean[woi]
#initial.thick.woi <- thickness.init.mean[woi]
initial.thick.woi <- apply(thickness.init,2,median)[woi]

for(i in 1:length(woi)){
  for(j in 1:4){
    Date <- c(Date, JulianDate[j])
    TML <- c(TML, woi[i])
    Thickness <- c(Thickness, thick.no.init[i,j])
  }
}

thick.no.init.data.frame <- data.frame(Date = Date,Thickness = Thickness,TML = as.factor(TML))


fit.date <- c("1Jan2004","1Jan2006","1Jan2008","1Jan2010","1Jan2012","1Jan2014")
fit.date.date <- date::as.date(fit.date)
fit.day <- julian(as.Date(fit.date.date))
xyplot(data = thick.no.init.data.frame, Thickness~Date|TML, ylab = "Fitted Wall Thickness",
       xlab = "Date", main = "TMLs without Initiation",
       xlim = c(11000, 16300), ylim = c(0, 0.6),
       scales = list(x = list(at = fit.day, labels = fit.date, rot = 90, cex = 0.6)),
       panel = function(x,y){
         panel.lines(x,y)
         panel.points(x,y, col = c("black","red","green","blue"),pch = 16, cex = 0.6)
         panel.points(initial.time.woi[panel.number()],
                      initial.thick.woi[panel.number()],pch = 4, col = "deeppink", cex = 0.6)
         panel.lines(c(x[4],initial.time.woi[panel.number()]),c(y[4],initial.thick.woi[panel.number()]))
         panel.abline(h = 0.2*initial.thick.woi[panel.number()], lty = 3, col = "blue")
         lastfit <- initial.thick.woi[panel.number()] - (16071 - initial.time.woi[panel.number()])/365 * corrosion.rate.woi[panel.number()]
         panel.lines(c(initial.time.woi[panel.number()],16071),c(initial.thick.woi[panel.number()], lastfit),
                     lty = 3,col = "grey")
       },
       key=list(space="right",lines=list(col=c("black","red","green","blue","deeppink"),type = "p", pch = c(16,16,16,16,4),cex = 0.75),
                text=list(c("Apr2000","Nov2000","Jul2001","Jan2004","Estimated Initiation"),cex = 0.75)
                )
       )

# With Initiation
thick.with.init <- thickness.fitted.measure[wi,]
Date = NULL;Thickness = NULL;TML = NULL;InitTime = NULL; 
CorrRate = NULL; InitThick = NULL;

initial.time.wi <- initial.time[wi]
corrosion.rate.wi <- beta.mean[wi]
initial.thick.wi <- apply(thickness.init,2,median)[wi]

for(i in 1:length(wi)){
  for(j in 1:4){
    Date <- c(Date, JulianDate[j])
    TML <- c(TML, wi[i])
    Thickness <- c(Thickness, thick.with.init[i,j])
  }
}

thick.with.init.data.frame <- data.frame(Date = Date,Thickness = Thickness,TML = as.factor(TML))

fit.date <- c("1Jan2004","1Jan2006","1Jan2008","1Jan2010","1Jan2012","1Jan2014")
fit.date.date <- date::as.date(fit.date)
fit.day <- julian(as.Date(fit.date.date))

xyplot(data = thick.with.init.data.frame, Thickness~Date|TML, ylab = "Fitted Wall Thickness",
       xlab = "Date", main = "TMLs with Initiation",
       xlim = c(11000, 16300), ylim = c(0, 0.6),
       scales = list(x = list(at = fit.day, labels = fit.date, rot = 90, cex = 0.6)),
       panel = function(x,y){
         itime <- initial.time.wi[panel.number()]
         ithick <- initial.thick.wi[panel.number()]
         panel.lines(c(x[-4],itime),c(y[-4],ithick))
         lastfit <- ithick - (16071 - itime)/365 * corrosion.rate.wi[panel.number()]
         panel.lines(c(itime,16071),c(ithick, lastfit),
                     lty = 3,col = "grey")
         panel.points(x,y,col = c("black","red","green","blue"), pch = 16, cex = 0.6)
         panel.points(itime,ithick,pch = 4, col = "deeppink", cex = 0.6)
         panel.abline(h = 0.2*initial.thick.wi[panel.number()], lty = 3, col = "blue")
         },
       key=list(space="right",lines=list(col=c("black","red","green","blue","deeppink"),type = "p", pch = c(16,16,16,16,4),cex = 0.75),
                text=list(c("Apr2000","Nov2000","Jul2001","Jan2004","Estimated Initiation"),cex = 0.75)
       )
)

# Page 22, Figure 12
#****************************
date.axis <- c("1Apr2000","1Nov2000","1Jun2001","1Jan2004","1Jan2006",
               "1Jan2008")
day.axis <- date::as.date(date.axis)
day.axis <- julian(as.Date(day.axis))

density.model <- apply(init[sample(dim(init)[1],1000),],2,density,width = 600)

par(las = 1, mgp = c(4,1,0), mar = c(6,7,4,2))
plot(density.model[[1]]$x,density.model[[1]]$y, xlim = c(10500,14000), ylim = c(0,0.002),type = "n",
     xlab = "Date", ylab = "Density",xaxt = "n")
axis(1, at=day.axis, labels = c("1Apr00","1Nov00","1Jun01","1Jan04","1Jan06",
                                "1Jan08"), las = 2,cex = 0.5)
legend(x = 12418, y = 0.002, c("TMLs with Initiation", "TMLs without Initiation"), lty = c(1,3), col = c("blue","red"), cex = 0.6)
for(i in 1:dim(init)[2]){
  type = i %in% wi
  if(type){
    ltype = 1
    color = "blue"
  }
  else {
    ltype = 3
    color = "red"
  }
  lines(density.model[[i]]$x,density.model[[i]]$y, lty = ltype,
        col = color)
}
