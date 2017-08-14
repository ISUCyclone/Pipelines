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
Df.elbow <- mean(Y[1:12,1])*.3; Df.elbow 
Df.pipe <- mean(Y[13:28,1])*.3; Df.pipe
Df.tee <- mean(Y[29:33,1])*.3; Df.tee

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

# Page 26, Figure 16
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

# Page 27, Figure 17
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

#################


t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[2])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar2.elbow <- t.out


#################
t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[3])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar3.elbow <- t.out

#################
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

#################


t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[2])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar2.pipe <- t.out


#################
t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[3])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar3.pipe <- t.out

#################
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

#################


t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[2])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar2.tee <- t.out


#################
t.out  <- NULL
for (i in 1: dim(remaining.life.time)[1]) {
  diff0<- abs(remaining.life.time[i,]-p.star[3])
  min.diff0 <- min(diff0)
  index<- median(which(diff0==min.diff0))
  t.out[i] <- t.seq.new[index]
}

t.out.pstar3.tee <- t.out

#################
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
