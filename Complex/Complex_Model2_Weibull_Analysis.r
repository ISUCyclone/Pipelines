post4 <- as.array(fit4)
post4 <- rbind(post4[,1,], post4[,2,],post4[,3,],post4[,4,])

beta <- post4[,43:75]
init <- post4[,76:108]
thickness.int <- post4[,10:42]

beta.mean <- apply(beta,2, mean)
init.mean <- apply(init,2, mean)
thickness.int.mean <- apply(thickness.int,2, mean)

index<- sample(dim(post4)[1], 1000)
beta.sub <- as.data.frame(beta[index,])
boxplot(beta.sub)

data <- post4[,c(5,6,7,1,2,3,110,109,8,4)]

data<- data[index,]

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

Df.elbow <- mean(Y[1:12,1])*.3; Df.elbow 
Df.pipe <- mean(Y[13:28,1])*.3; Df.pipe
Df.tee <- mean(Y[29:33,1])*.3; Df.tee 

library(date)
date.tmpt <- as.date(c("1apr2000","1nov2000","1jul2001","1jan2004", "1jan2025"),order="dmy")
JulianDate <- julian(as.Date(date.tmpt, "%d%b%Y"))

t.lower <- min(JulianDate)
t.upper <- max(JulianDate)

t.seq <- signif(seq(t.lower, t.upper, length=100),5)

failure.prob.elbow <- matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq))
failure.prob.pipe <-  matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq))
failure.prob.tee <- matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq))

N <-1000
for (i in 1: dim(data)[[1]]) {
  mcmc.run <- data[i,]
  ##############
  # Step 1
  ##############
  # For each type of the TML (elbow, pipe and tee), generate N=100,000
  # y0.elbows, y0.pipes, y0.tees,  beta_1s and T_Is.
  set.seed(198411)
  y0.elbow <- rnorm(N, mean = mcmc.run[4], sd =  mcmc.run[10])
  y0.pipe <- rnorm(N, mean = mcmc.run[5], sd =  mcmc.run[10])
  y0.tee <- rnorm(N, mean =  mcmc.run[6],  sd =  mcmc.run[10])
  beta1 <- rweibull(N, shape = mcmc.run[7], scale = mcmc.run[8])
  TI <- exp(rnorm(N, mean = mcmc.run[3], sd = mcmc.run[9]))
  
  for (j in 1: length(t.seq)){
    t <- rep(t.seq[j],N)
    failure.prob.elbow[i,j] <- sum(y0.elbow-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.elbow )/N
    failure.prob.pipe[i,j] <- sum(y0.pipe-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.pipe )/N
    failure.prob.tee[i,j] <- sum(y0.tee-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.tee )/N
  }
}

write.table(failure.prob.elbow, file = "failure.prob.elbow.weibull.csv", sep=",")
write.table(failure.prob.pipe, file = "failure.prob.pipe.weibull.csv", sep=",")
write.table(failure.prob.tee, file = "failure.prob.tee.weibull.csv", sep=",")

failure.prob.elbow.mean <- colMeans(failure.prob.elbow, na.rm = FALSE, dims = 1)
failure.prob.pipe.mean <- colMeans(failure.prob.pipe, na.rm = FALSE, dims = 1)
failure.prob.tee.mean <- colMeans(failure.prob.tee, na.rm = FALSE, dims = 1)

failure.prob.elbow.lower<- apply(failure.prob.elbow, 2, quantile, probs=0.025)
failure.prob.elbow.upper<- apply(failure.prob.elbow, 2, quantile, probs=0.975)
failure.prob.elbow.lower10<- apply(failure.prob.elbow, 2, quantile, probs=0.1)
failure.prob.elbow.upper10<- apply(failure.prob.elbow, 2, quantile, probs=0.9)

par(las=1)
plot((t.seq-min(t.seq))/365, failure.prob.elbow.mean, type="n", xlab="Years Since Installation", 
     ylab="Probability", xlim=c(min((t.seq-min(t.seq))/365),max((t.seq-min(t.seq))/365)),
     ylim=c(0.0, 0.2))

lines((t.seq-min(t.seq))/365, failure.prob.elbow.mean, lty=1, col="blue", lwd=2.5)
lines((t.seq-min(t.seq))/365, failure.prob.elbow.lower, lty=2, col="deeppink", lwd=2.5)
lines((t.seq-min(t.seq))/365, failure.prob.elbow.upper, lty=2, col="deeppink", lwd=2.5)
lines((t.seq-min(t.seq))/365, failure.prob.elbow.lower10, lty=3, col="green", lwd=2.5)
lines((t.seq-min(t.seq))/365, failure.prob.elbow.upper10, lty=3, col="green", lwd=2.5)

legend("bottomright", legend = c("cdf Estimate","95% Credible Intervals","80% Credible Intervals" ), col=c("blue","deeppink","green"), 
       lty=c(1,2,3), merge=TRUE , border=c(0,0,1,0), lwd=c(2.5, 2.5, 2.5),cex=0.7)

mtext("TML Feature: Elbow", line=1)

# Remaining Lifetime
tc <- 12418
t.seq.new <- seq(from=tc, to=tc+365*20,by=10)

remaining.life.time.elbow.w <- matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq.new))
remaining.life.time.pipe.w <- matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq.new))
remaining.life.time.tee.w <- matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq.new))

failure.prob.elbow.w <- matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq.new))
failure.prob.pipe.w <- matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq.new))
failure.prob.tee.w <- matrix(NA, nrow=dim(data)[[1]], ncol=length(t.seq.new))

failure.prob.tc.elbow.w <- NULL
failure.prob.tc.pipe.w <- NULL
failure.prob.tc.tee.w <- NULL

for (i in 1: dim(data)[1]) {
  mcmc.run <- data[i,]
  set.seed(198411)
  y0.elbow <- rnorm(N, mean = mcmc.run[4], sd =  mcmc.run[10])
  y0.pipe <- rnorm(N, mean = mcmc.run[5], sd =  mcmc.run[10])
  y0.tee <- rnorm(N, mean =  mcmc.run[6],  sd =  mcmc.run[10])
  beta1 <- rweibull(N, shape = mcmc.run[7], scale = mcmc.run[8])
  TI <- exp(rnorm(N, mean = mcmc.run[3], sd = mcmc.run[9]))
  
  failure.prob.tc.elbow.w <-  sum(y0.elbow-beta1*(tc-TI)/365*ifelse(t>TI, 1, 0)< Df.elbow )/N
  failure.prob.tc.pipe.w <-  sum(y0.elbow-beta1*(tc-TI)/365*ifelse(t>TI, 1, 0)< Df.elbow )/N
  failure.prob.tc.tee.w<-  sum(y0.elbow-beta1*(tc-TI)/365*ifelse(t>TI, 1, 0)< Df.elbow )/N
  for (j in 1: length(t.seq.new)){
    t <- rep(t.seq.new[j],N)
    failure.prob.elbow.w[i,j] <- sum(y0.elbow-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.elbow )/N
    failure.prob.pipe.w[i,j] <- sum(y0.pipe-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.pipe )/N
    failure.prob.tee.w[i,j] <- sum(y0.tee-beta1*(t-TI)/365*ifelse(t>TI, 1, 0)< Df.tee )/N
  }
  remaining.life.time.elbow.w[i,] <- (failure.prob.elbow.w[i,]-failure.prob.tc.elbow.w)/(1-failure.prob.tc.elbow.w)
  remaining.life.time.pipe.w[i,] <- (failure.prob.pipe.w[i,]-failure.prob.tc.pipe.w)/(1-failure.prob.tc.pipe.w)
  remaining.life.time.tee.w[i,] <- (failure.prob.tee.w[i,]-failure.prob.tc.tee.w)/(1-failure.prob.tc.tee.w)
}

t.lower <- min(t.seq.new)
t.upper <- max(t.seq.new)

write.table(remaining.life.time.elbow.w, file = "remaining.life.time.elbow.weibull.csv", sep=",")
write.table(remaining.life.time.pipe.w, file = "remaining.life.time.pipe.weibull.csv", sep=",")
write.table(remaining.life.time.tee.w, file = "remaining.life.time.tee.weibull.csv", sep=",")


remaining.life.time.mean.elbow.w <- colMeans(remaining.life.time.elbow.w , na.rm = FALSE, dims = 1)
remaining.life.time.lower.elbow.w <- apply(remaining.life.time.elbow.w , 2, quantile, probs=0.025)
remaining.life.time.upper.elbow.w <- apply(remaining.life.time.elbow.w , 2, quantile, probs=0.975)
remaining.life.time.lower10.elbow.w <- apply(remaining.life.time.elbow.w , 2, quantile, probs=0.1)
remaining.life.time.upper10.elbow.w <- apply(remaining.life.time.elbow.w , 2, quantile, probs=0.9)

