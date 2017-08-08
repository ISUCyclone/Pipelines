# Page 9, Figure 4
#******************************
pmean <- get_posterior_mean(fit2)
slope <- pmean[7:94]
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