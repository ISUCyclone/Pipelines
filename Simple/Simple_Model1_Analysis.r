# Page 7, Figure 3
#******************************
pmean <- get_posterior_mean(fit1)[,5]
slope <- pmean[4:91]
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

par(las = 1,mgp = c(4,1,0),mar = c(4,7,4,2))
plot(days,thick[1,],type = "n", xlim = c(x.min,x.max),
     ylim = c(0, 0.26), xlab = "Date", ylab = "Fitted Wall Thickness",
     xaxt = "n")

display.days <- days[-c(4,5,7)]
display.dates <- dates[-c(4,5,7)]
axis(1,display.days,display.dates,las=2,cex.axis = 0.75)

for(i in 1:88) {
  lines(days, thick[i,], type = "l", lwd = 1, lty = 1)
}

abline(h = 0.05, lwd = 1, col = "blue")
abline(v = 12072, lwd = 1.2, col = "blue", lty = 2)
abline(v = 9172, lwd = 1.2, col = "blue", lty = 2)

arrows(8500, 0.08, 9000, 0.08,col = "blue",length = 0.08)
text(7800,0.08, "First\nInspection", cex = 0.9, col = "blue")

arrows(11400, 0.08, 11900, 0.08, col = "blue", length = 0.08)
text(10700, 0.08, "Last\nInspection",cex = 0.9, col = "blue")

text(13800,00.03,"Critical Level",cex = 0.9, col = "blue")
