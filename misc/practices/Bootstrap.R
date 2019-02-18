set.seed(23434)

x1 <- rnorm(40,0,1)
x2 <- rnorm(40,0,1)
y <- 10 + x1*2 - 5*x2 + rnorm(40,0,1)

reg.1 <- lm(y ~ x1 + x2)
summary(reg.1)

library(boot)

my.data <- as.data.frame(cbind(y,x1,x2)) 

bootstrap <- function(formula, data, regressors) {
  dat <- data[regressors,]	# grabs the sample
  reg <- lm(formula, data = dat) # runs the regression
  return(coef(reg)) # we need these coefficients
}

bs.res <- boot(R = 1000, formula = y~x1 + x2 , data = my.data, statistic = bootstrap)
summary(bs.res)  # note everything that is saved

bs.res$t0 # index


#Basic bootstrap

bs.basic.x1 <- boot.ci(bs.res, type = "basic", index = 2) # 95% for variable x1
bs.basic.x1

bs.basic.x2 <- boot.ci(bs.res, type = "basic", index = 3) # 95% for variable x2
bs.basic.x2


#Adjusted bootstrap percentile (BCa)

bs.bca.x1 <- boot.ci(bs.res, type = "bca", index = 2)
bs.bca.x1


# Normal

bs.norm.x1 <- boot.ci(bs.res, type = "norm", index = 2)
bs.norm.x1

bs.norm.x2 <- boot.ci(bs.res, type = "norm", index = 3)
bs.norm.x2


#Percentile interval

bs.perc.x1 <- boot.ci(bs.res, type = "perc", index = 2)
bs.perc.x1

bs.perc.x2 <- boot.ci(bs.res, type = "perc", index = 3)
bs.perc.x2


#For X1:

plot(NULL, type = "l", xlim = c(1,3),ylim = c(0,5), ylab = NA, axes = FALSE, xlab = NA)
lines(c(bs.basic.x1$basic[4], bs.basic.x1$basic[5]), c(1,1))  # add the 95% confidence intervals
text(2, 1.2, "Basic", xpd = T, cex = .8)     # add the names
lines(c(bs.norm.x1$norm[2], bs.norm.x1$norm[3]), c(2,2))
text(2, 2.2, "Normal", xpd = T, cex = .8)
lines(c(bs.bca.x1$bca[4], bs.bca.x1$bca[5]), c(3,3))
text(2, 3.2, "BCa", xpd = T, cex = .8)
lines(c(bs.perc.x1$perc[4], bs.perc.x1$perc[5]), c(4,4))
text(2, 4.2, "Percentile", xpd = T, cex = .8)
abline(v = 2.0, lty = 3, col = "black") #  add vertical line at x1
axis(side = 1)  # add x axis
mtext(side = 2, "Bootstrap Method", line = -3)    # label side
mtext(side = 3, "Performance of Various Bootstrapping Methods for X1", line = 1)   # add a title

#For X2:

plot(NULL, type = "l", xlim = c(-6,-4),ylim = c(0,5), ylab = NA, axes = FALSE, xlab = NA)
lines(c(bs.basic.x2$basic[4], bs.basic.x2$basic[5]), c(1,1))  # add the 95% confidence intervals
text(-5, 1.2, "Basic", xpd = T, cex = .8)     # add the names
lines(c(bs.norm.x2$norm[2], bs.norm.x2$norm[3]), c(2,2))
text(-5, 2.2, "Normal", xpd = T, cex = .8)
lines(c(bs.bca.x2$bca[4], bs.bca.x2$bca[5]), c(3,3))
text(-5, 3.2, "BCa", xpd = T, cex = .8)
lines(c(bs.perc.x2$perc[4], bs.perc.x2$perc[5]), c(4,4))
text(-5, 4.2, "Percentile", xpd = T, cex = .8)
abline(v = -5.0, lty = 3, col = "black") #  add vertical line at x2
axis(side = 1)  # add x axis
mtext(side = 2, "Bootstrap Method", line = -3)    # label side
mtext(side = 3, "Performance of Various Bootstrapping Methods for X2", line = 1)   # add a title
