library(bootstrap) # load up the bootstrap package

set.seed(23434)

x1 <- rnorm(40,0,1)
x2 <- rnorm(40,0,1)
y <- 10 + x1*2 - 5*x2 + rnorm(40,0,1)

reg.1 <- lm(y ~ x1 + x2)
summary(reg.1)

x1[10] <- 10 # replace the 10th x1 observation with 10
my.data <- as.data.frame(cbind(y,x1,x2)) 

plot(y,x1)

reg.jack <- lm(y ~ x1 + x2)
summary(reg.jack) 


theta <- function(x, dat, coefficient){
  coef(lm(reg.jack , data = dat[x,]))[coefficient]
}


res <- jackknife(1:length(x1), theta, dat = my.data, coefficient = "x1")

summary(res)

res

plot(res$jack.values,y)