x <-c(8.26, 6.33, 10.4, 5.27, 5.35, 5.61, 6.12, 6.19, 5.2,
      7.01, 8.74, 7.78, 7.02, 6, 6.5, 5.8, 5.12, 7.41, 6.52, 6.21,
      12.28, 5.6, 5.38, 6.6, 8.74)

CV <- function(x){
  sqrt(var(x))/mean(x)
}


jack <- numeric(length(x)-1)
pseudo <- numeric(length(x))
for (i in 1:length(x)){
  jack<-x[-i]
  pseudo[i] <- length(x)*CV(x) -(length(x)-1)*CV(jack)
}

mean(pseudo)
var(pseudo)
hist(pseudo)
var(pseudo)/length(x)
mean(pseudo) + qt(0.975,length(x)-1)*sqrt(var(pseudo)/length(x))
mean(pseudo) - qt(0.975,length(x)-1)*sqrt(var(pseudo)/length(x))