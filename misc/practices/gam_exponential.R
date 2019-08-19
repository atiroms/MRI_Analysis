# Parameters
alpha <- 20
beta <- -0.05
theta <- 10

# Sample some points along x axis
n <- 100
x <- seq(n)

# Make  y = f(x) + Gaussian_noise 
data.df <- data.frame(x = x,
                      y = alpha * exp(beta * x) + theta + rnorm(n))

# plot data
plot(data.df$x, data.df$y)


# Select an approximate $\theta$, since theta must be lower than min(y), and greater than zero
#theta.0 <- min(data.df$y) * 0.5  
theta.0 <- min(data.df$y)*0.8


# Estimate the rest parameters using a linear model
model.0 <- lm(log(y - theta.0) ~ x, data=data.df)  
alpha.0 <- exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]

# Starting parameters
start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
start


#model <- nls(y ~ alpha * exp(beta * x) + theta , data = data.df, start = start)
model <- nls(y ~ alpha * exp(beta * x) + theta , data = data.df,start=list(alpha=alpha.0))

# Plot fitted curve
plot(data.df$x, data.df$y)
lines(data.df$x, predict(model, list(x = data.df$x)), col = 'skyblue', lwd = 3)

summary(model)
