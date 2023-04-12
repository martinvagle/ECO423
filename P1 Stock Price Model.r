# 1. Initilize the parameters
S0 <- 100
mu <- 0.1
sigma <- 0.25
T <- 5
N <- 5000

# 2. Box-Muller transform
set.seed(1)
U1 <- runif(N) # vector of N uniform random numbers
U2 <- runif(N) # vector of N uniform random numbers
X1 <- sqrt(-2 * log(U1)) * cos(2 * pi * U2)
X2 <- sqrt(-2 * log(U1)) * sin(2 * pi * U2)

# 3. Histogram
hist(X1)
hist(X2)

# 4. Compare emirical density with normal density
hist(X1, 10/0.25, probability = TRUE) #Normalizing histogram
lines(x1, dnorm(x1)) #implements normal density function

#Find sample mean and var
mean(X1); var(X1)

# 6. Simulate the future stock prices
ST <- S0 * exp((mu - sigma^2/2) * T + sigma * sqrt(T) * X1) 
#Since X1 is a vector ST becomes the reult of the function on the vector
ST[1:10] #Peek

# 7. Compute sample mean and standard deviation
mean(ST)
sd(ST)

# 8. compare the stats with theoretical/true values
mST <- S0 * exp(mu * T); mST #Expected S at time T
vST <- S0^2 * exp(2 * mu * T) * (exp(sigma^2 * T) - 1) #Var
sdST <- sqrt(vST); sdST #sd

# 9. compare empirical density function of St
summary(ST)
hist(ST) #Lognormal behavior
hist(ST, 100, probability = TRUE)
s1 <- seq(0, 800, 8) #caps observations at 0-800
points(s1, dnorm(s1, mST, sdST)) #will be wrong comparison

hist(ST, 100, probability = TRUE) #Refresh chart
#Calculates lognormal density function
s <- sqrt(log(vST/mST^2 + 1))
m <- log(mST) - s^2/2
points(s1, exp(-(log(s1) - m)^2/2/s^2) / s / sqrt(2 * pi) / s1) #Better fit

# 11. compute continuously compounded returns and compare
hist(log(ST) - log(S0), probability = TRUE)
mRT <- (mu-sigma^2 / 2) * T
vRT <- sigma^2 * T
points(x1, dnorm(x1, mRT, sqrt(vRT)))


# If we where to resample and edit values of sigma and mu this would happen:
# Increasing mu shifts distribution to the right, without
# changing its shape. Increasing sigma makes the tails fatter and the center values
# lower, without shifting the point of symmetry (its mean). For all values of sigma >
# 0, there are possible observations of values at the far ends of the interval for
# ’frequency’, corresponding to the normal distribution having the whole real line as
# support.

#Stock movement
plot(ST[1:365], type="l", xlab="day") 

