# Set initial values
S0 <- 12; r <- 0.02; v <- 0.3; T <- 0.65
M <- 365; dt <- T/(M-1)  #???M??? no longer = no.days

# Set seed for reproducibility and generate random numbers
set.seed(1)
u <- rnorm(M+1)

# Calculate values of Z using a loop
Z <- S0
for (m in 1:(M-1)) {
  Z <- c(Z, Z[m] + r*Z[m]*dt + v*Z[m]*sqrt(dt)*u[m+1])
}

# Generate sequence of time values and plot Z against t
t <- seq(from = 0, to = T, by = dt)
plot(t, Z, type = "l", xlab = "t", ylab = "Z(t)", sub = "set.seed(1)")
