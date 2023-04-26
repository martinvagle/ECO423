set.seed(1) # Groups 1, 2, 3 use seeds 1, 10, 100
S0 = 100; g = 0.05; theta = 110; v = 0.5; cy = -0.1
r = 0.05; K = 110; T = 1.5
states = 2; periods = 2; dt = T/periods
z = matrix(rnorm(states*periods), ncol=periods, byrow=T)
z
# Q-dynamics first state
WT = apply(z, 1, sum) * sqrt(dt)
WT
ST = S0 * exp((r-cy-v^2/2) * T + v*WT)
ST
# MC estimate of present value
gT = pmax(K - ST, 0); gT
p0 = mean(gT) * exp(-r * T);p0

