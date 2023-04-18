### ECO423 Advanced Derivatives - Greeks and discrete, dynamic delta hedging

### Illustrations of Greeks
# v = volatility

# Graph of delta, non-dividend paying asset
delta = function(x,r,v,T)
{
   d = (log(x)+(r+v^2/2)*T)/v/sqrt(T)
   return(pnorm(d))
}
curve(delta(x,.02,.4,1.0), 0.2, 1.8, xlab="S/K", ylab="delta")

# Graph of vega, non-dividend paying asset
vega = function(x,K,r,v,T)
{
   d = (log(x)+(r+v^2/2)*T)/v/sqrt(T)
   
   return(x*K*dnorm(d)*sqrt(T))
}
curve(vega(x,10,.02,.4,1.0), 0.2, 1.8, xlab="S/K", ylab="vega")

# Graph of theta, non-dividend paying asset
theta = function(x,r,v,T)
{
   d = (log(x)+(r+v^2/2)*T)/v/sqrt(T)
   N2 = pnorm(d-v*sqrt(T))
   
   return(-x*dnorm(d)*v/2/sqrt(T) - r*exp(-r*T)*N2)
}
curve(theta(x,.02,.4,1.0), 0.2, 1.8, xlab="S/K", ylab="theta")
curve(theta(1.0,.02,.4,x), 0.05, 5.0, xlab="T", ylab="theta")

### Analytic versus numeric Greeks
# Analytic delta
S = 100; K = 90; r = 0.05; v = 0.4; T = 0.5
d = (log(S/K)+(r+v^2/2)*T)/v/sqrt(T)
N1 = pnorm(d); N1

# Numeric delta
N2 = pnorm(d-v*sqrt(T))
c = S*N1-K*exp(-r*T)*N2

h = 10.0 # 10 bad, 1 decent, 0.1 great
dh = (log((S+h)/K)+(r+v^2/2)*T)/v/sqrt(T)
N1h = pnorm(dh)
N2h = pnorm(dh-v*sqrt(T))
ch = (S+h)*N1h-K*exp(-r*T)*N2h

(ch-c)/h

### Option beta ~ option elasticity
# Graph of elasticity, to illustrate variation in beta
ela = function(x,r,v,T)
{
   d = (log(x)+(r+v^2/2)*T)/v/sqrt(T)
   N1 = pnorm(d)
   N2 = pnorm(d-v*sqrt(T))
   c = x*N1-exp(-r*T)*N2
   
   return(N1/c)	
}
curve(ela(x,.02,.4,1.0), 0.2, 1.8, xlab="S/K", ylab="elasticity")

   
### ------------------- Effectiveness of dynamic hedging -------------------
BSM = function(S,K,r,sigma,T) {
   d = (log(S/K)+(r+sigma^2/2)*T)/sqrt(T)/sigma
   return(S*pnorm(d) - K*exp(-r*T)*pnorm(d-sqrt(T)*sigma))
}

# Parameter values from Hull's example in Section 19.1 (11e)
S0 = 49; K = 50; r = 0.05; sigma = 0.2; mu = 0.13

cost = function(weeks, M, N) {
   # 'weeks' to expiration, 'M' periods, 'N' paths
   T = weeks/52; dt = T/M
   date = dt*c(0:M) # M+1 dates 0, dt, 2*dt, ..., M*dt = T
   
   dW = matrix(rnorm(M*N, sd = sqrt(dt)), ncol = M)
   W = cbind(0,t(apply(dW, 1, cumsum))) # Must transpose apply (cf 'Value' in documentation)
   S = matrix(-1, nrow = N, ncol = M+1) # Initialize to -1, for easy error checking
   for(n in 1:N)
      S[n,] = S0*exp((mu-sigma^2/2)*date + sigma*W[n,]) # Compute price path n
   Delta = matrix(rep(-1, N*(M+1)), nrow = N) # Initialize to -1 as above
   for(n in 1:N) {
      Delta[n,] = delta(S[n,]/K,r,sigma,T-date) # Compute delta path n
   }
   
   # Store cost across all dates and states, in case I want to e.g. graph paths (which 
   # I don't do in this particular implementation; not a memory efficient approach)
   cost = matrix(rep(-1, N*(M+1)), nrow = N)
   cost[,1] = Delta[,1]*S[,1]
   for(m in 2:(M+1)) {
      cost[,m] = (Delta[,m]-Delta[,m-1])*S[,m] + cost[,m-1]*exp(r*dt)
   }
   
   # Reduce cost by receipt of strike price if ITM
   for(n in 1:N) {
      if(S[n,M+1] > K)
         cost[n,M+1] = cost[n,M+1] - K
   }
   
   cost.avg = exp(-r*T)*mean(cost[,M+1])
   cost.sd = exp(-r*T)*sd(cost[,M+1])
   
   bsm.ce.value = BSM(S0,K,r,sigma,T)
   
   return(c(bsm.ce.value, cost.avg, cost.sd/bsm.T.value))
}

# 20 weeks to expiration, as in Hull's example
set.seed(1)
res5 = cost(20,4,1000000)
res4 = cost(20,5,1000000)
res2 = cost(20,10,1000000)
res1 = cost(20,20,1000000)
res05 = cost(20,40,1000000)
res025 = cost(20,80,1000000)

res5 # etc to report results