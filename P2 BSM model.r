#Initilize the parameters
S0=100; mu=0.1; sigma =0.25; N=1000; K=165 #K 
T=5; dt = T/N; r =0.05


#Then append to the start of the vector a 0 and a 0
set.seed( 1 )
X = rnorm(N) #random numbers drawn from normal distribution rnorm(n, mean = 0, sd = 1)

X = append(X, 0 , 0 )

#Take the cumulative sum of the vector X and multiply by the square root of dt
#This is the Brownian motion, or "Standard" wiener process
Wt = cumsum(X)*sqrt(dt)

#We initiate a sequence of "dates" from 0 to TIME T with a step size of dt
t = seq(0,T,dt)

#Here we calculate the 95% confidence interval for each simulated path
lower = -1.96*sqrt(t); upper = 1.96*sqrt(t)

#We plot the simulated paths WT, for the different time intervals
plot(t, Wt, type="l", ylim = c(-4.5 ,4.5))
#Lines that show confidence interval for different vlaues of Wt
lines(t, lower)
lines(t, upper)

#Now we simulate the stockprice S(t) for each time interval
St = S0*exp((mu-0.5*sigma^2)*t + sigma*Wt)

#The expected value of the stock price at each future date is
#IE mean value
ESt = S0*exp(mu*t)

#We can generate 95% confidence intervals for the stock price process:
#First we have the beginning stockprice and drift, then we distort the Wiener process by the amount corresponding to 95% confidence interval
lSt = S0*exp((mu-0.5*sigma^2)*t - 1.96*sqrt(t)*sigma)
uSt = S0*exp((mu-0.5*sigma^2)*t + 1.96*sqrt(t)*sigma)

#Want to make a graph of stock price path, expected value and confidence interval
plot(t, St, type="l", ylim = c(46,422))
lines(t, ESt)   #Expected value
lines(t, lSt)   #Lower confidence interval
lines(t, uSt)   #Upper confidence interval

#This is the black scholes formula, why T-t?
d1 = (log(St/K)+(r+sigma^2/2)*(T-t)/(sigma*sqrt(T-t)))
d2 = d1-sigma*sqrt(T)

#Want to compute N(d1) for all dates but the last one
Delta = pnorm(d1[1:N])
#Appends 1 to the end of the vector if the stock price is above the strike price
#If not append 0
Delta = append(Delta, if(St[N+1]>=K) 1 else 0)

# In summary, the dnorm() function calculates the probability density at a specific point or set of points, 
# while the pnorm() function calculates the probability that a random variable from the distribution is
# less than or equal to a specific point or set of points.

#Value of stock position in replicating portfolio
#equals stock price times amount of stock held
stockVal = St*Delta

#Find position in riskless by computing BSM at date 0 
#S0 multiplied 
#Because i am waking up at day 0, and want to simulate,
#What the value of the porftolio would be at day T. this is because the decision is made in day 0
BSM = S0*pnorm(d1[1])-exp(-r*T)*K*pnorm(d2[1])

#Compute the position in the riskless asset at an arbitraty date u:
rlp = function(u) {
  if (u==1) return (BSM-Delta[1]*St[1])
  else return (rlp(u-1)*exp(r*dt)-(Delta[u] - Delta [u-1])*St[u])
}

#Value of the replicating portfolio
rp = function(u) return (Delta[u]*St[u] + rlp(u))

rp(N+1) #Value day T replicating portfolio
max(St[N+1]-K,0) #Value day T call option
rp(N+1)-max(St[N+1]-K,0) #Diff



