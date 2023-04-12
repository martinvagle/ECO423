
#Initilize the parameters
S0=100; mu=0.1; sigma =0.25; N=1000; K=165
T=5; dt = T/N; r =0.05

#I set the seed to 1 to predetermine the random numbers, X is 1000 random numbers drawn with mean 0 sigma = 1?
#Then append to the start of the vector a 0 and a 0
set.seed( 1 )
X = rnorm(N)
View(X)
X = append(X, 0 , 0 )
View(X)

#Take the cumulative sum of the vector X and multiply by the square root of dt
#This is the Brownian motion, or "Standard" wiener process
Wt = cumsum(X)*sqrt(dt)

#We initiate a sequence from 0 to TIME T with a step size of dt
t = seq(0,T,dt)

#Here we calculate the 95% confidence interval for each simulated path
lower = -1.96*sqrt(t) 
upper = 1.96*sqrt(t)
lower
upper

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
lSt
uSt

#Want to make a graph of stock price path, expected value and confidence interval
plot(t, St, type="l", ylim = c(46,422))
lines(t, ESt)   #Expected value
lines(t, lSt)   #Lower confidence interval
lines(t, uSt)   #Upper confidence interval
#In one instance stock price above upper confidence interval, very interesting

#This is the black scholes formula, why T-t?
d1 = (log(St/K)+(r+sigma^2/2)*(T-t)/(sigma*sqrt(T-t)))
d1[1:10]
d1[N]
d1[N+1]

#Want to compute N(d1) for all dates but the last one
#Pnorm fetches the likelihood for a particular observation or lower
#DNORM the cumulative distribution function for a particular observation 
#for the standard normal distribution at value d1
#Do this for all values of d1

Delta = pnorm(d1[1:N])
Delta

#Appends 1 to the end of the vector if the stock price is above the strike price
#If not append 0
Delta = append(Delta, if(St[N+1]>=K) 1 else 0)
Delta

#Value of stock position in replicating portfolio
#equals stock price times amount of stock held
stockVal = St*Delta

#Find position in riskless by computing BSM date 0 
#S0 multiplied 
#Because i am waking up at day 0, and want to simulate,
#What the value of the porftolio would be at day T
BSM = S0*pnorm(d1[1])-exp(-r*T)*K*pnorm(d1[1]-sigma*sqrt(T))
BSM

#to compute the position in the riskless asset at an arbitraty date u, we’ll define a function rlp(u),
#that recursively computes the position in the riskless asset at date u, given the position at date u-1. 
#The function rlp(u) is defined as follows:
#(recursive) function in R, using the function command
rlp = function (u) {
    if (u == 1) return (BSM-Delta[1]*St[1])
    else return (rlp(u-1)*exp(r*dt)-(Delta[u]-Delta[u-1])*St[u])
}