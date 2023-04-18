#Variables
r = 0.05; sigma = 0.25; S0 = 100; T = 5; K = 165; N = 50000;

# Black Scholes components
d1 = (log(S0/K) + (r + sigma^2/2)*T)/(sigma*sqrt(T))
d2 = d1 - sigma*sqrt(T)

# Black Scholes price
BSM = S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2); BSM

# generate N standard normal random variates
set.seed(1)
X = rnorm(N)

# generate N possible date T stock prices under the equivalent martingale measure
ST = S0 * exp((r - sigma^2/2) * T + sigma * sqrt(T) * X) #using r instead of mu

hist(ST, breaks = 100, col = "lightblue", 
main = "Histogram of Stock Prices at T", 
xlab = "Stock Price at T", ylab = "Frequency")

mean(ST)

# generate call option payoffs that correspond to the N possible stock prices
CT = pmax(ST - K, 0)

# approximate expectation of the call with discounted mean of CT
BSMmc = mean(CT) * exp(-r * T)
mean(BSMmc)

# standard error of the mean
SEM = sd(CT * exp(-r * T))/sqrt(N); SEM

# 95% confidence interval
ci = BSMmc + c(-1, 1) * 1.96 * SEM
# assign lower limit to lci
lci = ci[1]
# assign upper limit to uci
uci = ci[2]

uci; lci

# indicator function
if (BSM < uci & BSM > lci) 1 else 0

### Nice, checks out ###

# Antithetic Variates
aX = -X
sum(aX + X) # should be 0

# generate N possible date T stock prices under the equivalent martingale measure
aST = S0 * exp((r - sigma^2/2) * T + sigma * sqrt(T) * aX)

# generate call option payoffs that correspond to the N possible stock prices
aCT = pmax(aST - K, 0)

# approximate expectation of the call with discounted mean of CT
aBSMmc = mean(aCT) * exp(-r * T)

# standard error of the mean
aSEM = sd(aCT * exp(-r * T))/sqrt(N)

# 95% confidence interval
aci = aBSMmc + c(-1, 1) * 1.96 * aSEM
alci = aci[1]
auci = aci[2]
auci; alci

# indicator function
if (BSM < auci & BSM > alci) 1 else 0
### Nice, checks out ###

# Average of the normal and antithetic variates
BSMav = (BSMmc + aBSMmc)/2 #Price
CTav = (CT + aCT)/2 #Payoff

# standard error of the mean CTav
avSEM = sd(CTav * exp(-r * T))/sqrt(N); avSEM

# 95% confidence interval
avSD = BSMav + c(-1, 1) * 1.96 * avSEM
avuci = BSMav + 1.96 * avSEM
avlci = BSMav - 1.96 * avSEM
avuci; avlci

# indicator function
if (BSM < avuci & BSM > avlci) 1 else 0
# The % improvement in SEM from antithetic variates technique:
(SEM/sqrt(2)-avSEM) / (SEM/sqrt(2))

# Monte Carlo simulation
# Step 15
ctrBSMmc = 0; ctraBSMmc = 0; ctrBSMav = 0

for(i in 1:100) {
  X = rnorm(N)
  ST = S0*exp((r-sigma^2/2)*T + sigma*sqrt(T)*X)
  CT = pmax(ST-K,0)
  BSMmc = exp(-r*T)*mean(CT)
  sdBSMmc = sd(CT*exp(-r*T))/sqrt(N)
  uci = BSMmc + qnorm(1-0.05/2)*sdBSMmc
  lci = BSMmc + qnorm(0.05/2)*sdBSMmc
  if(BSM < uci & BSM > lci) {
    ctrBSMmc = ctrBSMmc + 1
  }
  
  aX = -X
  aST = S0*exp((r-sigma^2/2)*T + sigma*sqrt(T)*aX)
  aCT = pmax(aST-K,0)
  aBSMmc = exp(-r*T)*mean(aCT)
  asdBSMmc = sd(aCT*exp(-r*T))/sqrt(N)
  auci = aBSMmc + qnorm(1-0.05/2)*asdBSMmc
  alci = aBSMmc + qnorm(0.05/2)*asdBSMmc
  if(BSM < auci & BSM > alci) {
    ctraBSMmc = ctraBSMmc + 1
  }
  
  CTav = (CT+aCT)/2
  BSMav = exp(-r*T)*mean(CTav) # = (BSMmc+aBSMmc)/2
  sdBSMav = sd(CTav*exp(-r*T))/sqrt(N)
  uciAv = BSMav + qnorm(1-0.05/2)*sdBSMav
  lciAv = BSMav + qnorm(0.05/2)*sdBSMav
  if(BSM < uciAv & BSM > lciAv) {
    ctrBSMav = ctrBSMav + 1
  }
}

ctrBSMmc; ctraBSMmc; ctrBSMav
