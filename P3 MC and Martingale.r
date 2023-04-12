# initialise the environment
rm(list = ls())

# Variables
r = 0.05
sigma = 0.25
S0 = 100
T = 5
K = 165
N = 50000

# Black Scholes components
d1 = (log(S0/K) + (r + sigma^2/2)*T)/(sigma*sqrt(T))
d2 = d1 - sigma*sqrt(T)

# Black Scholes price
BSM = S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
BSM

# setting seed for reproducibility
set.seed(1)

# generate N standard normal random variates
X = rnorm(N)

# generate N possible date T stock prices under the equivalent martingale measure
ST = S0 * exp((r - sigma^2/2) * T + sigma * sqrt(T) * X)

hist(ST, breaks = 100, col = "lightblue", 
main = "Histogram of Stock Prices at T", 
xlab = "Stock Price at T", ylab = "Frequency")

mean(ST)

# generate call option payoffs that correspond to the N possible stock prices
CT = pmax(ST - K, 0)
mean(CT)
hist(CT, breaks = 100, col = "lightblue",
main = "Histogram of Call Option Payoffs at T",
xlab = "Call Option Payoff at T", ylab = "Frequency")

# approximate expectation of the call with discounted mean of CT
BSMmc = mean(CT) * exp(-r * T)
mean(BSMmc)

# standard error of the mean
SEM = sd(CT * exp(-r * T))/sqrt(N)
SEM

# 95% confidence interval
BSMmc + c(-1, 1) * 1.96 * SEM

# assign upper limit to uci
uci = BSMmc + 1.96 * SEM

# assign lower limit to lci
lci = BSMmc - 1.96 * SEM

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
aBSMmc + c(-1, 1) * 1.96 * aSEM

# assign upper limit to uci
auci = aBSMmc + 1.96 * aSEM

# assign lower limit to lci
alci = aBSMmc - 1.96 * aSEM

auci; alci

# indicator function
if (BSM < auci & BSM > alci) 1 else 0

### Nice, checks out ###

# antithetic variate mean
BSMav = (BSMmc + aBSMmc)/2

# average option payoff
CTav = (CT + aCT)/2
head(CTav)

# standard error of the mean CTav
avSEM = sd(CTav * exp(-r * T))/sqrt(N)
avSEM

# 95% confidence interval
BSMav + c(-1, 1) * 1.96 * avSEM

# assign upper limit to uci
avuci = BSMav + 1.96 * avSEM

# assign lower limit to lci
avlci = BSMav - 1.96 * avSEM

avuci; avlci

# indicator function
if (BSM < avuci & BSM > avlci) 1 else 0

### Nice, checks out ###

# Monte Carlo simulation
# initialise the environment
rm(list = ls())
set.seed(1)

# Variables
r = 0.05
sigma = 0.25
S0 = 100
T = 5
K = 165
N = 50000
ctrBSMmc = 0 ; ctraBSMmc = 0 ; ctrBSMav = 0

# Black Scholes components
d1 = (log(S0/K) + (r + sigma^2/2)*T)/(sigma*sqrt(T))
d2 = d1 - sigma*sqrt(T)

# Black Scholes price
BSM = S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)

for (i in 1:100) {
    X = rnorm(N)
    ST = S0 * exp((r - sigma^2/2) * T + sigma * sqrt(T) * X)
    CT = pmax(ST - K, 0)
    BSMmc = mean(CT) * exp(-r * T)
    SEM = sd(CT * exp(-r * T))/sqrt(N)
    uci = BSMmc + qnorm(1 - 0.05 / 2) * SEM
    lci = BSMmc + qnorm(0.05 / 2) * SEM
    if (BSM < uci & BSM > lci) {
        ctrBSMmc = ctrBSMmc + 1
    }

    aX = -X
    aST = S0 * exp((r - sigma^2/2) * T + sigma * sqrt(T) * aX)
    aCT = pmax(aST - K, 0)
    aBSMmc = mean(aCT) * exp(-r * T)
    aSEM = sd(aCT * exp(-r * T))/sqrt(N)
    auci = aBSMmc + qnorm(1 - 0.05 / 2) * aSEM
    alci = aBSMmc + qnorm(0.05 / 2) * aSEM
    if (BSM < auci & BSM > alci) {
        ctraBSMmc = ctraBSMmc + 1
    }

    CTav = (CT + aCT)/2
    BSMav = (BSMmc + aBSMmc)/2
    avSEM = sd(CTav * exp(-r * T))/sqrt(N)
    avuci = BSMav + qnorm(1 - 0.05 / 2) * avSEM
    avlci = BSMav + qnorm(0.05 / 2) * avSEM
    if (BSM < avuci & BSM > avlci) {
        ctrBSMav = ctrBSMav + 1
    }
}

ctrBSMmc; ctraBSMmc; ctrBSMav
