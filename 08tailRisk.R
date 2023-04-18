# ECO423 Advanced Derivatives - Simple bootstrapped VaR
# (You must install package yfR if not already installed)
require(yfR)

# Model for future call prices
BSM = function(S,K,r,sigma,T) {
   d = (log(S/K)+(r+sigma^2/2)*T)/sqrt(T)/sigma
   return(S*pnorm(d) - K*exp(-r*T)*pnorm(d-sqrt(T)*sigma))
}

# Sample to be used in empirical distributions
first.date = "2000-01-01"
last.date = "2020-12-31" # Set to Sys.Date() if you want to use latest observations
freq.data = 'daily'
tickers = c('MSFT')      # Feel free to try Exxon, General Electric,
                         # Crude Oil Futures: 'XOM','GE','CL=F'
PriceData = yf_get(tickers = tickers, first_date = first.date,
                   last_date = last.date, freq_data = freq.data)
R = PriceData$ret_adjusted_prices[-1] # Exclude first observation, which equals NA

# Parameters; naked position in one European call
S0 = last(PriceData$price_close); v = sd(R[(length(R)-251):length(R)])*sqrt(252)
K = 220; r = 0.02; T = 0.5
# Set current price to (latest) market price; "cheating" here, using BSM price
c0 = BSM(S0,K,r,v,T)

# Simple historical VaR - ignoring autocorrelation structure
N = 10000; set.seed(1)
Rsample = sample(R, N, replace = TRUE)
summary(Rsample)

# Scenarios for call prices tomorrow, using BSM assumptions. This is the step where we need
# a model; here the BSM formula. If simulation-based model, remember to simulate under Q
# using S0*(1+Rsample) as initial conditions for the simulations.
cScenarios = BSM(S0*(1+Rsample), K, r, v, T-1/252)
hist(c0 - cScenarios)
VaR = -quantile(c0 - cScenarios, probs = 0.01); VaR

# 10-day VaR
VaR10 = VaR*sqrt(10); VaR10
