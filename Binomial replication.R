#Binomial Model

#Parameters
S_0 = 9; u=1.25; d=1/u; R=1.05; p = 0.55; K = 10 #Note that p is not needed

#First we find High/Low value for put
H=max(K-S_0*u, 0); L=max(K-S_0*d,0)

#Compute some more params
delta = (H-L)/(S_0*(u-d)) #How many stocks to hold/short for negative
theta = (H-delta*u*S_0)/R #Savings in riskless rate
q = (R-d)/(u-d)

#Payoffs are the same between put and replicating portfolio
portfolio_H = delta*S_0*u+theta*R
portfolio_L = delta*S_0*d+theta*R

#Date 0 prices
portfolio_price = S_0*delta + theta
put_price = (q*H+(1-q)*L)/R

#q-expected return should equal R
S_R = (q*S_0*u+(1-q)*S_0*d)/S_0
put_R = (q*H+(1-q)*L)/put_price

#Find price off call
C_H = max(S_0*u-K, 0); C_L=max(S_0*d-K,0)
call_price = (q*C_H+(1-q)*C_L)/R #Can use same q (Not depended on direction long/short)

#Using put-call parity
p_c = S_0+put_price-K/R

