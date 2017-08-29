# Model: 
# max_y psi1*log(y+1) + psi0(I-y)
# s.t. y>=0 and I-y>=0

# psi1 = exp(x*beta + epsilon1)
# psi0 = exp(epsilon0)
# where epsilon1 and epsilon0 are iid Gumbel distribution

library(evd)
library(maxLik)
set.seed(66)

# Set parameters 
beta	<- c(-1, .4)

# Simualte data 
N		<- 2000						# Number of observations
I 		<- rep(100, N)				# Income
X		<- cbind(1, runif(N, -1, 1))
eps		<- matrix(rgev(2*N), N, 2)
psi1	<- exp(X%*%beta + eps[,1])
psi0	<- exp(eps[,2])

exp.fn	<- function(psi1, psi0, I){
# This function returns optimal expenditure y 	
	y			<- psi1/psi0 - 1		# y* = psi1/psi0 - 1 from FOC
	sel.z		<- y < 0 
	y[sel.z]	<-0						# y* = 0 if psi1 < psi0
	sel.c		<- y > I				
	y[sel.c]	<- I[sel.c]				# y* = I if y*> I
	return(y)
}

y	<- exp.fn(psi1, psi0, I)
table(y==0)
table(y==I)							# A few observations hit income
# hist(y)

# Estimation 
ll.fn	<- function(beta, y, X){
# This function returns log likelihood

# Write v = -X1*beta + log(y+1)
# P(y = 0) = P(epsilon1 - epsilon0 < -X1*beta + log(y+1)) = exp(v)/(1+exp(v))
# For y> 0, density of y is 
# g(y)  = |J|f(-X1*beta + log(y+1)) = 1/(y+1)*f(v) = 1/(y+1)* exp(-v)/(1 + exp(-v))^2, 
# where f(.)  is the density of logistic distribution
	xbeta	<- X %*% beta
	v		<- -xbeta + log(y+1)
	l0		<- v - log(1 + exp(v))					# log(p(y=0))
	# l0		<- -log(1 + exp(v))					# log(p(y=0))
	lJ		<- -log(1 + y)							# Jacobean
	lp		<- lJ - v - 2*log(1 + exp(-v))			# log(g(y))
	ll		<- 1*(y==0)*l0 + 1*(y>0)*lp
	return(ll)
}

system.time(tmp <- ll.fn(beta, y, X))
beta.init	<- c(-5, 1)
sol		<- maxLik(ll.fn, start = beta.init, y = y, X = X)
summary(sol)
cbind(beta, coef(sol))

