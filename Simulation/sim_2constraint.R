# A simple simulation of allocation problem with 2 constraints

# Model: 
# max_b \sum_k psi_k*log(b_k + 1)
# s.t.	b_k >=0, sum_k b_k <= I, sum_{k<K} b_k >= MO

library(evd)

set.seed(111)
K 		<- 3						# Number of alternatives, assuming that the third one is outside good.
N 		<- 200						# Number of simulations
log_psi <- c(-0.2, 0.5, 0)			# Parameters in the marginal utility
log_psi	<- c(-2, -1, 0)						# Under this parameter set, we will see a larger chance that MO constriant binds
I 		<- 100						# Income constraint
MO 		<- 10						# Minimum order 
eps 	<- matrix(rgumbel(N*K, scale = 1), N, K)	# Draw random draws from extreme value distribution
psi 	<- exp(rep(1, N) %*% t(log_psi) + eps) 		# Maringal utility

# Define (negative) utility function
uf<- function(b, psi){
  -sum(psi*log(b + 1))
}

# Gradient of utility function
ugrad <- function(b, psi){
  -psi/(b+1)
}

# First simulate allocation without accounting for shipping cost 
# Define constriants: ui %*% b - ci >= 0
# 4 constraints: b_k >= 0, sum_k b_k <= I
ui		<- rbind(diag(K), c(-1,-1,-1))
ci 		<- c(rep(0, K), -I)
b0 		<- matrix(NA, N, K)
u0		<- rep(NA, N)
for(i in 1:N){
	tmp <- constrOptim(b.init, uf, grad = ugrad, ui = ui, ci= ci, psi = psi[i,])
	if(tmp$convergence != 0){
		print(i)
	}
	b0[i,] <- tmp$par
	u0[i]	<- -tmp$value
}
summary(b0)
# Sum up the total expenditure including outside good, should be exactly =  I. 
summary(rowSums(b0))
y0 <- rowSums(b0[,1:2])
summary(b0)

# Simulate allocation when using topping up strategy
# Define constriants: ui %*% b - ci >= 0
# 5 constraints: b_k >= 0, sum_k b_k <= I, b_1 + b_2 >= MO
ui.top <- rbind(diag(K), c(-1,-1,-1), c(1,1,0))
ci.top <- c(rep(0, K), -I, MO)
b.init <- c(10,5, 30)
uf(b.init, psi[1,])  
ui.top %*% b.init - ci.top

# Simulate optimal expenditure by solving the constrained optimization
b.top 	<- matrix(NA, N, K)
u.top	<- rep(NA, N)
for(i in 1:N){
  tmp <- constrOptim(b.init, uf, grad = ugrad, ui = ui.top, ci= ci.top, psi = psi[i,])
  if(tmp$convergence != 0){
    print(i)
  }
  b.top[i,] <- tmp$par
  u.top[i]	<- -tmp$value
}
summary(b.top)
# Sum up the basket expenditure excluding outside good. 
y.top <- rowSums(b.top[,1:2])
summary(y.top)

# Examine the cases where y0 < MO
sel	<- b0 < MO
summary(y.top[sel])
summary(y.top[!sel])
summary(y.top[!sel] - y0[!sel])
summary(u0[sel] - u.top[sel])
summary(u0[!sel] - u.top[!sel])

