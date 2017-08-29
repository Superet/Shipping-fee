# Model: 
# max_y psi1*log(y+1) + psi0(I-y-c_ship(y)) - c_search(y)
# s.t. y>=0 and I-y>=0
# c_ship(y) = I(y < MO)*c_f
# c_search(y)  = I(y>=MO)*alpha/(y-MO+1)

# psi1 = exp(x*beta + epsilon1)
# psi0 = exp(epsilon0)
# where epsilon1 and epsilon0 are iid Gumbel distribution

# This script includes covariates X

library(evd)
library(maxLik)
library(nloptr)
library(np)
library(reshape2)
library(mgcv)
library(hdrcde)
set.seed(66)

# Set parameters 
beta	<- c(-1.5, .5)
alpha	<- 1				
N		<- 2000							# Number of observations
I 		<- 100					
MO		<- 20
cf		<- 8
X		<- cbind(1, runif(N))
eps		<- matrix(rgev(N*2), N, 2)		# Randon epsilon draws for psi
psi1	<- exp(X%*%beta + eps[,1])
psi0	<- exp(eps[,2])
nx		<- length(beta)

# ---------------------------- #
# Simulate optimal expenditure #

#### Functions ####
# Utility funciton 
u		<- function(y, psi1, psi0, alpha, MO, cf, I){
	psi1*log(y+1) + psi0*(I- y - 1*(y<MO)*cf) - 1*(y>=MO)*alpha/(y - MO + 1)
}

# y* = psi1/psi0 - 1 from FOC and KKT conditions
y1star	<- function(psi1, psi0, I){ 
	out	<- psi1/psi0 - 1 
	return(ifelse(out > I, I, out))
}

# y** is such that psi1/(y**+1) - psi0 + alpha/(y** - MO+ 1)^2 = 0 for y** >= MO
y2star	<- function(psi1, psi0, alpha, MO, I){
	foc <- function(y){psi1/(y+1) - psi0 + alpha/(y - MO+ 1)^2 }
	sol	<- try(uniroot(foc, lower = MO - 1e-2, upper = I+1e-2), silent = T)		# Solve for the univariate equation
	if(class(sol) == "try-error"){											# When sign(foc(MO)) = sign(foc(I)), it does not give a solution
		u1	<- u(MO, psi1, psi0, alpha, MO, cf, I)							# Then the utility is maximized at the two corners.
		u2	<- u(I, psi1, psi0, alpha, MO, cf, I)
		out	<- ifelse(u1>u2, MO, I)
	}else{
		out <- sol$root
	}
	return(out)
}

# Find optimal expenditure 
# This funciton gives the optimal expenditure based on our efficient algorithm
exp.fn	<- function(psi1, psi0, alpha, MO, cf, I, return.u = FALSE){
	# Find y* 
	y1		<- y1star(psi1, psi0, I)

	# Find y** only for psi1 > psi0
	sel		<- which(y1 > 0)					
	y2		<- rep(0, length(y1))
	for(i in 1:length(sel)){
		y2[sel[i]]	<- y2star(psi1[sel[i]], psi0[sel[i]], alpha, MO, I)
	}
	
	# Calculate utility at critical values -- y*, y** and MO
	u3.y1	<- u(y1, psi1, psi0, alpha, MO, cf, I)
	u3.MO	<- u(MO, psi1, psi0, alpha, MO, cf, I)
	u3.y2	<- u(y2, psi1, psi0, alpha, MO, cf, I)

	# Optimal expenditure
	y		<- y1
	y[y1<0]	<- 0	
	y[u3.MO >= u3.y1 & u3.MO >= u3.y2]	<- MO
	y[u3.y2	>= u3.y1 & u3.y2 >= u3.MO]	<- y2[u3.y2	>= u3.y1 & u3.y2 >= u3.MO]
		
	if(return.u){
		maxu	<- u(y, psi1, psi0, alpha, MO, cf, I)			# Utilty at optimal value
		return(cbind(y = y, u = maxu, y1 = y1, y2 = y2))
	}else{
		return(y)
	}	
}

### Test simulation ### 
y 	<- exp.fn(psi1, psi0, alpha, MO, cf, I)

# Check the simulated data
hist(y)
table(y == 0)
table(y == MO)
table(y > MO)

#################
# Simulated MLE # 
#################
K			<- 500									# Number of simulation draws
eps.draw	<- matrix(rgev(K*2), K, 2)				# Draws for estimation

# ------------------------------------------------------#
# Let d = X*beta
# Simulate the distribution of (y|d, alpha)
# This step can be slow, but we only need to simulate once 
d.seq 		<- seq(-5, 5, length = 15)
alpha.seq 	<- seq(.1, 3, length = 15)
par.grid	<- expand.grid(d.seq, alpha.seq); colnames(par.grid)	<- c("d", "alpha")
ysim		<- matrix(NA, K, nrow(par.grid))
for(i in 1:nrow(par.grid)){
	# Given the parameters d and alpha, simulate optimal expenditure
	psi1.draw	<- exp(par.grid[i,1] + eps.draw[,1])
	psi0.draw	<- exp(eps.draw[,2])
	ysim[,i]	<- exp.fn(psi1.draw, psi0.draw, par.grid[i,2], MO, cf, I)
}
summary(c(ysim))

# --------------------------------------- #
# Model the conditional distribution of y #
# Bandwidth selection for conditional distribution takes time. 
ysim.df		<- melt(ysim)
colnames(ysim.df)	<-	c("draw.idx", "par.idx", "value")
ysim.df		<- cbind(ysim.df, par.grid[ysim.df$par.idx,])
ysim.df$y0	<- 1*(ysim.df$value == 0)
ysim.df$y.mo<- 1*(ysim.df$value == MO)

# Conditional probabity of y == 0
# We model logit(P(y==0)) = cbind(1, d, alpha) %*% Delta1
p0.mod		<- glm(y0 ~ d + alpha, data = ysim.df, family = binomial(link = "logit"))
summary(p0.mod)
Delta1		<- coef(p0.mod)
p0.fn		<- function(d, alpha){
	xbeta	<- cbind(1, d, alpha) %*% Delta1
	out		<- exp(xbeta)/(1 + exp(xbeta))
	return(out)
} 

# Conditional probability of y == MO
# We model logit(P(y==MO)) = cbind(1, d, alpha) %*% Delta2
pmo.mod		<- glm(y.mo ~ d + alpha, data = ysim.df, family = binomial(link = "logit"))
summary(pmo.mod)
Delta2		<- coef(pmo.mod)
pmo.fn		<- function(d, alpha){
	xbeta	<- cbind(1, d, alpha) %*% Delta2
	out		<- exp(xbeta)/(1 + exp(xbeta))
	return(out)
}

# Conditional density of (y|d, alpha) for y > 0 and y!= MO
sel		<- ysim.df$value > 0 & ysim.df$value != MO
ysimp.df<- ysim.df[sel,]
fhat	<- cde(x = as.matrix(ysimp.df[,c("d", "alpha")]), y = ysimp.df[,"value"])
plot(fhat)
# plot(fhat, plot.fn = "hdr")

fhat.fn	<- function(ynew, xnew, ysimp.df, a = fhat$a, b = fhat$b){
# Conditional density estimation in Equation (1.1) from Bashtannyk and Hyndman (2001)
	N		<- length(ynew)
	n		<- nrow(ysimp.df)
	x0		<- as.matrix(ysimp.df[,c("d", "alpha")])
	xtb		<- table(ysimp.df$par.idx)
	x00		<- par.grid[as.numeric(names(xtb)),]
	nux		<- length(xtb)
	y0		<- ysimp.df[,"value"]
	out		<- rep(NA, N)
	for(i in 1:N){
		kn.x	<- dnorm( sqrt(rowSums( (matrix(rep(xnew[i,], nux), ncol = 2, byrow = T) - x00)^2))/a ) 
		f.x		<- sum(xtb*kn.x)/(n*a)
		kn.x	<- kn.x[ysimp.df[,"par.idx"]]
		kn.y	<- dnorm( (ynew[i] - y0)/b )
		f.xy	<- mean(kn.x*kn.y)/(a*b)
		out[i]	<- f.xy/f.x 
	}
	
	# We need to multiple (1-p0-p.mo) so that the probability across the entire region sums up to 1. 
	out		<- out * (1 - p0.fn(xnew[,1], xnew[,2]) - pmo.fn(xnew[,1], xnew[,2]) )
	return(out)
}
# Test function
xnew	<- cbind(c(-2, 1, 0, 1), c(1, 1.2, 1.5, 1))
fhat.fn(rep(5, 4), xnew, ysimp.df)

# Conditional density estimation from package np has better bandwidth selection, but is TOO slow. 
# np.fhat		<- npcdens(value ~ d + alpha, data = ysimp.df)

# ------------------------#
# Log likelihood function # 
ll.fn <- function(param){
	beta	<- param[1:nx]
	alpha	<- exp(param[length(param)])
	
	### Simulation step ###
	# Given the parameters, simulate optimal expenditure y and have the distribution 
	d		<- X %*% beta
	p0		<- p0.fn(d, alpha)
	p.mo	<- pmo.fn(d, alpha)
	g		<- fhat.fn(y, cbind(d, alpha), ysimp.df)
	e		<- -100									# Small negative values for not well-define or -Inf when taking log
	
	### Evaluation step ###
	# Calculate the log likelihood function
	l0 		<- log(p0)
	lp.mo	<- log(p.mo)
	lp.mo	<- ifelse(lp.mo == -Inf | is.na(lp.mo), e, lp.mo)
	l.cont	<- log(g)
	l.cont	<- ifelse(l.cont == -Inf | is.na(l.cont), e, l.cont)
	ll		<- 	1*(y == 0)*l0 + 
				1*(y==MO)*lp.mo + 
				1*(y>0 & y!= MO)*l.cont
	return(ll)
}

# Likelihood evaluation is slow...
system.time(tmp <- ll.fn(c(beta, log(alpha))))
beta.init	<- c(beta, log(alpha))
sol		<- maxLik(ll.fn, start = beta.init, method = "NM")
summary(sol)
matrix(cbind(c(beta, log(alpha)), coef(sol)), 3,2, dimnames = list(c("beta1", "beta2","log(alpha)"), c("true", "estimate")))

# beta.init	<- c(beta - .2, log(.9))
# sol		<- maxLik(ll.fn, start = beta.init, method = "NM")
# summary(sol)
# matrix(cbind(c(beta, log(alpha)), coef(sol)), 3,2, dimnames = list(c("beta1", "beta2","log(alpha)"), c("true", "estimate")))

# Check the shape of log likelihood function
# Fix alpha
tmp	<- seq(-1, 1, .05) + beta[2]
tmp1	<- sapply(tmp, function(x) sum(ll.fn(c(beta[1], x, log(alpha)))))
plot(tmp, tmp1, xlab = "beta2", ylab = "Loglikelihood")
lines(tmp, tmp1)
abline(v = beta, col = "red")

# Fix beta
tmp	<- seq(-1, 1, .05) + log(alpha)
tmp1	<- sapply(tmp, function(x) sum(ll.fn(c(beta, x))))
plot(tmp, tmp1, xlab = "log(alpha)", ylab = "Loglikelihood")
lines(tmp, tmp1)
abline(v = log(alpha), col = "red")



