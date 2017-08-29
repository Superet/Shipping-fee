# Model: 
# max_y psi1*log(y+1) + psi0(I-y-c_ship(y)) - c_search(y)
# s.t. y>=0 and I-y>=0
# c_ship(y) = I(y < MO)*c_f
# c_search(y)  = I(y>=MO)*alpha/(y-MO+1)

# psi1 = exp(x*beta + epsilon1)
# psi0 = exp(epsilon0)
# where epsilon1 and epsilon0 are iid Gumbel distribution

# This script does not include covariates in x

library(evd)
library(maxLik)
library(nloptr)
library(ggplot2)
library(fastGHQuad)

set.seed(66)

# Set parameters 
beta	<- 2
I 		<- 100					
MO		<- 20
cf		<- 8
alpha.mu	<- .8						# Coefficients search cost
lalpha.mu	<- log(alpha.mu)			# alpha_h ~ lognormal(mu, sigma)
lalpha.sigma<- .3
nhh		<- 100							# Number of households 
nT		<- 50							# Number of time periods per household
N		<- nhh*nT						# Number of observations
hh.idx	<- rep(1:nhh, each = nT)		# Index of households
eps		<- matrix(rgev(N*2), N, 2)		# Randon epsilon draws for psi
lalpha	<- rlnorm(nhh, lalpha.mu, lalpha.sigma) 					
alpha	<- exp(lalpha)[hh.idx]
summary(alpha)
psi1	<- exp(beta + eps[,1])
psi0	<- exp(eps[,2])

# ---------------------------- #
# Simulate optimal expenditure #
#### Functions ####
# Utility funciton 
u		<- function(y, psi1, psi0, alpha, MO, cf, I){
	out <- psi1*log(y+1) + psi0*(I- y - 1*(y>0 & y<MO)*cf) - 1*(y>=MO)*alpha/(y - MO + 1)
	return(out)
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
# psi1, psi0, and MO are vectors of parameters of the same dimension

	# Find y* 
	y1		<- y1star(psi1, psi0, I)
	y1[y1<0]<- 0

	# Find y** only for psi1 > psi0
	sel		<- which(y1 > 0)					
	y2		<- rep(0, length(y1))
	for(i in 1:length(sel)){
		y2[sel[i]]	<- y2star(psi1[sel[i]], psi0[sel[i]], alpha[sel[i]], MO, I)
	}
	
	# Calculate utility at critical values -- y*, y** and MO
	u3.y1	<- u(y1, psi1, psi0, alpha, MO, cf, I)
	u3.MO	<- u(MO, psi1, psi0, alpha, MO, cf, I)
	u3.y2	<- u(y2, psi1, psi0, alpha, MO, cf, I)

	# Optimal expenditure
	y		<- y1	
	y[u3.MO >= u3.y1 & u3.MO >= u3.y2]	<- MO
	y[u3.y2	>= u3.y1 & u3.y2 >= u3.MO]	<- y2[u3.y2	>= u3.y1 & u3.y2 >= u3.MO]
		
	if(return.u){
		maxu	<- u(y, psi1, psi0, alpha, MO, cf, I)			# Utilty at optimal value
		return(cbind(y = y, u = maxu, y1 = y1, y2 = y2))
	}else{
		return(y)
	}	
}

# Simulate data 
y 	<- exp.fn(psi1, psi0, alpha, MO, cf, I, FALSE)
hist(y, breaks = 50)
table(y == 0)
table(y == MO)
table(y > MO)

#################
# Simulated MLE # 
#################
K			<- 5000									# Number of simulation draws
eps.draw	<- matrix(rgev(K*2), K, 2)				# Draws for estimation
num.nodes	<- 8

SimDist.fn	<- function(beta, alpha, esp.draw){	
	# Given the parameters, simulate optimal expenditure
	psi1.draw	<- exp(beta + eps.draw[,1])
	psi0.draw	<- exp(eps.draw[,2])
	y			<- exp.fn(psi1.draw, psi0.draw, alpha, MO, cf, I)
	
	# Derive the empirical distribution 
	p0		<- mean(y == 0)							# P(y == 0 )
	p.mo	<- mean(y == MO)						# P(y == MO)
	g		<- try(density(y[y>0 & y!= MO], from = 0, to = I), silent = T) # g(y) for continuous y
	if(class(g) == "try-error"){
		# When no values are such that y>0 & y!= MO, then we have g(y)
		g		<- function(y){0}
 	}else{
		g$y		<- g$y * (1 - p0 - p.mo)				# For continuous values, g(y) = g(y|y>0 & y!= MO)*P(y>0 & y!= MO)
		g		<- approxfun(g)							# Linear interpolation											# Directly write the Kernel density function
	}
	return(list(p0 = p0, p.mo = p.mo, g = g))
}

ll.fn <- function(param, eps.draw){
	beta	<- param[1]
	lalpha.mu	<- param[2]
	lalpha.sigma<- exp(param[3])

	# Quadrature nodes 
	q.grid	<- gaussHermiteData(num.nodes)$x
	q.grid	<- sqrt(2)*lalpha.sigma*q.grid  + lalpha.mu
	q.wt	<- gaussHermiteData(num.nodes)$w
	
	# Simulate expenditure at each nodes
	p0		<- p.mo	<- rep(NA, num.nodes)
	g		<- vector("list", length = num.nodes)
	es		<- 1e-8
	for(i in 1:num.nodes){
		dist0	<- SimDist.fn(beta, rep(exp(q.grid[i]), K), eps.draw)
		p0[i]	<- dist0$p0
		p.mo[i] <- ifelse(dist0$p.mo==0 | is.na(dist0$p.mo), es, dist0$p.mo)
		g[[i]]	<- dist0$g
	}
		
	# Calculate likelihood for each household 
	# LL = sum_h log( E[p(y_{ht})] ), where E is taken over alpha
	# Given log(alpha) ~ Normal(mu, sigma), and Gaussian-Hermite quadrature notes nd_j and weight wt_j,
	# 		E[p(y_{ht})] = 1/sqrt(pi)*sum_j wt_j*(prod_t p(y_{ht}| sqrt(2)*sigma*nd_j+mu)
	ll.hh		<- rep(NA, nhh)
	for(i in 1:nhh){
		sel		<- hh.idx == i
		phh.grid	<- matrix(NA, sum(sel), num.nodes)
		for(j in 1:num.nodes){
			p.cont	<- g[[j]](y[sel])
			p.cont	<- ifelse(p.cont == 0 | is.na(p.cont), es, p.cont)
			phh.grid[,j]	<- 	1*(y[sel] == 0)*p0[j] + 
								1*(y[sel] == MO)*p.mo[j] + 
								1*(y[sel]>0 & y[sel]!= MO)*p.cont
		}
		Ephh		<- 1/sqrt(pi)*sum(q.wt * apply(phh.grid, 2, prod))
		ll.hh[i]	<- log(Ephh)
	}
	
	# Integrate the likelihood over the quadrature grids
	ll	<- sum(ll.hh)
	
	return(ll)
}

# NOTE: the estimation is VERY sensitive to initial values if K is not large enough
max.method 	<- "NM"
beta.init	<- c(beta, lalpha.mu, log(lalpha.sigma))
system.time(tmp <- ll.fn(beta.init, eps.draw))
sol		<- maxLik(ll.fn, start = beta.init, eps.draw= eps.draw, method = max.method)
summary(sol)
matrix(cbind(c(beta, lalpha.mu, log(lalpha.sigma)), coef(sol)), 3,2, dimnames = list(c("beta", "lalpha.mu", "lalpha.sigma"), c("true", "estimate")))

beta.init	<- c(1.5, log(.9), .1)
sol1		<- maxLik(ll.fn, start = beta.init, eps.draw= eps.draw, method = max.method)
summary(sol1)
matrix(cbind(c(beta, lalpha.mu, lalpha.sigma), coef(sol1)), 3,2, dimnames = list(c("beta", "lalpha.mu", "lalpha.sigma"), c("true", "estimate")))


# See the shape of log likelihood function 
tmp	<- seq(-1, 1, .1)
ggtmp	<- data.frame(lalpha.mu = lalpha.mu + tmp)				
ggtmp$ll<- sapply(ggtmp$lalpha.mu, function(x) ll.fn(c(beta, x, lalpha.sigma), eps.draw) )
ggtmp$true	<- lalpha.mu

ggplot(ggtmp, aes(lalpha.mu, ll)) + geom_line() + geom_point() +
 		geom_vline(aes(xintercept = true), col = "red", size = .25, linetype = 2) 
# 		facet_wrap(~ vary.var, scales = "free_x")

