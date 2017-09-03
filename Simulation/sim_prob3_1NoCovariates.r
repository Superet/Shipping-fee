# This script tests the efficient simulation method and tests estimation 

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
set.seed(66)

# Set parameters 
beta	<- 6
alpha	<- 1500				
N		<- 3000							# Number of observations
I 		<- 250000/12					
MO		<- 1000
cf		<- 100
eps		<- matrix(rgev(N*2), N, 2)		# Randon epsilon draws for psi
psi1	<- exp(beta + eps[,1])
psi0	<- exp(eps[,2])

# ---------------------------- #
# Simulate optimal expenditure #
# NOTE: 
# We compare the optimal solution among three algorithems 
# (i) the efficient algorithm we derive from the first order condition and the KKT condition
# (ii) global optimization from the package NLOPT
# (iii) Brent opimization

# The global optimization should give us more accurate results, but it is slow. 
# The following test shows that our efficient algorithem gives pretty accurate solutions. 

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
	y1[y1<0]<- 0

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

# This function gives the optimal expenditure using global optimizor
exp.nlopt	<- function(psi1, psi0, alpha, MO, cf, I, return.u = FALSE){
	nloptr.opts	<- list("algorithm"="NLOPT_GN_DIRECT", "xtol_abs" = 1e-3, "maxeval" = 200)
	y			<- rep(NA, length(psi1))
	maxu		<- rep(NA, length(psi1))
	for(i in 1:length(psi1)){
		obj.fn	<- function(y){ - u(y, psi1[i], psi0[i], alpha, MO, cf, I)}				# Objective function
		sol		<- try(nloptr(x0= MO, eval_f = obj.fn, lb = 0, ub = I, opts = nloptr.opts), silent=TRUE)		
		y[i]	<- sol$solution
		if(y[i]< 1e-3){ y[i] <- 0}
		maxu[i]	<- sol$objective
	}
	return(cbind(y = y, u = -maxu))
}

# This function gives the optimal expenditure using Brent optimizer
exp.opt	<- function(psi1, psi0, alpha, MO, cf, I, return.u = FALSE){
	y			<- rep(NA, length(psi1))
	maxu		<- rep(NA, length(psi1))
	for(i in 1:length(psi1)){
		obj.fn	<- function(y){ - u(y, psi1[i], psi0[i], alpha, MO, cf, I)}
		sol		<- try(optimize(f = obj.fn, lower = 0, upper = I), silent=TRUE)
		if(class(sol) == "try-error"){
			y[i]	<- NA
		}else{
			y[i]	<- sol$minimum
			if(y[i]< 1e-3){ y[i] <- 0}			
			maxu[i]	<- sol$objective
		}
	}
	return(cbind(y = y, u = -maxu))
}

### Test simulation ### 
y 	<- exp.fn(psi1, psi0, alpha, MO, cf, I, TRUE)
y1	<- exp.nlopt(psi1, psi0, alpha, MO, cf, I, TRUE)
y2	<- exp.opt(psi1, psi0, alpha, MO, cf, I, TRUE)

# See if three algorithems give similar results 
maxu<- cbind(y[,"u"], y1[,"u"], y2[,"u"])
cat("Correlation of the maximized utlity:\n"); print(cor(maxu)); cat("\n")		# They all reach the same utility level 							
y	<- y[,"y"]; y1 <- y1[,"y"]; y2	<- y2[,"y"]
cat("Correlation of optimal solution:\n"); print(cor(cbind(y, y1, y2))); cat("\n")
sel <- y==0
summary(y[sel] - y1[sel]); summary(y[sel] - y2[sel])			# Have the same solution for y = 0 
sel <- y==MO													# Have the very close solution for y = MO. 
summary(y[sel] - y1[sel]); summary(y[sel] - y2[sel])
maxu[sel,]

# Plot the simulated data
par(mfrow = c(2,1))
hist(y, breaks = 100)
hist(y[y>0], breaks = 100, xlim = c(0, 5000))

table(y == 0)
table(y == MO)
table(y > MO)

# Plot utility function
tmp1	<- seq(-1, 1, .5) + beta
tmp2	<- seq(-10, 20, 5) + alpha
tmpy	<- seq(0, 2*MO, 20)
ggtmp	<- data.frame(NULL)
for(i in 1:length(tmp1)){
	for(j in 1:length(tmp2)){
		tmpu <- u(tmpy, exp(tmp1[i]), psi0 = 1, alpha = tmp2[j], MO, cf, I)
		ggtmp	<- rbind(ggtmp, data.frame(beta = tmp1[i], alpha = tmp2[j], y = tmpy, u = tmpu))
	}
}
ggplot(ggtmp, aes(y, u)) + geom_line() + 
		facet_grid(beta ~ alpha, labeller = label_both)

#################
# Simulated MLE # 
#################
K			<- 5000									# Number of simulation draws
eps.draw	<- matrix(rgev(K*2), K, 2)				# Draws for estimation
ly.min		<- -5									# Lower bound of log(y) in simulation 

SimDist.fn	<- function(param, eps.draw, plotit = FALSE){
	beta	<- param[1]
	alpha	<- param[2]
	
	# Given the parameters, simulate optimal expenditure
	psi1.draw	<- exp(beta + eps.draw[,1])
	psi0.draw	<- exp(eps.draw[,2])
	y			<- exp.fn(psi1.draw, psi0.draw, alpha, MO, cf, I)
	
	# Derive the empirical distribution 
	p0		<- mean(y == 0)							# P(y == 0 )
	p.mo	<- mean(y == MO)						# P(y == MO)
	g		<- try(density(log(y)[y>0 & y!= MO], from = ly.min, to = log(I)), silent = T) # g(y) for continuous y
	if(class(g) == "try-error"){
		# When no values are such that y>0 & y!= MO, then we have g(y)
		g		<- function(y){0}
 	}else{
		g$y		<- g$y * (1 - p0 - p.mo)				# For continuous values, g(y) = g(y|y>0 & y!= MO)*P(y>0 & y!= MO)
		g		<- approxfun(g)							# Linear interpolation
# 		g		<- splinefun(g)							# Spline interpolation, here it does not make a difference
# 		g		<- function(v, x = y[y>0 & y!= MO], bw = bw.nrd0(y[y>0 & y!= MO])){
# 			sapply(v, function(vv)  mean(exp(-.5*((vv-x)/bw)^2)/sqrt(2*pi))/bw* (1 - p0 - p.mo) )
# 		}												# Directly write the Kernel density function
	}
	if(plotit){
		tmp	<- seq(exp(ly.min), 2*MO, length = 50)
		tmpd<- g(log(tmp))	
		par(mfrow = c(3, 1))
		print(plot(density(y), xlab="y", ylab = "density", main = "Density of simulated expenditure"))
		print(plot(density(y), xlim = range(tmp), xlab="y", ylab = "density", main = "Density of simulated expenditure"))		
		print(plot(tmp, tmpd, type = "l", xlab = "y", ylab = "density", main = "Approximated density function"))
	}
	return(list(p0 = p0, p.mo = p.mo, g = g))
}

# Check the simulated function with true value
out	<- SimDist.fn(c(beta, alpha), eps.draw)
out$p0; exp(-beta)/(1 + exp(-beta))						# The derived P(y == 0) = exp(-beta)/(1 + exp(-beta))
out$p.mo
out$p0 + out$p.mo + integrate(out$g, ly.min, log(I))$value		# Should equal to 1

ll.fn <- function(param, eps.draw){
	param[2]<- exp(param[2])
	
	### Simulation step ###
	# Given the parameters, simulate optimal expenditure y and have the distribution 
	dist0	<- SimDist.fn(param, eps.draw)
	p0		<- dist0$p0
	p.mo	<- dist0$p.mo
	g		<- dist0$g
	e		<- -100									# Small negative values for not well-define or -Inf when taking log
	psmall	<- 1e-15 
	
	### Evaluation step ###
	# Calculate the log likelihood function
	l0 		<- log(p0)
	lp.mo	<- log(p.mo)
	lp.mo	<- ifelse(lp.mo == -Inf | is.na(lp.mo), e, lp.mo)
	l.cont	<- e
	sel		<- y>0 & y!=MO
	l.cont[sel] <- log(g(log(y[sel])))
	l.cont	<- ifelse(l.cont == -Inf | is.na(l.cont), e, l.cont)
	ll		<- 	1*(y == 0)*l0 + 
				1*(y==MO)*lp.mo + 
				1*(y>0 & y!= MO)*l.cont
	return(ll)
}

# NOTE: the estimation is VERY sensitive to initial values if K is not large enough
system.time(tmp <- ll.fn(c(beta, log(alpha)), eps.draw))
theta.init	<- c(beta, log(alpha))
sol		<- maxLik(ll.fn, start = theta.init, eps.draw= eps.draw, method = "NM")
summary(sol)
matrix(cbind(c(beta, log(alpha)), coef(sol)), 2,2, dimnames = list(c("beta", "log(alpha)"), c("true", "estimate")))

theta.init	<- c(1, log(4))
sol1		<- maxLik(ll.fn, start = theta.init, eps.draw= eps.draw, method = "NM")
summary(sol1)
matrix(cbind(c(beta, log(alpha)), coef(sol1)), 2,2, dimnames = list(c("beta", "log(alpha)"), c("true", "estimate")))

# Check the shape of log likelihood function
# Fix alpha
tmp	<- seq(-1, 1, .05) + beta
tmp1	<- sapply(tmp, function(x) sum(ll.fn(c(x, log(alpha)), eps.draw)))
plot(tmp, tmp1, xlab = "beta", ylab = "Loglikelihood")
lines(tmp, tmp1)
abline(v = beta, col = "red")

# Fix beta
tmp	<- seq(-5, 2, .2) + log(alpha)
tmp1	<- sapply(tmp, function(x) sum(ll.fn(c(beta, x), eps.draw)) )
plot(tmp, tmp1, xlab = "log(alpha)", ylab = "Loglikelihood")
lines(tmp, tmp1)
abline(v = log(alpha), col = "red")



