# Model: 
# max_y psi1*log(y+1) + psi0(I-y-c_ship(y)) - c_search(y)
# s.t. y>=0 and I-y>=0
# c_ship(y) = I(y < MO)*c_f
# c_search(y)  = I(y>=MO)*alpha/(y-MO+1)

# psi1 = exp(x*beta + epsilon1)
# psi0 = exp(epsilon0)
# where epsilon1 and epsilon0 are iid Gumbel distribution

# This script does not include covariates in x

### NOTE ###
# 1. Both alpha1 and alpha2 can be identified, if enough observations around the region (y - epsilon, y + epsilon); 
# 2. When beta is small, yielding few observations that y > MO, alpha2 is not identified well. 
# 3. The number of simulation draws needs to be large enougth, otherwise the likelihood function is not smooth. 


library(evd)
library(maxLik)
library(nloptr)
library(ggplot2)
library(reshape2)
library(plot3D)

set.seed(66)

# Set parameters 
beta	<- 2
alpha	<- c(.5, 1)						# Coefficients for shipping cost, and search cost
nhh		<- 100							# Number of households 
nT		<- 50							# Number of time periods per household
N		<- nhh*nT						# Number of observations
I 		<- 100					
MO		<- 20
cf		<- 8
eps		<- matrix(rgev(N*2), N, 2)		# Randon epsilon draws for psi
psi1	<- exp(beta + eps[,1])
psi0	<- exp(eps[,2])

# ---------------------------- #
# Simulate optimal expenditure #
#### Functions ####
# Utility funciton 
u		<- function(y, psi1, psi0, alpha, MO, cf, I){
	psi1*log(y+1) + psi0*(I- y - 1*(y<MO)*alpha[1]*cf) - 1*(y>=MO)*alpha[2]/(y - MO + 1)
}

# y* = psi1/psi0 - 1 from FOC and KKT conditions
y1star	<- function(psi1, psi0, I){ 
	out	<- psi1/psi0 - 1 
	return(ifelse(out > I, I, out))
}

# y** is such that psi1/(y**+1) - psi0 + alpha/(y** - MO+ 1)^2 = 0 for y** >= MO
y2star	<- function(psi1, psi0, alpha, MO, I){
	foc <- function(y){psi1/(y+1) - psi0 + alpha[2]/(y - MO+ 1)^2 }
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

SimDist.fn	<- function(param, esp.draw){
	beta	<- param[1]
	alpha	<- param[c(2,3)]
	
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
		g		<- approxfun(g)							# Linear interpolation
		# g		<- splinefun(g)							# Spline interpolation, here it does not make a difference
		# g		<- function(v, x = y[y>0 & y!= MO], bw = bw.nrd0(y[y>0 & y!= MO])){
		# 	sapply(v, function(vv)  mean(exp(-.5*((vv-x)/bw)^2)/sqrt(2*pi))/bw* (1 - p0 - p.mo) )
		# }												# Directly write the Kernel density function
	}
	return(list(p0 = p0, p.mo = p.mo, g = g))
}

# Check the simulated function with true value
out	<- SimDist.fn(c(beta, alpha), eps.draw)
out$p0; exp(-beta)/(1 + exp(-beta))						# The derived P(y == 0) = exp(-beta)/(1 + exp(-beta))
out$p.mo
out$p0 + out$p.mo + integrate(out$g, 0, I)$value		# Should equal to 1

ll.fn <- function(param, eps.draw){
	param[2]<- exp(param[2])
	param[3]<- exp(param[3])
	
	### Simulation step ###
	# Given the parameters, simulate optimal expenditure y and have the distribution 
	dist0	<- SimDist.fn(param, eps.draw)
	p0		<- dist0$p0
	p.mo	<- dist0$p.mo
	g		<- dist0$g
	e		<- -100									# Small negative values for not well-define or -Inf when taking log
	
	### Evaluation step ###
	# Calculate the log likelihood function
	l0 		<- log(p0)
	lp.mo	<- log(p.mo)
	lp.mo	<- ifelse(lp.mo == -Inf | is.na(lp.mo), e, lp.mo)
	l.cont	<- log(g(y))
	l.cont	<- ifelse(l.cont == -Inf | is.na(l.cont), e, l.cont)
	ll		<- 	1*(y == 0)*l0 + 
				1*(y==MO)*lp.mo + 
				1*(y>0 & y!= MO)*l.cont
	return(ll)
}


sum(ll.fn(c(beta, log(alpha)), eps.draw))
sum(ll.fn(c(beta, log(c(.6, 1.79))), eps.draw))
sum(ll.fn(c(beta, log(c(.7, 1.8))), eps.draw))

# Check the likelihood function
tmp	<- seq(.1, 2.5, by = .1)
ggtmp	<- data.frame(NULL)
for(i in 1:length(tmp)){
	for(j in 1:length(tmp)){
		tmpl	<- ll.fn(c(beta, log(c(tmp[i], tmp[j]))), eps.draw)
		ggtmp	<- rbind(ggtmp, data.frame(alpha1 = tmp[i], alpha2 = tmp[j], y = sum(tmpl)))
	}
}

ggplot(ggtmp, aes(alpha1, y)) + geom_line() + facet_wrap(~alpha2)
ggplot(ggtmp, aes(alpha1, alpha2, z = y)) + geom_contour(aes(colour = ..level..))

z <- dcast(ggtmp, alpha1 ~ alpha2)
z <- as.matrix(z[,-1])
persp3D(tmp, tmp, z = z, theta=150, phi=20, r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, 
  nticks=5, ticktype="detailed", col="cyan", xlab="alpha1", 
  ylab="alpha2", zlab="ll")


scatter3d(ggtmp$alpha1, y = ggtmp$alpha2, z = ggtmp$y, surface = T)



# NOTE: the estimation is VERY sensitive to initial values if K is not large enough
system.time(tmp <- ll.fn(c(beta, log(alpha)), eps.draw))
beta.init	<- c(beta, log(alpha))
sol		<- maxLik(ll.fn, start = beta.init, eps.draw= eps.draw, method = "NM")
summary(sol)
matrix(cbind(c(beta, log(alpha)), coef(sol)), 3,2, dimnames = list(c("beta", "log(alpha1)", "log(alpha2)"), c("true", "estimate")))

beta.init	<- c(-.8, log(.9), log(2))
sol1		<- maxLik(ll.fn, start = beta.init, eps.draw= eps.draw, method = "NM")
summary(sol1)
matrix(cbind(c(beta, log(alpha)), coef(sol1)), 3,2, dimnames = list(c("beta", "log(alpha1)", "log(alpha2)"), c("true", "estimate")))

# Check the shape of log likelihood function
# Fix beta and alpha[2]
tmp	<- seq(-1, 1, .05)
ggtmp	<- rbind(data.frame(vary.var = "alpha1", alpha1 = log(alpha[1]) + tmp, alpha2 = log(alpha[2])), 
				 data.frame(vary.var = "alpha2", alpha1 = log(alpha[1]), alpha2 = log(alpha[2]) + tmp))				
ggtmp$ll<- sapply(1:nrow(ggtmp), function(i) sum(ll.fn(c(beta, ggtmp[i,"alpha1"], ggtmp[i,"alpha2"]), eps.draw)) )
ggtmp$value	<- ifelse(ggtmp$vary.var == "alpha1", ggtmp$alpha1, ggtmp$alpha2)
ggtmp$true 	<- ifelse(ggtmp$vary.var == "alpha1", log(alpha[1]), log(alpha[2]))

ggplot(ggtmp, aes(value, ll)) + geom_line() + geom_point() +
 		geom_vline(aes(xintercept = true), col = "red", size = .25, linetype = 2) + 
		facet_wrap(~ vary.var, scales = "free_x")


# Check the shape of log likelihood function
# Fix alpha
tmp	<- seq(-1, 1, .05) + beta
tmp1	<- sapply(tmp, function(x) sum(ll.fn(c(x, log(alpha)), eps.draw)))
plot(tmp, tmp1, xlab = "beta", ylab = "Loglikelihood")
lines(tmp, tmp1)
abline(v = beta, col = "red")

# Fix beta
tmp	<- seq(-1, 1, .05) + log(alpha)
tmp1	<- sapply(tmp, function(x) sum(ll.fn(c(beta, x), eps.draw)))
plot(tmp, tmp1, xlab = "log(alpha)", ylab = "Loglikelihood")
lines(tmp, tmp1)
abline(v = log(alpha), col = "red")



