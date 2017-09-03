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
library(reshape2)
library(fastGHQuad)
library(plot3D)
set.seed(66)

# Set parameters 
beta	<- 6			
I 		<- 250000/12					
MO		<- 1000
cf		<- 100
sigma	<- 1
alpha.mu	<- c(.8, 1500)						
lalpha.mu	<- log(alpha.mu)			# alpha_h ~ lognormal(mu, sigma)
lalpha.sigma<- c(.3, 0)
nhh		<- 500							# Number of households 
nT		<- 30							# Number of time periods per household
N		<- nhh*nT						# Number of observations
hh.idx	<- rep(1:nhh, each = nT)		# Index of households
eps		<- matrix(rgev(N*2, scale = sigma), N, 2)		# Random epsilon draws for psi
alpha	<- sapply(1:2, function(i) rlnorm(nhh, lalpha.mu[i], lalpha.sigma[i]) )					
summary(alpha)
alpha	<- alpha[hh.idx,]
psi1	<- exp(beta + eps[,1])
psi0	<- exp(eps[,2])

# ---------------------------- #
# Simulate optimal expenditure #
#### Functions ####
# Utility funciton 
u		<- function(y, psi1, psi0, alpha, MO, cf, I){
	if(length(alpha) > 2){
		out <- psi1*log(y+1) + psi0*(I- y - 1*(y>0 & y<MO)*alpha[,1]*cf) - 1*(y>=MO)*alpha[,2]/(y - MO + 1)
	}else{
		out <- psi1*log(y+1) + psi0*(I- y - 1*(y>0 & y<MO)*alpha[1]*cf) - 1*(y>=MO)*alpha[2]/(y - MO + 1)
	}
	return(out)
}

# y* = psi1/psi0 - 1 from FOC and KKT conditions
y1star	<- function(psi1, psi0, I){ 
	out	<- psi1/psi0 - 1 
	return(ifelse(out > I, I, out))
}

# y** is such that psi1/(y**+1) - psi0 + alpha/(y** - MO+ 1)^2 = 0 for y** >= MO
y2star	<- function(psi1, psi0, alpha, MO, I){
	foc <- function(y){psi1/(y+1) - psi0*alpha[1] + alpha[2]/(y - MO+ 1)^2 }
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
		y2[sel[i]]	<- y2star(psi1[sel[i]], psi0[sel[i]], alpha[sel[i],], MO, I)
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
hist(y, breaks = 100)
hist(y, xlim = c(0, 2*MO), breaks = 300)
table(y == 0)
table(y == MO)
table(y > MO)

#################
# Simulated MLE # 
#################
K			<- 5000									# Number of simulation draws
eps.draw	<- matrix(rgev(K*2), K, 2)				# Draws for estimation
num.nodes	<- 8									# Number of quadrature nodes
ly.min		<- -5									# Lower bound of log(y) in simulation 

SimDist.fn	<- function(beta, alpha, sigma, esp.draw, mod.log = TRUE){	
	# Given the parameters, simulate optimal expenditure
	psi1.draw	<- exp(beta + eps.draw[,1]*sigma)
	psi0.draw	<- exp(eps.draw[,2]*sigma)
	y			<- exp.fn(psi1.draw, psi0.draw, alpha, MO, cf, I)
	
	# Derive the empirical distribution 
	p0		<- mean(y == 0)							# P(y == 0 )
	p.mo	<- mean(y == MO)						# P(y == MO)
	if(mod.log){
		g		<- try(density(log(y[y>0 & y!= MO]), from = ly.min, to = log(I)), silent = T) # g(y) for continuous y
	}else{
		g		<- try(density(y[y>0 & y!= MO], from = 0, to = I), silent = T) # g(y) for continuous y
	}
	if(class(g) == "try-error"){
		# When no values are such that y>0 & y!= MO, then we have g(y)
		g		<- function(y){0}
 	}else{
		g$y		<- g$y * (1 - p0 - p.mo)				# For continuous values, g(y) = g(y|y>0 & y!= MO)*P(y>0 & y!= MO)
		g		<- approxfun(g)							# Linear interpolation											# Directly write the Kernel density function
	}
	return(list(p0 = p0, p.mo = p.mo, g = g, y = y))
}

# Check the simulated function with true value
out	<- SimDist.fn(beta, alpha, sigma, eps.draw)
out$p0; exp(-beta)/(1 + exp(-beta))						# The derived P(y == 0) = exp(-beta)/(1 + exp(-beta))
out$p.mo
out$p0 + out$p.mo + integrate(out$g, ly.min, log(I))$value		# Should equal to 1
hist(out$y, main = "Simulated y")

ll.fn <- function(param, eps.draw, mod.log = TRUE){
	beta	<- param[1]
	lalpha.mu	<- param[2:3]
	lalpha.sigma<- exp(param[4:5])
	sigma	<- exp(param[6])
	
	# Replace y with log(y)
	if(mod.log){
		sel <- y>0 & y!=MO
		y[sel]	<- log(y[sel])
	}

	# Quadrature nodes 
	q.grid	<- gaussHermiteData(num.nodes)$x
	sel 	<- which(lalpha.sigma < 1e-10)
	if(length(sel) == 0){
		tmp1	<- sqrt(2)*lalpha.sigma[1]*q.grid  + lalpha.mu[1]
		tmp2	<- sqrt(2)*lalpha.sigma[2]*q.grid  + lalpha.mu[2]
		q.wt	<- gaussHermiteData(num.nodes)$w
		q.grid	<- expand.grid(tmp1, tmp2)															
		q.wt	<- kronecker(q.wt, q.wt) 	# Tensor product of grids in two dimensions		
	}else{
		q.grid	<- sqrt(2)*lalpha.sigma[-sel]*q.grid  + lalpha.mu[-sel]
		q.grid	<- cbind(q.grid, lalpha.mu[sel])
		q.wt	<- gaussHermiteData(num.nodes)$w
	}
		
	# Simulate expenditure at each nodes
	p0		<- p.mo	<- rep(NA, nrow(q.grid))
	g		<- vector("list", length = nrow(q.grid))
	es		<- 1e-8
	for(i in 1:nrow(q.grid)){
		dist0	<- SimDist.fn(beta, rep(1, K) %*% matrix(exp(q.grid[i,]), nrow = 1), sigma, eps.draw)
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
		phh.grid	<- matrix(NA, sum(sel), nrow(q.grid))
		for(j in 1:nrow(q.grid)){
			p.cont	<- g[[j]](y[sel])
			p.cont	<- ifelse(p.cont == 0 | is.na(p.cont), es, p.cont)
			phh.grid[,j]	<- 	1*(y[sel] == 0)*p0[j] + 
								1*(y[sel] == MO)*p.mo[j] + 
								1*(y[sel]>0 & y[sel]!= MO)*p.cont
		}
		phh.grid[phh.grid < es]	<- es
		Ephh		<- 1/sqrt(pi)*sum(q.wt * apply(phh.grid, 2, prod))
		ll.hh[i]	<- log(Ephh)
	}
	
	# Integrate the likelihood over the quadrature grids
	ll	<- sum(ll.hh)
	
	return(ll)
}

# -------------------------------------- #
# Check the shape of likelihood function # 
# We focus on alpha1.mu and alpha1.sigma here 
tmp1	<- seq(-.6, 1, by = .1)
tmp2	<- seq(-2.5, -.5, by = .5)	#log(seq(.1, .8, by = .1))
ggtmp	<- data.frame(NULL)
for(i in 1:length(tmp1)){
	for(j in 1:length(tmp2)){
		tmpl	<- ll.fn(c(beta, tmp1[i], lalpha.mu[2], tmp2[j], -100, log(sigma)), eps.draw)
		ggtmp	<- rbind(ggtmp, data.frame(lalpha1.mu = tmp1[i], lalpha1.sigma = tmp2[j], ll = sum(tmpl)))
	}
}

# Marginal log likelihood function
# the optimal alpha1 moves together alpha2
ggplot(ggtmp, aes(lalpha1.mu, ll)) + geom_line() + facet_wrap(~lalpha1.sigma, labeller = "label_both")
quartz()
ggplot(ggtmp, aes(lalpha1.sigma, ll)) + geom_line() + facet_wrap(~lalpha1.mu, labeller = "label_both")

# Plots of joint log likelihood function 
z <- dcast(ggtmp, lalpha1.mu ~ lalpha1.sigma)
z <- as.matrix(z[,-1])
persp3D(tmp1, tmp2, z = z, theta=150, phi=20, r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, 
  nticks=5, ticktype="detailed", col="cyan", xlab="lalpha1.mu", 
  ylab="lalpha1.sigma", zlab="ll")

# ---------------#
# Run estimation # 
theta.init	<- c(beta, lalpha.mu, log(lalpha.sigma[1]), -100, log(sigma))				
names(theta.init) <- c("beta", "lalpha1.mu", "lalpha2.mu", "lalpha1.sigma", "lalpha2.sigma", "log(sigma)")
system.time(tmp <- ll.fn(theta.init, eps.draw))
pct		<- proc.time()

# Fix scale parameter and lalpha1.sigma, lalpha2.mu, lalpha2.sigma
sol1	<- maxLik(ll.fn, start = theta.init, eps.draw= eps.draw, method = "NM", fix = c(3,4,5,6))
summary(sol1)

# Fix scale parameter and lalpha2.mu, lalpha2.sigma
sol2	<- maxLik(ll.fn, start = theta.init, eps.draw= eps.draw, method = "NM", fix = c(4,5,6))
summary(sol2)

# Fix scale parameter and lalpha2.sigma
sol3	<- maxLik(ll.fn, start = theta.init, eps.draw= eps.draw, method = "NM", fix = c(5,6))
summary(sol3)

# Fix lalpha2.sigma
sol4	<- maxLik(ll.fn, start = theta.init, eps.draw= eps.draw, method = "NM", fix = 5)
summary(sol4)

# Fix scale parameter
sol5	<- maxLik(ll.fn, start = theta.init, eps.draw= eps.draw, method = "NM", fix = 6)
summary(sol5)

# Have a reasonable guess of parameters 
theta.init1	<- c(log(2000), log(.1), log(1200), log(.5), -100, 0)
names(theta.init1) <- c("beta", "lalpha1.mu", "lalpha2.mu", "lalpha1.sigma", "lalpha2.sigma", "log(sigma)")
sol6	<- maxLik(ll.fn, start = theta.init1, eps.draw= eps.draw, method = "NM", fix = 5)
summary(sol6)
cat("All estimation finishes with", (proc.time() - pct)[3]/60, "min.\n")

# Summarize estimation 
# See which gets the highest log likelihood
sol.nm 	<- c("true", "fix a2, a.sig, sigma", "fix a2, a2.sig, sigma",  
			 "fix a2.sig, sigma", "fix a2.sig", "fix sigma", "good init")
tmp1	<- c(sum(tmp), sol1$maximum, sol2$maximum, sol3$maximum, sol4$maximum, sol5$maximum, sol6$maximum)
names(tmp1) <- sol.nm
round(tmp1, 3)
round(tmp1[order(tmp1, decreasing = T)], 3)

tmp		<- cbind(theta.init, coef(sol1), coef(sol2),  coef(sol3), coef(sol4),  coef(sol5), coef(sol6))
colnames(tmp) 	<- sol.nm
round(tmp, 3)


