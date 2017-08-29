# Model: 
# max_y psi1*log(y+1) + psi0(I-y-c_ship(y)) - c_search(y)
# s.t. y>=0 and I-y>=0
# c_ship(y) = I(y < MO)*c_f
# c_search(y)  = I(y>=MO)*alpha/(y-MO+1)

# psi1 = exp(x*beta + epsilon1)
# psi0 = exp(epsilon0)
# where epsilon1 and epsilon0 are iid Gumbel distribution

# Estimation is within-loop simulated MLE. 

library(evd)
library(maxLik)
library(nloptr)
library(reshape2)
library(mgcv)
library(mixtools)
# setwd("~/Documents/Research/ShippingFees/Processed_data")
setwd("U:/Users/ccv103/Documents/Research/ShippingFee")
load("init_est.rdata")

# A random sub-sample
# set.seed(111)
# idx	<- sample(1:nrow(init_est), 5000)
# init_est	<- init_est[idx,]

###################
#### Functions ####
###################
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
	u3.y1	<- ifelse(is.na(u3.y1), 0, u3.y1)
	u3.MO	<- ifelse(is.na(u3.MO), 0, u3.MO)
	u3.y2	<- ifelse(is.na(u3.y2), 0, u3.y2)

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

# Log likelihood function # 
ll.fn <- function(param){
	beta	<- param[1:nx]
	alpha	<- exp(param[length(param)])			# Make sure alpha > 0 

	### Simulation step ###
	# Given the parameters, first simulate expenditure for each unique X and each draw
	d		<- X.unq %*% beta
	ysim	<- matrix(NA, K, nX.unq)
	for(i in 1:nX.unq){
		# Given the parameters d and alpha, simulate optimal expenditure
		psi1.draw	<- exp(d[i] + eps.draw[,1])
		psi0.draw	<- exp(eps.draw[,2])
		ysim[,i]	<- exp.fn(psi1.draw, psi0.draw, alpha, MO, cf, I)
	}
	
	# Calculate probability and pdf of expenditure under the current parameter
	p0.unq	<- colMeans(ysim == 0)
	pmo.unq	<- colMeans(ysim == MO)
	g.ls	<- vector("list", length = nX.unq)		# A density function for each X
	for(i in 1:nX.unq){
		g	<- try(density(ysim[,i][ysim[,i]>0 & ysim[,i]!=MO], from = 0, to = I), silent = T) # g(y) for continuous y
		if(class(g) == "try-error"){
			# When no values are such that y>0 & y!= MO, then we have g(y)
			g		<- function(y){0}
	 	}else{
			# Density function conditional on y>0
			g$y		<- g$y * (1 - p0.unq[i] - pmo.unq[i])/(1 - p0.unq[i])		# For continuous values, g(y) = g(y|y>0 & y!= MO)*P(y>0 & y!= MO)
			g		<- approxfun(g)							# Linear interpolation									
		}
		g.ls[[i]]	<- g
	}
	pmo.unq	<- pmo.unq/(1 - p0.unq)
	
	### Evaluation step ###
	# Calculate the log likelihood function
	eps.e	<- -100
	lp.mo	<- log(pmo.unq[x.idx])
	lp.mo	<- ifelse(lp.mo == -Inf | is.na(lp.mo), eps.e, lp.mo)
	l.cont	<- rep(eps.e, nrow(X))
	for(i in 1:nX.unq){
		sel	<- x.idx == i
		l.cont[sel]	<- log(g.ls[[i]](y[sel]))
	}
	l.cont	<- ifelse(l.cont == -Inf | is.na(l.cont), eps.e, l.cont)
	ll		<- 	1*(y==MO)*lp.mo + 
				1*(y>0 & y!= MO)*l.cont
	return(ll)
}

#####################
### Organize data ###
#####################
head(init_est)
dim(init_est)

N		<- nrow(init_est)						# Number of observations
I 		<- 250000/12							# Monthly income
MO		<- 1000									# Minimum order of free shipping
cf		<- 100									# Shipping fee
X		<- as.matrix(cbind(1, log(init_est[,c("items_count", "mean_price")])))
summary(X)
nx		<- ncol(X)
y		<- init_est[,"preshipping"]
table(y==0)
sum(y>=I)							
hist(y, breaks = 50); abline(v = MO, col = "red")

# Simple correlation 
cor(cbind(y, X[,-1]))
summary(lm(y ~ X-1))

# Unique X 
X.unq	<- unique(X)
x.idx	<- rep(NA, nrow(X))
nX.unq	<- nrow(X.unq)
for(i in 1:nX.unq){
	sel <- rowSums(X - rep(1, nrow(X)) %*% matrix(X.unq[i,], nrow = 1)) == 0 
	if(sum(sel) > 0){
		x.idx[sel]	<- i
	}
}

#####################
### Simulated MLE ###
#####################
K			<- 1000									# Number of simulation draws, need to be large, 1000 + gives better estimates
eps.draw	<- matrix(rgev(K*2), K, 2)				# Draws for estimation

# ------------------------#
# Run MLE optimization
par.init	<- c(Intercept = -10, iterms_count = 0, mean_price= 2, log_alpha=log(2))			# beta, log(alpha)
system.time(tmp <- ll.fn(par.init))
# par.initls	<- list(par.init, 
# 					c(Intercept = -1, iterms_count = .6, mean_price= 0.1, log_alpha=log(2)),
# 					c(Intercept = 3, iterms_count = -1, mean_price= -2, log_alpha=log(10)),
# 					c(Intercept = 1, iterms_count = .1, mean_price= -.2, log_alpha=log(3)), 
# 					c(Intercept = 0, iterms_count = -.1, mean_price= -.1, log_alpha=log(2))
# 					)
par.initls	<- list(par.init, 
					c(Intercept = -1, iterms_count = .6, mean_price= 0.1, log_alpha=log(20)),
					c(Intercept = -3, iterms_count = -.01, mean_price= .5, log_alpha=log(10)),
					c(Intercept = -1, iterms_count = .1, mean_price= .2, log_alpha=log(50)), 
					c(Intercept = 0, iterms_count = -.05, mean_price= .1, log_alpha=log(60))
					)
pct		<- proc.time()
sol.ls	<- lapply(par.initls, function(x) maxLik(ll.fn, start = x, method = "NM", control = list(iterlim = 1000)) )
(sel		<- which.max(sapply(sol.ls, function(x) x$maximum)))
sol		<- sol.ls[[sel]]
summary(sol)
par.est	<- coef(sol)
use.time<- proc.time() - pct
cat("Estimation uses", use.time[3]/60,"min.\n")

# Check the shape of log likelihood function
idx	<- 1
tmp	<- seq(-1, 1, .2) 
tmp1	<- sapply(tmp, function(x){p1 <- par.est; p1[idx] <- p1[idx]+x; sum(ll.fn(p1))} )
plot(tmp + par.est[idx], tmp1, xlab = paste("par[",idx,"]", sep=""), ylab = "Loglikelihood")
lines(tmp + par.est[idx], tmp1)
abline(v = par.est[idx], col = "red")

save.image(file = paste("init_est_inloop_", Sys.Date(), ".rdata",sep=""))
