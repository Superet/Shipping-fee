# This script explores the roles of parameters to have identification intuition.

# Model: 
# max_y psi1*log(y+1) + psi0(I-y-c_ship(y)) - c_search(y)
# s.t. y>=0 and I-y>=0
# c_ship(y) = I(y < MO)*alpha1*c_f
# c_search(y)  = I(y>=MO)*alpha2/(y-MO+1)

# psi1 = exp(x*beta + epsilon1)
# psi0 = exp(epsilon0)
# where epsilon1 and epsilon0 are iid Gumbel distribution

# This script does not include covariates in x

library(evd)
library(maxLik)
library(nloptr)
library(ggplot2)
set.seed(66)

#############
# Functions # 
#############
# Utility funciton 
u		<- function(y, psi1, psi0, alpha, MO, cf, I){
	psi1*log(y+1) + psi0*(I- y - 1*(y<MO & y>0)*alpha[1]*cf) - 1*(y>=MO)*alpha[2]/(y - MO + 1)
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

##############
# Simulation # 
##############
# Set parameters 
beta	<- 7
alpha	<- c(1, 1000)				
N		<- 3000							# Number of observations
I 		<- 250000/12					# Income		
MO		<- 1000							# Minimum order
cf		<- 100							# Shipping fee
eps		<- matrix(rgev(N*2), N, 2)		# Randon epsilon draws for psi
psi1	<- exp(beta + eps[,1])
psi0	<- exp(eps[,2])

# Baseline: alpha1 = 0, alpha2  = 0
y0 	<- exp.fn(psi1, psi0, c(0, 0), MO, cf, I)

# Scenario 1: alpha1 > 0, alpha2 = 0
y1	<- exp.fn(psi1, psi0, c(alpha[1], 0), MO, cf, I)

# Scenario 2: alpha1 = 0, alpha2 > 0
y2	<- exp.fn(psi1, psi0, c(0, alpha[2]), MO, cf, I)

# Scenario 3: alpha1 > 0, alpha2 > 0
y3	<- exp.fn(psi1, psi0, alpha, MO, cf, I)

# Scenario 4: alpha1 > 0, alpha2 > 0, large beta
delta.beta <- 2
y4	<- exp.fn(psi1 = exp(beta+delta.beta+eps[,1]), psi0, alpha, MO, cf, I)

# Scenario 5: alpha1 > 0, alpha2 > 0, large scale parameter
eps.scale <- 2
y5	<- exp.fn(psi1 = exp(beta+ eps.scale*eps[,1]), exp(eps.scale*eps[,2]), alpha, MO, cf, I)

# Plot the distribution of y 
# Baseline: continuous and smooth density				
# Scenario 1: There is a gap [gbar,MO) in distribution as long as alpha1 > 0; gbar < MO - alpha1*cf
# Scenario 2: There is a dip in distribution around MO. 
# Scenario 3: similar to scenario 1, but the spike is shorter 
# Scenario 4: distribution shifted towards right
# Scenario 5: distribution spread out to both sides. 

scn.par <- rbind(c(beta, 0, 0, 1), c(beta, alpha[1], 0, 1), c(beta, 0, alpha[2], 1), 
				 c(beta, alpha, 1), c(beta+delta.beta, alpha, 1), c(beta, alpha, eps.scale))
dimnames(scn.par)	<- list(0:5, c("beta", "alpha1", "alpha2", "scale"))				 
tmp		<- apply(scn.par, 1, function(x) paste(paste(names(x), x, sep="="), collapse = ","))
tmp	
ggtmp	<- data.frame(scenario = rep(0:5, each = N), y = c(y0, y1, y2, y3, y4, y5))
# ggtmp$scn.lab <- tmp[as.character(ggtmp$scenario)]
ggtmp$scn.lab <- factor(ggtmp$scenario, levels = 0:5, 
						labels = c("free shipping", "pure topping up", "only search cost", 
								    "topping up + search cost", "high utility", "high variance"))
ggplot(ggtmp, aes(y)) + geom_histogram() + facet_wrap(~scn.lab)
quartz()
ggplot(ggtmp, aes(y)) + geom_histogram() + facet_wrap(~scn.lab) + xlim(c(0, 2*MO))

