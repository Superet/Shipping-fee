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
library(data.table)
library(stargazer)
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
	psi1*log(y+1) + psi0*(I- y - 1*(y>0 & y<MO)*alpha[1]*cf) - 1*(y>=MO)*alpha[2]/(y - MO + 1)
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

#######################
# Identification test # 
#######################
# Look at the distribution by individual parameter
my.alpha1 	<- c(.1, .2, .3, .4, .5, 1, 3, 5)
y.ls1		<- sapply(my.alpha1, function(x) exp.fn(psi1, psi0, c(x, alpha[2]), MO, cf, I, FALSE))
colnames(y.ls1) 	<- my.alpha1
ggtmp1		<- melt(y.ls1)
ggplot(ggtmp1, aes(value)) + geom_histogram(binwidth = 1) + 
		facet_wrap( ~ Var2)

my.alpha2 	<- c(.1, .5, 1, 2, 3, 4, 5)
y.ls2		<- sapply(my.alpha2, function(x) exp.fn(psi1, psi0, c(alpha[1], x), MO, cf, I, FALSE))
colnames(y.ls2) 	<- my.alpha2
ggtmp2		<- melt(y.ls2)
quartz()
ggplot(ggtmp2, aes(value)) + geom_histogram(binwidth = 1) + 
		facet_wrap( ~ Var2)

# Expenditure distribution across a two-dimensional grid
tmp	<- c(.01, .1, .2, .3, .4, .5, 1, 3, 5)
ggtmp	<- data.frame(NULL)
for(i in 1:length(tmp)){
	for(j in 1:length(tmp)){
		tmpy	<- exp.fn(psi1, psi0, c(tmp[i], tmp[j]), MO, cf, I, FALSE)
		ggtmp	<- rbind(ggtmp, data.frame(alpha1 = tmp[i], alpha2 = tmp[j], y = tmpy))
	}
}

pdf("~/Desktop/dist_exp.pdf", width = 12, height = 10)
ggplot(ggtmp, aes(y)) + geom_histogram(binwidth = 1) + 
		facet_grid(alpha1 ~ alpha2, labeller = label_both, scales = "free_y")
dev.off()

# Probability of expenditure regions 
tmp.tbl <- data.table(ggtmp)
tmp.tbl	<- tmp.tbl[,list(p0 = mean(y==0), p1 = mean(y>0 & y<MO), 
						pMO	= mean(y==MO), p2 = mean(y>MO)), 
			by = list(alpha1, alpha2)]
tmp.tbl	<- data.frame(tmp.tbl)			

colv	<- c("p0","p1","pMO","p2")			
tmpls	<- vector("list", 4); names(tmpls) <- colv
for(i in 1:4){
	tmpls[[i]] 	<- dcast(tmp.tbl[,c("alpha1", "alpha2", colv[i])], alpha1 ~ alpha2) 
}
			
stargazer(tmpls, type = "text", summary = rep(FALSE,4), title = c("P(y=0)", "P(0<=y<MO)", "P(y=MO)", "P(y>MO)"))
stargazer(tmpls, type = "html", summary = rep(FALSE,4), title = c("P(y=0)", "P(0<y<MO)", "P(y=MO)", "P(y>MO)"), 
		out = "~/Desktop/proby.html", rownames = NULL)


