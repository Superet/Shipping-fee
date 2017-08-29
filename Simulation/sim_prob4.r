# Model: 
# max_{b1, b2} psi1*log(b1+1) + psi2*log(b2+1) + psi0*log(I - y- c_ship(y)) - c_search(y)
# s.t. b1 >= 0, b2 >= 0 and I- y - c_ship(y)>=0, b1 + b2 = y
# c_ship(y) = I(y < MO)*c_f
# c_search(y)  = I(y>=MO)*alpha/(y-MO+1)

# psi_i= exp(x_i*beta + epsilon_i), for i = 1, 2
# psi0 = exp(epsilon0)
# where epsilon's are iid Gumbel distribution

library(evd)
library(maxLik)
library(nloptr)
library(reshape2)
library(mgcv)
library(mixtools)
library(AlgDesign)
set.seed(66)

# Set parameters 
C		<- 2							# Number of categories
beta	<- c(-2.5, -2, 1.5)				# Model parameters, psi_i = X_i*beta, the first two are intercepts for each category, the third one is for covariate
alpha	<- 1							# Search parameter
N		<- 2000							# Number of observations
I 		<- 100							# Income
MO		<- 20							# Minimum order
cf		<- 8							# Shipping fee
X		<- matrix(runif(N*C, -1, 1), N, C)	# A single covariate
X_list	<- lapply(1:C, function(i){out<- matrix(0, N, C); out[,i] <- 1; cbind(out, X[,i])})		
eps		<- matrix(rgev(N*(C+1)), N, C+1)		# Randon epsilon draws for psi
psi		<- sapply(1:C, function(i) exp(X_list[[i]] %*% beta + eps[,i]))
psi0	<- exp(eps[,(C+1)])
nx		<- length(beta)

# ---------------------------- #
# Simulate optimal expenditure #

#### Functions ####
# Utility funciton 
u		<- function(e_vec, psi_vec, psi0, alpha, MO, cf, I){
# psi1*log(b1+1) + psi2*log(b2+1) + psi0*log(I - y- c_ship(y)) - c_search(y)	
	y	<- sum(e_vec)
	sum(psi_vec*log(e_vec+1)) + psi0*(I- y - 1*(y<MO & y>0)*cf) - 1*(y>=MO)*alpha/(y - MO + 1)
}

# Function of solving the optimal expenditures
exp.fn	<- function(psi, psi0,alpha, MO, cf, I, opt.method = "nlopt"){
# Two optimization methods: constrOptim and nloptr
# We need multiple initial values for constrOpitm, 
# nloptr tend to give better optimization results. 	
	N		<- nrow(psi)
	e_mat	<- matrix(0, N, C)
	
	if(opt.method == "constr"){
		# Constraints: ui %*% theta - ci >= 0.
		ui		<- matrix(c(1, 0, 0, 1, -1, -1), 3, C, byrow = T)
		ci		<- c(0, 0, -I)
		e_init	<- c(MO, MO)*.25				# Initial values

		for(i in 1:N){
			obj.fn	<- function(e_vec){ - u(e_vec, psi[i,], psi0[i],alpha, MO, cf, I)}
			sol		<- try(constrOptim(e_init, f = obj.fn, grad = NULL, ui = ui, ci = ci), silent = T)
			if(class(sol) == "try-error" | sol$convergence != 0){
				e_mat[i,]	<- NA
			}else{
				e_mat[i,]	<- sol$par
			}
		}
	}else if(opt.method == "nlopt"){
		# Set up nloptr controls 
		# Possible non-gradient-based algorithems: NLOPT_LN_COBYLA, NLOPT_GN_DIRECT_L. 
		# 200 iterations usually give very close results, although the algorithm gives a warning that says it might not converge yet. 
		nloptr.opts	<- list("algorithm"="NLOPT_LN_COBYLA", "xtol_abs" = 1e-2, "ftol_abs" = 1e-4, "maxeval" = 200)
		
		# Inequality constraints: b1 + b2 - I <= 0
		eval_g 		<- function(e_vec){ e_vec[1] + e_vec[2] - I }					 
		
		for(i in 1:N){
			# When solving for continuous b1, b2, it is almost impossible to get exact y = MO, 
			# Therefore, we also find the maximum utility when y = MO, i.e., solving for b1, MO - b1, 
			# 			then find the solution with larger utility value. 
			
			# Objective function u(b1, b2)
			obj.fn	<- function(e_vec){ - u(e_vec, psi[i,], psi0[i],alpha, MO, cf, I)}
			sol		<- try(nloptr(x0= c(0.5*MO, MO), eval_f = obj.fn, lb = rep(0, C), ub = rep(I, C), eval_g_ineq = eval_g,
								opts = nloptr.opts), silent=TRUE)
											
			# Objective function u(b1, MO-b1)
			obj.fn1	<- function(b1)	{- u(c(b1, MO-b1), psi[i,], psi0[i],alpha, MO, cf, I)}
			sol1	<- try(optimize(f = obj.fn1, lower = 0, upper = MO), silent=TRUE)			# Univariate optimaization -- Brent algorithem
			
			if(sol$objective < sol1$objective){	
				e_mat[i,]	<- ifelse(sol$solution < 1e-3, 0, sol$solution)
			}else{		# When u(b1, b2) < u(b1, MO - b1)
				b1			<- ifelse(sol1$minimum < 1e-3, 0, sol1$minimum)
				e_mat[i,]	<- c(b1, MO - b1)
			}				
		}
	}

	return(e_mat)
}

opt.method 	<- "nlopt"
e_mat 	<- exp.fn(psi, psi0, alpha, MO, cf, I, opt.method)
summary(e_mat)
y		<- rowSums(e_mat)
summary(y)
hist(y)
table(y == 0)
table(y == MO)
table(y == I)

#################
# Simulated MLE # 
#################
# Steps: 
# 1. Set up simulation -- have K set of epsilon draws and a set of pararameter grids (d1, d2, alpha)
# 	1a. For large C, (d1, ..., dC, alpha) has too many combinations, and thus simulation takes long time. 
#		We use fractional factorial design to reduce the number of evaluation points, i.e., ngrid.  
# 2. Solve for optimal expenditure for each set of parameters and each set of random draws. 
# 3. Model conditional distribution of y	
# 	3a: P(y=0|d1, ..., dC, alpha) and P(y=MO|d1, ..., dC, alpha) are modeled with logit regressions. 
# 	3b: g(y|d1, ..., dC, alpha) is approximated with finite mixture model 

### NOTE: #### 
# The current parameters are set for a (relatively) quick run but not give us good estiamtes, 
# It's necessary to have a large number of draws and fine parameter grids to have accurate estiamtes. 
# E.g., the simulation parameters: 
# K			... 300+, or even 1000
# The number of the grids for d.seq and alpha.seq is 50+, or even larger
# ngrid		... 2500+, or even more 

K			<- 100									# Number of simulation draws, need to be large, 1000 + gives better estimates
eps.draw	<- matrix(rgev(K*(C+1)), K, C+1)		# Draws for estimation
summary(sapply(1:C, function(i) X_list[[i]] %*% beta ))

# ------------------------------------------------------#
# Let d_i = X_i*beta
# Simulate the distribution of (y|d1, d2, alpha)
# This step can be slow, but we only need to simulate once 
d.seq 		<- seq(-4, 2, length = 40)					# Need to have "dense" points around true values. 
alpha.seq 	<- exp(seq(log(.1), log(2.5), length = 30))
par.grid	<- as.matrix(expand.grid(d.seq, d.seq, alpha.seq))
cat("Before reducing the dimension of evaluation grid, dim(par.grid) =", dim(par.grid), "\n")

# Notice that par.grid has a large dimension, so to reduce computation time, 
# We use fractional factorial design to reduce the number of evaluation points.
# Below we use find exact design using Federov's exchange algorithm.
pct			<- proc.time()
ngrid		<- 1000									# The number of evaluation grids
# Run exact design. It's faster when setting approximate = TRUE
des			<-optFederov(~quad(.),par.grid,nTrials=ngrid,eval=TRUE, approximate=FALSE)			
par.grid	<- as.matrix(des$design[,c("Var1", "Var2", "Var3")])
cat("After reducing the dimension of evaluation grid, dim(par.grid) =", dim(par.grid), "\n")
colnames(par.grid)	<- c("d1", "d2", "alpha")
use.time	<- proc.time() - pct
cat("Fractional factorial design uses", use.time[3]/60, "min.\n")

# Run simulation -- solving for opitmal expenditure under each set of parameters and each set of draws 
esim		<- array(NA, c(K, nrow(par.grid), C))
ysim		<- matrix(NA, K, nrow(par.grid))
pct			<- proc.time()
for(i in 1:nrow(par.grid)){
	# Given the parameters d and alpha, simulate optimal expenditure
	psi.draw	<- exp(rep(1, K) %*% t(par.grid[i,1:C])+ eps.draw[,1:C])
	psi0.draw	<- exp(eps.draw[,(C+1)])
	esim[,i,]	<- exp.fn(psi.draw, psi0.draw, par.grid[i,(C+1)], MO, cf, I, opt.method)
	ysim[,i]	<- rowSums(esim[,i,])
}
summary(c(esim))
use.time	<- proc.time() - pct
cat("Simulation uses", use.time[3]/60, "min.\n")

# --------------------------------------- #
# Model the conditional distribution of y #
# Estimating conditional distribution is costly. 
# But we only need to estimate it once outside the MLE optimization. 

ysim.df		<- melt(ysim)
colnames(ysim.df)	<-	c("draw.idx", "par.idx", "value")
ysim.df		<- cbind(ysim.df, par.grid[ysim.df$par.idx,])
ysim.df$y0	<- 1*(ysim.df$value == 0)
ysim.df$y.mo<- 1*(ysim.df$value == MO)
mean(ysim.df$y0 == 1)
mean(ysim.df$y.mo == 1)

# Conditional probabity of y == 0
# We model logit(P(y==0)) = cbind(1, d, alpha) %*% Delta1
# We can let the logistic model to be polynomial of (d, alpha). 
p0.mod		<- glm(y0 ~ d1 + d2 + alpha, data = ysim.df, family = binomial(link = "logit"))
summary(p0.mod)
Delta1		<- coef(p0.mod)
p0.fn		<- function(xnew){
	xbeta	<- cbind(1, xnew) %*% Delta1
	out		<- exp(xbeta)/(1 + exp(xbeta))
	return(out)
} 

# Conditional probability of y == MO
# We model logit(P(y==MO)) = cbind(1, d, alpha) %*% Delta2
pmo.mod		<- glm(y.mo ~ d1 + d2 + alpha, data = ysim.df, family = binomial(link = "logit"))
summary(pmo.mod)
Delta2		<- coef(pmo.mod)
pmo.fn		<- function(xnew){
	xbeta	<- cbind(1, xnew) %*% Delta2
	out		<- exp(xbeta)/(1 + exp(xbeta))
	return(out)
}

# Conditional density of (y|d, alpha) for y > 0 and y!= MO
# Normal mixture: p(y|X) = \sum_j lambda_j * Normal(y; X*beta_j, sigma_j)
# We set the number of mixture to be 2 for illustration, 3 mixture components produce better estimates but EM converges slowly. 
num.comp<- 2					# Number of mixture components
sel		<- ysim.df$value > 0 & ysim.df$value != MO
ysimp.df<- ysim.df[sel,]
plot(density(ysimp.df$value))
x1		<- as.matrix(ysimp.df[,c("d1", "d2","alpha")])

# Take initial values for the mixture model 	
# Takes possible modes of the distribution of y as cutoff points
cut.pt	<- c(0, MO, Inf)					# For num.comp = 2				
# cut.pt	<- c(0, MO, I - cf, Inf)				# For num.comp = 3
tmp		<- as.numeric(cut(ysimp.df$value, cut.pt, labels = 1:num.comp))			# Cut y into a few bins 
tmp1	<- prop.table(table(tmp))												# The probability of y falling into each bin 
tmp2	<- lapply(1:num.comp, function(i) lm(ysimp.df[tmp==i, "value"] ~ x1[tmp==i,]) )			# Fit linear regression of y ~ d,alpha for each bin
b.init	<- sapply(tmp2, coef)		
s.init	<- sapply(tmp2, function(x) summary(x)$sigma)

pct		<- proc.time()
mix.out <- regmixEM(y= ysimp.df$value, x= x1, verb = TRUE, epsilon = 1e-03, k = num.comp, 
					lambda = as.vector(tmp1), beta = b.init, sigma = s.init) 	
summary(mix.out)
use.time	<- proc.time() - pct
cat("Mixutre model uses", use.time[3]/60, "min.\n")

# This gives the conditional density function (y|d, alpha)
fhat.fn	<- function(ynew, xnew, lambda, beta, sigma){
	mu	<- apply(beta, 2, function(b) cbind(1,xnew) %*% b)				# Location_j = X %*% beta_j
	out	<- rowSums(sapply(1:num.comp, function(i) lambda[i]*dnorm(ynew, mu[,i], sigma[i])) )
	# Scale the conditional density so that the sum of probability across y equals to 1
	out		<- out * (1 - p0.fn(xnew) - pmo.fn(xnew) )	
	return(out)
}

# ------------------------#
# Log likelihood function # 
# ------------------------#
ll.fn <- function(param){
	beta	<- param[1:nx]
	alpha	<- exp(param[length(param)])			# Make sure alpha > 0 

	### Simulation step ###
	# Given the parameters, calculate the theoretical distribution 
	d		<- sapply(X_list, function(x) x%*% beta)
	p0		<- p0.fn(cbind(d, alpha))
	p.mo	<- pmo.fn(cbind(d, alpha))
	g		<- fhat.fn(y, cbind(d, alpha), lambda = mix.out$lambda, beta = mix.out$beta, sigma = mix.out$sigma)
	e		<- -100									# Small negative values for not well-define or -Inf when taking log
	
	### Evaluation step ###
	# Calculate the log likelihood function for y
	l0 		<- log(p0)
	lp.mo	<- log(p.mo)
	lp.mo	<- ifelse(lp.mo == -Inf | is.na(lp.mo), e, lp.mo)
	l.cont	<- log(g)
	l.cont	<- ifelse(l.cont == -Inf | is.na(l.cont), e, l.cont)
	ll		<- 	1*(y == 0)*l0 + 
				1*(y==MO)*lp.mo + 
				1*(y>0 & y!= MO)*l.cont
	
	# Log liklihood of b1, ..., b_{C-1} conditional on y
	# P(b1,...,b_{C-1}|y) = |J| * prod_{c:b_c=0}exp(V_c)/(1 + exp(V_c)) * prod_{c:b_c>0}f(V_c)
	# where |J| = prod_{c:b_c>0} 1/(b_c + 1) is the Jacobean matrix 
	#		V_c = -X_c*beta + log(b_c+1) + log(1 + 1(y>MO)c_search'(y)), 
	# 		f(x) = exp(-x)/(1 + exp(-x))^2, pdf of logistic distribution 
	e_mat1	<- as.matrix(e_mat[,-C], ncol = C-1)
	e_sgn1	<- 1*(e_mat1 >0)
	dsdy	<- 1*(y > MO) * alpha/(y-MO+1)^2				# 1(y>MO)c_search'(y)
	V		<- -d[,-C] + log(e_mat1 +1) + log(1 + dsdy)		# V_c = -X_c*beta + log(b_c+1) + log(1 + 1(y>MO)c_search'(y))
	logJ	<- rowSums(-log(e_mat1 + 1)*e_sgn1)				# log(Jacobean) = sum_{c: e_c>0} -log(b_c + 1)
	ll		<- ll + logJ + 
				rowSums((1 - e_sgn1) * (V - log(1 + exp(V))) ) + 					
			    rowSums( e_sgn1 	 * ( -V - 2*log(1 + exp(-V))) )
	return(ll)
}

# Check the evaluation time
# Although simulation takes a long time, log likelihood evaluation is very fast 
system.time(tmp <- ll.fn(c(beta, log(alpha))))

# Run MLE optimization
beta.init	<- c(beta, log(alpha))
sol		<- maxLik(ll.fn, start = beta.init, method = "NM")
summary(sol)
matrix(cbind(c(beta, log(alpha)), coef(sol)), 4,2, dimnames = list(c("beta1", "beta2", "beta3", "log(alpha)"), c("true", "estimate")))

beta.init	<- c(beta - .2, log(.9))
sol		<- maxLik(ll.fn, start = beta.init, method = "NM")
summary(sol)
matrix(cbind(c(beta, log(alpha)), coef(sol)), 4,2, dimnames = list(c("beta1", "beta2", "beta3","log(alpha)"), c("true", "estimate")))

# Check the shape of log likelihood function
# Fix alpha, and beta[1] and beta[3]
tmp	<- seq(-1, 1, .05) + beta[3]
tmp1	<- sapply(tmp, function(x) sum(ll.fn(c(beta[1], beta[2], x, log(alpha)))))
plot(tmp, tmp1, xlab = "beta3", ylab = "Loglikelihood")
lines(tmp, tmp1)
abline(v = beta[3], col = "red")

# Fix beta
tmp	<- seq(-2, 2, .05) + log(alpha)
tmp1	<- sapply(tmp, function(x) sum(ll.fn(c(beta, x))))
plot(tmp, tmp1, xlab = "log(alpha)", ylab = "Loglikelihood")
lines(tmp, tmp1)
abline(v = log(alpha), col = "red")



