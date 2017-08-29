# Conditional density estimation 
# y_i ~ sum_j lambda_j*Normal(X_i*beta_j, sigma_j)

library(np)
library(mixtools)
library(hdrcde)

# Set parameters
K 	<- 3
lambda	<- c(.1, .5, .4)
beta	<- matrix(c(0, -1, -3, 0, 2, -1, 0, 5, -.5), ncol= K)
sigma	<- c(.5, 1, 1)

# Simulate data
N		<- 2000
X		<- cbind(1, matrix(runif(N*2, -1, 1), N, 2))
mu		<- apply(beta, 2, function(b) X %*%b)
idx		<- sample(1:K, N, prob = lambda, replace = T)
y 		<- rep(NA, N)
for(i in 1:N){
	y[i]<- rnorm(1, mu[i,idx[i]], sigma[idx[i]])
}
hist(y)
simdat	<- data.frame(X[,-1], y)
newdat	<- simdat[1:100,]

# Mixture normal estimation 
mix.out <- regmixEM(y= y, x= X[,-1], verb = TRUE, epsilon = 1e-04, k = K) 
beta
mix.out$beta

# Nonparametric estimation 
np.out	<- npcdens(y ~ X1 + X2, data = simdat)
plot(np.out)

# Compare the density estimation 
muhat		<- apply(mix.out$beta, 2, function(b) X %*%b)
sigmahat	<- mix.out$sigma
lambda.hat	<- mix.out$lambda
ggtmp	<- cbind(simdat, d.np = np.out$condens)
ggtmp$d.mix	<- rowSums(sapply(1:K, function(i) lambda.hat[i]*dnorm(y, muhat[,i], sigmahat[i])))
ggtmp$d.true<- sapply(1:N, function(i) dnorm(y[i], mu[i,idx[i]], sigma[idx[i]]))
cor(ggtmp[,c("d.true", "d.mix", "d.np")])

