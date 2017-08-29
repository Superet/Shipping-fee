# Plot utility function 

# Model: 
# U1 = psi_1*log(y+1) + psi_0 (I - y - c_f)
# s.t. y>= 0, I-y >= 0

# U2 = psi_1*log(y+1) + psi_0 (I - y) - c_ship(y)
# s.t. y>= 0, I-y-c_ship(y) >= 0
# and c_ship(y) = c_f if y<MO and c_ship(y) = 0 if y >= MO

# U3 = U2 - c_search(y)
# c_search(y) = 0 if y<MO and alpha_s/(y-MO+1) if y>= MO

library(ggplot2)
library(reshape2)

# Utility functions
u1 <- function(psi0, psi1, y){
	psi1*log(y+1) + psi0*(I-y-cf)
}

u2 <- function(psi0, psi1, alpha, y){
	psi1*log(y+1) + psi0*(I-y) - c.search(alpha, y)
}

c.search <- function(alpha, y){
	alpha/(y - MO + 1)
}

u3 	<- function(psi0, psi1, alpha, y){
	psi1*log(y+1) + psi0*(I-y - 1*(y<MO)*cf) - 1*(y>=MO)*c.search(alpha, y)
}

# Set parameters 
psi0	<- 1
psi1	<- c(5, 10, 25)
alpha	<- c(1, 10)
I 		<- 100					
MO		<- 20
y		<- seq(0, 50, 1)
cf		<- 8

# Calculate utlity along y for each set of parameters
ggtmp	<- data.frame(expand.grid(y, alpha, psi1))
names(ggtmp)	<- c("y","alpha", "psi1")
ggtmp$u1	<- sapply(1:nrow(ggtmp), function(i) u1(psi0, ggtmp[i,"psi1"], ggtmp[i,"y"]))
ggtmp$u2	<- sapply(1:nrow(ggtmp), function(i) u2(psi0, ggtmp[i,"psi1"], ggtmp[i,"alpha"],ggtmp[i,"y"]))
ggtmp$u		<- with(ggtmp, ifelse(y >= MO, u2, u1))
ggtmp1		<- melt(ggtmp[,c("y", "alpha", "psi1","u1","u")], id.vars = c("y", "alpha", "psi1"))
ggtmp1$variable	<- factor(ggtmp1$variable, levels = c("u", "u1"))
ggtmp1$psi1	<- factor(ggtmp1$psi1, levels = psi1, labels = c("small", "large", "larger"))
ggtmp1$alpha<- factor(ggtmp1$alpha, levels = alpha, labels = c("small", "large"))

# Plot the utility function 
ggplot(ggtmp, aes(y, u)) + geom_line() + 
		geom_vline(xintercept = MO, col = "red") + 
		facet_wrap(alpha ~ psi1, labeller = label_both, scale = "free")

# Plot the (discontinuous) utiity function and (continuous) u1
ggplot(ggtmp1, aes(y, value, linetype = variable)) + geom_line() + 
		geom_vline(xintercept = MO, col = "red") + 
		facet_grid(alpha ~ psi1, labeller = label_both)

		

		