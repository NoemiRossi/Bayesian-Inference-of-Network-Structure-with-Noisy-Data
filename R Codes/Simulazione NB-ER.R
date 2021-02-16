library(network)
library(Matrix)
library(igraph)
library(networkD3)
library(plot.matrix)
library(lattice)

set.seed(2406)

### ------------------------------------------------------------------------------
### 1. Simulation
###    Suppose: X_ij|A_ij=1 ~ NegBin(mu=14, phi=50)
###                X_ij|A_ij=0 ~ Negbin(mu=2, phi=25)
### ------------------------------------------------------------------------------

# Generate the true adjacency matrix. Suppose: A_ij ~ Be(p), n = #nodi
n = 20
p = 0.2
A = forceSymmetric(Matrix(rbinom(n*n,1,p),n))
A = as.matrix(A)
isSymmetric(A)


# Plot
cols2 = colorRampPalette(c("white","black"))(256)
plot(x, asp=T, col=cols2)

# Simulation 1
sim1 = as.matrix(Matrix(0,n,n))
for(i in 1:(n-1)){
  for(j in (i+1):n){
    if(A[i,j]==1){
      sim1[i,j] = rnbinom(1,50,mu=14)
    }
    if(A[i,j]==0){
      sim1[i,j] = rnbinom(1,25,mu=2)
    }
    sim1[j,i] = sim1[i,j]
  }
}
sim1

# Plot simulation 1 and matrix
par(mfrow=c(1,2))
plot(A, asp=T, col=cols2, main='Underlying Network')
plot(sim1, asp=T, col=cols2, main='Simulation 1')

### ------------------------------------------------------------------------------
### 2. Stan for the posterior
### ------------------------------------------------------------------------------
library(rstan)

n = dim(sim1)[1]

X_data <- list(n = n, X = sim1, rates_std_prior = c(100,100), rho_prior = c(1,1))

fit = stan(
  file    = "Stan/NegBin-ER.stan",   # Stan program
  data    = X_data,         # named list of data
  chains  = 1,              # number of Markov chains
  warmup  = 1000,           # number of warmup iterations per chain
  iter    = 4000,           # total number of iterations per chains
  cores   = 1,              # number of cores (could use one per chain)
  refresh = 0               # no progress shown
)

observations = as.data.frame(fit)
mu_0      = mean(observations$`mu[1]`)
mu_1      = mean(observations$`mu[2]`)
phi_0      = mean(observations$`phi[1]`)
phi_1      = mean(observations$`phi[2]`)
rho          = mean(observations$rho)
Q            = matrix(0,n,n)

for(i in 1:n){
  for(j in 1:n){
    var = paste('Q[',i,',',j,']', sep='')
    Q[i,j] = mean(as.matrix(observations[var]))
  }
}

par(mfrow=c(1,3))
plot(sim1, asp=T, col=cols2, main='Dataset originale', border=NA)
plot(A, asp=T, col=cols2, main='Underlying Network', border=NA)
plot(Q, asp=T, col=cols2, main='Posterior', border=NA)

##################
#   Traceplots   #
##################

par(mfrow=c(2,3))
plot(ts(observations$`mu[1]`), ylab='', main='Mu_0')
abline(h=2, col='red', lwd=5)
plot(ts(observations$`mu[2]`), ylab='', main='Mu_1')
abline(h=14, col='red', lwd=5)
plot(ts(observations$`phi[1]`), ylab='', main='Phi_0')
abline(h=25, col='red', lwd=5)
plot(ts(observations$`phi[2]`), ylab='', main='Phi_1')
abline(h=50, col='red', lwd=5)
plot(ts(observations$rho), ylab='', main='Rho')
abline(h=0.2, col='red', lwd=5)

#################
#   Densplots   #
#################


par(mfrow=c(2,3))
densplot(mcmc(observations$`mu[1]`), main='Mu_0')
abline(v=2, col='red', lwd=5)
densplot(mcmc(observations$`mu[2]`), main='Mu_1')
abline(v=14, col='red', lwd=5)
densplot(mcmc(observations$`phi[1]`), main='Phi_0')
abline(v=25, col='red', lwd=5)
densplot(mcmc(observations$`phi[2]`), main='Phi_1')
abline(v=50, col='red', lwd=5)
densplot(mcmc(observations$rho), main='Rho')
abline(v=0.2, col='red', lwd=5)

######################
#   Autocorrelation  #
######################

par(mfrow=c(2,3))
acf(observations$`mu[1]`, main = 'Mu_0')
acf(observations$`mu[2]`, main = 'Mu_1')
acf(observations$`phi[1]`, main = 'Phi_0')
acf(observations$`phi[2]`, main = 'Phi_1')
acf(observations$rho, main = 'Rho')

###################
#  Post-Stan:CI   #
###################

par(mfrow=c(2,3))
plot(density(observations$`mu[1]`, adj=2), xlab='', main='Mu_0')
abline(v=quantile(observations$`mu[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mean(observations$`mu[1]`), lwd=2, col='red')

plot(density(observations$`mu[2]`, adj=2), xlab='', main='Mu_1')
abline(v=quantile(observations$`mu[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mean(observations$`mu[2]`), lwd=2, col='red')

plot(density(observations$`phi[1]`, adj=2), xlab='', main='Phi_0')
abline(v=quantile(observations$`phi[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mean(observations$`phi[1]`), lwd=2, col='red')

plot(density(observations$`phi[2]`, adj=2), xlab='', main='Phi_1')
abline(v=quantile(observations$`phi[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mean(observations$`phi[2]`), lwd=2, col='red')

plot(density(observations$rho, adj=2), xlab='', main='Rho')
abline(v=quantile(observations$rho, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mean(observations$rho), lwd=2, col='red')

