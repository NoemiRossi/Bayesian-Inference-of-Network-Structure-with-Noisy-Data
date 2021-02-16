library(network)
library(Matrix)
library(igraph)
library(networkD3)
library(plot.matrix)
library(lattice)
library(Matrix)

### ------------------------------------------------------------------------------
### Simulation
### Suppose: X_ij|A_ij=1 ~ Poisson(14)
###             X_ij|A_ij=0 ~ Poisson(0.5)
### ------------------------------------------------------------------------------


# Generate the true adjacency matrix. Suppose: A_ij ~ Be(p), n = #nodi
n = 20
p = 0.2
A = forceSymmetric(Matrix(rbinom(n*n,1,p),n))
A = as.matrix(A)
isSymmetric(A)

# Plot
cols2 = colorRampPalette(c("white","black"))(256)
plot(x, asp=T, col=c('white','black'))

# Simulation 1
set.seed(1231235)
sim1 = as.matrix(Matrix(0,n,n))
for(i in 1:(n-1)){
  for(j in (i+1):n){
    if(A[i,j]==1){
      sim1[i,j] = rpois(1,14)
    }
    if(A[i,j]==0){
      sim1[i,j] = rpois(1,0.5)
    }
    sim1[j,i] = sim1[i,j]
  }
}
sim1

# Plot simulation 1 and original
par(mfrow=c(1,2))
plot(A, asp=T, col=cols2, main='Underlying Network')
plot(sim1, asp=T, col=cols2, main='Simulation 1')

### ------------------------------------------------------------------------------
### Stan
### ------------------------------------------------------------------------------
library(rstan)

# Simulazione considerata:
sim = sim1

n = dim(sim)[1]

X_data <- list(
  n = n,
  X = sim,
  rates_std_prior = c(100,100),
  rho_prior = c(1,1)
)

fit = stan(
  file    = "Stan/Poisson-ER.stan",   # Stan program
  data    = X_data,         # named list of data
  chains  = 4,              # number of Markov chains
  warmup  = 1000,           # number of warmup iterations per chain
  iter    = 2000,           # total number of iterations per chains
  cores   = 1,              # number of cores (could use one per chain)
  refresh = 0               # no progress shown
)

observations = as.data.frame(fit)
rates_1      = mean(observations$`rates[1]`)
rates_2      = mean(observations$`rates[2]`)
rho          = mean(observations$rho)
Q            = matrix(0,n,n)

for(i in 1:n){
  for(j in 1:n){
    var = paste('Q[',i,',',j,']', sep='')
    Q[i,j] = mean(as.matrix(observations[var]))
  }
}

par(mfrow=c(1,3))
plot(x, asp=T, col=cols2, main='Underlying Network')
plot(sim, asp=T, col=cols2, main='Simulation')
plot(Q, asp=T, col=cols2, main='Posterior')

################
#  Traceplots  #
################

par(mfrow=c(1,3))
plot(ts(observations$`rates[1]`), ylab='', main='lambda_0')
plot(ts(observations$`rates[2]`), ylab='', main='lambda_1')
plot(ts(observations$rho), ylab='', main='rho')

#####################
#  Autocorrelation  #
#####################

par(mfrow=c(1,3))
acf(observations$`rates[1]`, main = 'lambda_0')
acf(observations$`rates[2]`, main = 'lambda_1')
acf(observations$rho, main = 'rho')

####################
#   Post-Stan: CI  #
####################

par(mfrow=c(1,3))
plot(density(observations$`rates[1]`, adj=2), xlab='', main='lambda_0')
abline(v=quantile(observations$`rates[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mean(observations$`rates[1]`), lwd=5, col='red')

plot(density(observations$`rates[2]`, adj=2), xlab='', main='lambda_1')
abline(v=quantile(observations$`rates[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mean(observations$`rates[2]`), lwd=5, col='red')

plot(density(observations$rho, adj=2), xlab='', main='rho')
abline(v=quantile(observations$rho, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mean(observations$rho), lwd=5, col='red')

