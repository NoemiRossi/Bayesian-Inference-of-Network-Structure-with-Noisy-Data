library(network)
library(Matrix)
library(igraph)
library(networkD3)
library(plot.matrix)
library(lattice)
library(rstan)
library(coda)

### ------------------------------------------------------------------------------
### 1. Simulation
###    Suppose: Xij|Aij=2 ~ Poisson(lambda2)
###                Xij|Aij=1 ~ Poisson(lambda1)
###                Xij|Aij=0 ~ Poisson(lambda0)
### ------------------------------------------------------------------------------

lamb_2  = 100
lamb_1  = 20
lamb_0  = 1

# Generate the true adjacency matrix. Suppose: A_ij|gruppo_j ~ Be(p_j), n = #nodi
n  = 20
p2 = 0.1
p1 = 0.3
p0 = 1-p1-p2

# Build nxn Matrix
x = forceSymmetric(Matrix(0,n,n))
x = as.matrix(x)
diag(x) = 0*1:length(diag(x))

for(i in 1:(n-1)){
  for(j in (i+1):n){
    x[i,j] = sample(c(0,1,2), 1, prob=c(p0, p1, p2), replace=TRUE)
    x[j,i] = x[i,j]
  }
}
isSymmetric(x)

# Plot of the matrix
cols2 = colorRampPalette(c("white","black"))(256)
plot(x, asp=T, col=cols2)

# Simulation 1
sim1 = as.matrix(Matrix(0,n,n))
for(i in 1:(n-1)){
  for(j in (i+1):n){
    if(x[i,j]==2){
      sim1[i,j] = rpois(1,lamb_2)
    }
    if(x[i,j]==1){
      sim1[i,j] = rpois(1,lamb_1)
    }
    if(x[i,j]==0){
      sim1[i,j] = rpois(1,lamb_0)
    }
    sim1[j,i] = sim1[i,j]
  }
}
sim1

# Plot
par(mfrow=c(1,2))
plot(x, asp=T, col=cols2, main='Underlying Network')
plot(sim1, asp=T, col=cols2, main='Simulation 1')

### ------------------------------------------------------------------------------
### 2. Stan for the posterior
### ------------------------------------------------------------------------------

X_data <- list(
  n = n, 
  X = sim1, 
  T = 3,
  rates_std_prior = c(100,100,100)
)

fit = stan(
  file   = "Stan/Poisson_multi.stan", # Stan program
  data    = X_data,              # named list of data
  chains  = 1,                   # number of Markov chains
  warmup  = 1000,                # number of warmup iterations per chain
  iter    = 4000,                # total number of iterations per chains
  cores   = 1,                   # number of cores (could use one per chain)
  refresh = 1                    # show progress
)

observations = as.data.frame(fit)
rates_1 = mean(observations$`rates[1]`)
rates_2 = mean(observations$`rates[2]`)
rates_3 = mean(observations$`rates[3]`)
rho_1   = mean(observations$`rho[1]`)
rho_2   = mean(observations$`rho[2]`)
rho_3   = mean(observations$`rho[3]`)

Q = matrix(0,n,n)

for(i in 1:n){
  for(j in 1:n){
    var1  = paste('Q[',i,',',j,',',1,']', sep='')
    var2  = paste('Q[',i,',',j,',',2,']', sep='')
    var3  = paste('Q[',i,',',j,',',3,']', sep='')
    temp0 = mean(as.matrix(observations[var1]))
    temp1 = mean(as.matrix(observations[var2]))
    temp2 = mean(as.matrix(observations[var3]))
    if(max(temp0, temp1, temp2)==temp0){
      Q[i,j] = 0
    }
    if(max(temp0, temp1, temp2)==temp1){
      Q[i,j] = 1
    }
    if(max(temp0, temp1, temp2)==temp2){
      Q[i,j] = 2
    }
  }
}

par(mfrow=c(1,3))
plot(x, asp=T, col=cols2, main='Underlying Network', border=TRUE)
plot(sim1, asp=T, col=cols2, main='Simulation', border=TRUE)
plot(Q-2*diag(dim(Q)[1]), asp=T, col=cols2, main='Posterior', border=TRUE)

##############
# Traceplots #
##############

par(mfrow=c(2,3))
plot(ts(observations$`rates[1]`), ylab='', main='lambda_0')
abline(h=lamb_0, lwd=5, col='red')
plot(ts(observations$`rates[2]`), ylab='', main='lambda_1')
abline(h=lamb_1, lwd=5, col='red')
plot(ts(observations$`rates[3]`), ylab='', main='lambda_2')
abline(h=lamb_2, lwd=5, col='red')
plot(ts(observations$`rho[1]`), ylab='', main='rho_0')
abline(h=p0, lwd=5, col='red')
plot(ts(observations$`rho[2]`), ylab='', main='rho_1')
abline(h=p1, lwd=5, col='red')
plot(ts(observations$`rho[3]`), ylab='', main='rho_2')
abline(h=p2, lwd=5, col='red')


##############
# Denseplot  #
##############

par(mfrow=c(2,3))
densplot(mcmc(observations$`rates[1]`), ylab='', main='lambda_0')
abline(v=lamb_0, col='red', lwd=2)
densplot(mcmc(observations$`rates[2]`), ylab='', main='lambda_1')
abline(v=lamb_1, col='red', lwd=2)
densplot(mcmc(observations$`rates[3]`), ylab='', main='lambda_2')
abline(v=lamb_2, col='red', lwd=2)
densplot(mcmc(observations$`rho[1]`), ylab='', main='rho_0')
abline(v=p0, col='red', lwd=2)
densplot(mcmc(observations$`rho[2]`), ylab='', main='rho_1')
abline(v=p1, col='red', lwd=2)
densplot(mcmc(observations$`rho[3]`), ylab='', main='rho_2')
abline(v=p2, col='red', lwd=2)



###################
# Autocorrelation #
###################

par(mfrow=c(2,3))
acf(observations$`rates[1]`, ylab='', main='Sim1: lambda_0')
acf(observations$`rates[2]`, ylab='', main='Sim1: lambda_1')
acf(observations$`rates[3]`, ylab='', main='Sim1: lambda_2')
acf(observations$`rho[1]`, ylab='', main='Sim1: rho_0')
acf(observations$`rho[2]`, ylab='', main='Sim1: rho_1')
acf(observations$`rho[3]`, ylab='', main='Sim1: rho_2')


#####################
#   Post-Stan: CI  #
####################

par(mfrow=c(2,3))
densplot(mcmc(observations$`rates[1]`), ylab='', main='lambda_0')
abline(v=quantile(observations$`rates[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=lamb_0, col='red', lwd=5)
densplot(mcmc(observations$`rates[2]`), ylab='', main='lambda_1')
abline(v=quantile(observations$`rates[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=lamb_1, col='red', lwd=5)
densplot(mcmc(observations$`rates[3]`), ylab='', main='lambda_2')
abline(v=quantile(observations$`rates[3]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=lamb_2, col='red', lwd=5)
densplot(mcmc(observations$`rho[1]`), ylab='', main='rho_0')
abline(v=quantile(observations$`rho[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=p0, col='red', lwd=5)
densplot(mcmc(observations$`rho[2]`), ylab='', main='rho_1')
abline(v=quantile(observations$`rho[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=p1, col='red', lwd=5)
densplot(mcmc(observations$`rho[3]`), ylab='', main='rho_2')
abline(v=quantile(observations$`rho[3]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=p2, col='red', lwd=5)


