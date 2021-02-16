library(network)
library(Matrix)
library(igraph)
library(networkD3)
library(plot.matrix)
library(lattice)
library(rstan)

#MODEL:
#Xij|Aij=1 ~ NegBin(mu=30, phi=0.9)
#Xij|Aij=0 ~ NegBin(mu=5, phi=0.1)

mu0=5
mu1=30
phi0=0.1
phi1=0.9

#FOR Aij|theta I want to implement Stochastic Block Model
#P(Aij=1|theta)=w_gi_gj
#P(Aij=0|theta)=1-w_gi_gj


#Consider 20 nodes- > 400 elements
n=20
A<-matrix(1:(n*n), nrow=n, ncol=n)

#Omega is symmetric
omega<-rbind(c(0.8,0.3,0.1),c(0.3,0.7,0.2),c(0.1,0.2,0.8))
isSymmetric(omega)


#Consider 3 balanced classes 
nodes_classes=rbind(c(1,1),c(2,1),c(3,1),c(4,1),c(5,1),c(6,1),c(7,1),c(8,2),c(9,2),c(10,2), c(11,2),c(12,2),c(13,2),c(14,2),c(15,3),c(16,3),c(17,3),c(18,3),c(19,3),c(20,3))
set.seed(1234)

#Generate adjecency matrix A 
for(i in 1:n){
  for(j in i:n){
    if(i==j)
      A[i,j]<-0
    else
      A[j,i]<-A[i,j]<-rbinom(1,1,omega[nodes_classes[i,2],nodes_classes[j,2]])
    
  }   
}

isSymmetric(A)

#Generate simulation1

sim1<-matrix(1:(n*n), nrow=n, ncol=n)
for(i in 1:n){
  for(j in i:n){
    if(i==j)
      sim1[i,j]<-0
    else{
    if(A[i,j]==1)
      sim1[j,i]<-sim1[i,j]<-rnbinom(1,phi1, mu=mu1)
    else
      sim1[j,i]<-sim1[i,j]<-rnbinom(1,phi0, mu=mu0)
    }
  }
}

isSymmetric(sim1)

cols2 = colorRampPalette(c("white","black"))(256)

par(mfrow=c(1,2))
plot(sim1, asp=T, col=cols2, border=NA, main='Simulated Data')
plot(A, asp=T, col=cols2, border=NA, main='Underlying Network')

X_data <- list(
  n=n, #number of nodes
  k=3, #number of classes 
  X=sim1, #data 
  C=nodes_classes, #node-classes
  rates_std_prior=c(10,10)
)


fit = stan(
  file   = "Stan/NegBin-SBM.stan",   # Stan program
  data    = X_data,         # named list of data
  chains  = 1,              # number of Markov chains
  warmup  = 1000,           # number of warmup iterations per chain
  iter    = 4000,           # total number of iterations per chains
  cores   = 1,              # number of cores (could use one per chain)
  refresh = 1               # show progress
)


observations = as.data.frame(fit)
mu_0      = mean(observations$`mu[1]`)
mu_1      = mean(observations$`mu[2]`)
phi_0      = mean(observations$`phi[1]`)
phi_1      = mean(observations$`phi[2]`)

Q            = matrix(0,n,n)

for(i in 1:n){
  for(j in 1:n){
    var = paste('Q[',i,',',j,']', sep='')
    Q[i,j] = mean(as.matrix(observations[var]))
  }
}

fit_estimated_omega=matrix(1:9, nrow=3)
fit_estimated_omega[1,1]=mean(observations$`omega[1,1]`)
fit_estimated_omega[2,2]=mean(observations$`omega[2,2]`)
fit_estimated_omega[3,3]=mean(observations$`omega[3,3]`)
fit_estimated_omega[1,2]<-fit_estimated_omega[2,1]<-mean(observations$`omega[1,2]`)
fit_estimated_omega[1,3]<-fit_estimated_omega[3,1]<-mean(observations$`omega[1,3]`)
fit_estimated_omega[2,3]<-fit_estimated_omega[3,2]<-mean(observations$`omega[2,3]`)
fit_estimated_omega


######################################
# Plot of Data and Adjecency Matrix  #
######################################
par(mfrow=c(1,3))
plot(sim1, asp=T, col=cols2, main='Original Dataset', border=NA)
plot(A, asp=T, col=cols2, main='Posterior', border=NA)
plot(round(Q), asp=T, col=cols2, main='Rounded-Posterior', border=NA)


##############
# Traceplots #
##############

par(mfrow=c(2,5))
plot(ts(observations$`mu[1]`), ylab='', main='Sim1: mu_0')
abline(h=mu0, col='red', lwd=5)
plot(ts(observations$`mu[2]`), ylab='', main='Sim1: mu_1')
abline(h=mu1, col='red', lwd=5)
plot(ts(observations$`phi[1]`), ylab='', main='Sim1: phi_0')
abline(h=phi0, col='red', lwd=5)
plot(ts(observations$`phi[2]`), ylab='', main='Sim1: phi_1')
abline(h=phi1, col='red', lwd=5)
plot(ts(observations$`omega[1,1]`), ylab='', main='Sim1: omega[1,1]')
abline(h=omega[1,1], col='red', lwd=5)
plot(ts(observations$`omega[2,2]`), ylab='', main='Sim1: omega[2,2]')
abline(h=omega[2,2], col='red', lwd=5)
plot(ts(observations$`omega[3,3]`), ylab='', main='Sim1: omega[3,3]')
abline(h=omega[3,3], col='red', lwd=5)
plot(ts(observations$`omega[1,2]`), ylab='', main='Sim1: omega[1,2]')
abline(h=omega[1,2], col='red', lwd=5)
plot(ts(observations$`omega[1,3]`), ylab='', main='Sim1: omega[1,3]')
abline(h=omega[1,3], col='red', lwd=5)
plot(ts(observations$`omega[2,3]`), ylab='', main='Sim1: omega[2,3]')
abline(h=omega[2,3], col='red', lwd=5)

##############
# Denseplot  #
##############
par(mfrow=c(2,5))
densplot(mcmc(observations$`mu[1]`), ylab='', main='Sim1: mu_0')
abline(v=mu0, col='red', lwd=5)
densplot(mcmc(observations$`mu[2]`), ylab='', main='Sim1: mu_1')
abline(v=mu1, col='red', lwd=5)
densplot(mcmc(observations$`phi[1]`), ylab='', main='Sim1: phi_0')
abline(v=phi0, col='red', lwd=5)
densplot(mcmc(observations$`phi[2]`), ylab='', main='Sim1: phi_1')
abline(v=phi1, col='red', lwd=5)
densplot(mcmc(observations$`omega[1,1]`), ylab='', main='Sim1: omega[1,1]')
abline(v=omega[1,1], col='red', lwd=5)
densplot(mcmc(observations$`omega[2,2]`), ylab='', main='Sim1: omega[2,2]')
abline(v=omega[2,2], col='red', lwd=5)
densplot(mcmc(observations$`omega[3,3]`), ylab='', main='Sim1: omega[3,3]')
abline(v=omega[3,3], col='red', lwd=5)
densplot(mcmc(observations$`omega[1,2]`), ylab='', main='Sim1: omega[1,2]')
abline(v=omega[1,2], col='red', lwd=5)
densplot(mcmc(observations$`omega[1,3]`), ylab='', main='Sim1: omega[1,3]')
abline(v=omega[1,3], col='red', lwd=5)
densplot(mcmc(observations$`omega[2,3]`), ylab='', main='Sim1: omega[2,3]')
abline(v=omega[2,3], col='red', lwd=5)

###################
# Autocorrelation #
###################

par(mfrow=c(2,5))
acf(observations$`mu[1]`, ylab='', main='Sim1: mu_0')
acf(observations$`mu[2]`, ylab='', main='Sim1: mu_1')
acf(observations$`phi[1]`, ylab='', main='Sim1: phi_0')
acf(observations$`phi[2]`, ylab='', main='Sim1: phi_1')
acf(observations$`omega[1,1]`, ylab='', main='Sim1: omega[1,1]')
acf(observations$`omega[2,2]`, ylab='', main='Sim1: omega[2,2]')
acf(observations$`omega[3,3]`, ylab='', main='Sim1: omega[3,3]')
acf(observations$`omega[1,2]`, ylab='', main='Sim1: omega[1,2]')
acf(observations$`omega[1,3]`, ylab='', main='Sim1: omega[1,3]')
acf(observations$`omega[2,3]`, ylab='', main='Sim1: omega[2,3]')

####################
#  Post-Stan: CI   #
####################

par(mfrow=c(2,5))
densplot(mcmc(observations$`mu[1]`), ylab='', main='Sim1: mu_0')
abline(v=quantile(observations$`mu[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mu0, col='red', lwd=5)
densplot(mcmc(observations$`mu[2]`), ylab='', main='Sim1: mu_1')
abline(v=quantile(observations$`mu[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mu1, col='red', lwd=5)
densplot(mcmc(observations$`phi[1]`), ylab='', main='Sim1: phi_0')
abline(v=quantile(observations$`phi[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=phi0, col='red', lwd=5)
densplot(mcmc(observations$`phi[2]`), ylab='', main='Sim1: phi_1')
abline(v=quantile(observations$`phi[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=phi1, col='red', lwd=5)
densplot(mcmc(observations$`omega[1,1]`), ylab='', main='Sim1: omega[1,1]')
abline(v=quantile(observations$`omega[1,1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=omega[1,1], col='red', lwd=5)
densplot(mcmc(observations$`omega[2,2]`), ylab='', main='Sim1: omega[2,2]')
abline(v=quantile(observations$`omega[2,2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=omega[2,2], col='red', lwd=5)
densplot(mcmc(observations$`omega[3,3]`), ylab='', main='Sim1: omega[3,3]')
abline(v=quantile(observations$`omega[3,3]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=omega[3,3], col='red', lwd=5)
densplot(mcmc(observations$`omega[1,2]`), ylab='', main='Sim1: omega[1,2]')
abline(v=quantile(observations$`omega[1,2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=omega[1,2], col='red', lwd=5)
densplot(mcmc(observations$`omega[1,3]`), ylab='', main='Sim1: omega[1,3]')
abline(v=quantile(observations$`omega[1,3]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=omega[1,3], col='red', lwd=5)
densplot(mcmc(observations$`omega[2,3]`), ylab='', main='Sim1: omega[2,3]')
abline(v=quantile(observations$`omega[2,3]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=omega[2,3], col='red', lwd=5)

