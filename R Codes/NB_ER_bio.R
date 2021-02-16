library(network)
library(Matrix)
library(igraph)
library(networkD3)
library(plot.matrix)
library(lattice)
library(coda)

setwd("")

### ------------------------------------------------------------------------------
### 1. Generate the matrix
### ------------------------------------------------------------------------------
data_bio    = read.csv('Datasets/new_dataset1_bio.csv')
data_bio$V4 = as.factor(data_bio$V4)
data_bio$V5 = as.factor(data_bio$V5)

# Delete NA
data_bio = data_bio[-which(is.na(data_bio$V1)),]
summary(data_bio)

# Index
indici = sort(unique(c(data_bio$V1, data_bio$V2)))
indici
n = length(indici)
n

# Matrix nxn
data = matrix(0, nrow=n, ncol=n)
data = as.data.frame(data)
colnames(data) = indici
rownames(data) = indici

for(k in 1:dim(data_bio)[1]){
  val_i = data_bio$V1[k]
  val_j = data_bio$V2[k]
  data[which(rownames(data)==val_i), which(colnames(data)==val_j)] = data_bio$V3[k]
  data[which(rownames(data)==val_j), which(colnames(data)==val_i)] = data_bio$V3[k]
}

data = as.matrix(data)

### ------------------------------------------------------------------------------
### 2. Plot of the Matrix
### ------------------------------------------------------------------------------
cols2 = colorRampPalette(c("white","black"))(256)
plot(data, asp=T, col=cols2, border=NA)

### ------------------------------------------------------------------------------
### 3. Stan for the Posterior
### ------------------------------------------------------------------------------
library(rstan)

X_data <- list(n = n, X = data, rates_std_prior = c(100,100), rho_prior = c(1,1))

fit = stan(
  file    = "Stan/NegBin-ER.stan",   # Stan program
  data    = X_data,         # named list of data
  chains  = 1,              # number of Markov chains
  warmup  = 1000,           # number of warmup iterations per chain
  iter    = 2000,           # total number of iterations per chains
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
plot(data, asp=T, col=cols2, main='Dataset originale', border=NA)
plot(Q, asp=T, col=cols2, main='Posterior', border=NA)

################
#  Traceplots  #
################

par(mfrow=c(2,3))
plot(ts(observations$`mu[1]`), ylab='', main='Mu_0')
plot(ts(observations$`mu[2]`), ylab='', main='Mu_1')
plot(ts(observations$`phi[1]`), ylab='', main='Phi_0')
plot(ts(observations$`phi[2]`), ylab='', main='Phi_1')
plot(ts(observations$rho), ylab='', main='Rho')

################
#  Denseplots  #
################

par(mfrow=c(2,3))
densplot(mcmc(observations$`mu[1]`), main='Mu_0')
densplot(mcmc(observations$`mu[2]`), main='Mu_1')
densplot(mcmc(observations$`phi[1]`), main='Phi_0')
densplot(mcmc(observations$`phi[2]`), main='Phi_1')
densplot(mcmc(observations$rho), main='Rho')

####################
# Autocorrelation  #
####################

par(mfrow=c(2,3))
acf(observations$`mu[1]`, main = 'Mu_0')
acf(observations$`mu[2]`, main = 'Mu_1')
acf(observations$`phi[1]`, main = 'Phi_0')
acf(observations$`phi[2]`, main = 'Phi_1')
acf(observations$rho, main = 'Rho')

###################
#  Post Stan: CI  #
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

### ------------------------------------------------------------------------------
### 4. Posterior Predictive Assessment with Discrepancy 
### ------------------------------------------------------------------------------

discrepancy<-function(X,Q,mu0,mu1,phi0,phi1){
  
  #Build X_star
  X_star=matrix(0,nrow=n,ncol=n)
  for(i in 1: n){
    for(j in i:n){     
      X_star[i,j]=Q[i,j]*mu1+(1-Q[i,j])*mu0
    }
  }
  
  #Compute the Discrepancy
  somma=0
  for(i in 1:n){
    for(j in i:n){
      if(X[i,j]!=0 && X_star[i,j]!=0)
        somma=somma+X[i,j]*log(X[i,j]/X_star[i,j])
    }
  }
  return(somma)
}

d_data=numeric(500)
d_artificial=numeric(500)

for (param_id in 1:500){
  print(param_id)
  Qq=matrix(nrow=n,ncol=n)
  A=matrix(nrow=n,ncol=n)
  
  #Build A
  for(i in 1:n){
    for(j in i:n){
      var = paste('Q[',i,',',j,']', sep='')
      Qq[i,j] = as.matrix(observations[var])[param_id,1]
      A[i,j]=as.numeric(runif(1,0,1)<Qq[i,j])
      A[j,i]=A[i,j]
    }
  }
  
  #Build X_tilde with replicated data
  X_tilde_0=matrix(rnbinom(n*n,observations["phi[1]"][param_id,1],mu=observations["mu[1]"][param_id,1]),nrow=n,ncol=n)
  X_tilde_1=matrix(rnbinom(n*n,observations["phi[2]"][param_id,1],mu=observations["mu[2]"][param_id,1]),nrow=n,ncol=n)
  X_tilde=(1-A)*X_tilde_0+A*X_tilde_1
  d_data[param_id]=discrepancy(data,Qq,observations["mu[1]"][param_id,1],observations["mu[2]"][param_id,1],observations["phi[1]"][param_id,1],observations["phi[2]"][param_id,1])
  d_artificial[param_id]=discrepancy(X_tilde,Qq,observations["mu[1]"][param_id,1],observations["mu[2]"][param_id,1],observations["phi[1]"][param_id,1],observations["phi[2]"][param_id,1])
}

#Scatter plot
par(mfrow=c(1,1))
plot(d_data[d_data<d_artificial],d_artificial[d_data<d_artificial], col='red', xlab='D(X|theta)', ylab='D(X_tilde|theta)', xlim=c(5000,56000), ylim=c(5000,56000))
points(d_data[d_data>d_artificial],d_artificial[d_data>d_artificial], col='blue')
abline(a=0,b=1)
pp_p_value<-length(d_artificial[d_data<d_artificial])/500
pp_p_value

### ------------------------------------------------------------------------------
### 5. Plot of the Network  
### ------------------------------------------------------------------------------

data_bio$V4 = as.factor(as.numeric(data_bio$V4))
data_bio$V5 = as.factor(as.numeric(data_bio$V5))
summary(data_bio)

classi = as.data.frame(indici)
classi[,'nuovo_numero'] = 1:length(indici)
classi[,'classe'] = 0*(1:length(indici))

for(k in 1:dim(data_bio)[1]){
  val_i = data_bio$V1[k]
  val_j = data_bio$V2[k]
  categoria_i = data_bio$V4[k]
  categoria_j = data_bio$V5[k]
  
  classi$classe[which(classi$indici==val_i)] = categoria_i
  classi$classe[which(classi$indici==val_j)] = categoria_j
}

classi$classe = as.factor(classi$classe)

num_elem = 110
matrice_plot = data[1:num_elem, 1:num_elem]
matrice_plot_posterior = round(Q[1:num_elem, 1:num_elem])

colori = as.numeric(classi$classe[1:num_elem])
colori[which(colori==1)] = 'tomato'
colori[which(colori==2)] = 'royalblue3'
colori[which(colori==3)] = 'springgreen3' 

g = graph.adjacency(matrice_plot_posterior, mode="undirected", weighted=NULL)
coords = layout_with_fr(g)

# Plot of the graph
plot(g, layout=coords, vertex.size=10, vertex.color = colori, vertex.label=NA, 
     edge.color='black')

plot(g, vertex.size=10, vertex.color = colori, layout = layout_in_circle,
     vertex.label=NA, edge.color='black')

plot(g, vertex.size=10, vertex.color = colori, layout = layout_nicely,
     edge.color='black', vertex.label=NA)

plot(g, vertex.size=10, vertex.color = colori,
     vertex.label=NA, edge.color='black')


plot(g, vertex.size=5, vertex.color = colori, layout = layout_nicely,
     edge.color='black', vertex.label=NA)

