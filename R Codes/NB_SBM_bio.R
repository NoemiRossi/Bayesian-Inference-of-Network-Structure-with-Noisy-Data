library(network)
library(Matrix)
library(igraph)
library(networkD3)
library(plot.matrix)
library(lattice)
library(coda)
library(igraph)
library(lsa)

#setwd("")

### ------------------------------------------------------------------------------
### 1. Genererate the matrix of the data
### ------------------------------------------------------------------------------
data_bio    = read.csv('Datasets/new_dataset1_bio.csv')

data_bio$V4 = as.factor(data_bio$V4)
data_bio$V5 = as.factor(data_bio$V5)

# Delete the NA if their are present
data_bio = data_bio[-which(is.na(data_bio$V1)),]
summary(data_bio)

# Create index vector with grouped classes
metadata <- read.delim("Datasets/metadata.txt", header=FALSE)

metadata_bio1<-metadata[metadata$V2=='2BIO1', ]
metadata_bio2<-metadata[metadata$V2=='2BIO2', ]
metadata_bio3<-metadata[metadata$V2=='2BIO3', ]

indici = c(metadata_bio1$V1, metadata_bio2$V1, metadata_bio3$V1)
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
### 2. Plot of the matrix
### ------------------------------------------------------------------------------
cols2 = colorRampPalette(c("white","black"))(256)
plot(data, asp=T, col=cols2, border=NA)

### ------------------------------------------------------------------------------
### 3. Classes
### ------------------------------------------------------------------------------
metadata_bio<-metadata[metadata$V2=='2BIO1' | metadata$V2=='2BIO2' | metadata$V2=='2BIO3',]
classes<-as.data.frame(indici)
for (i in 1:dim(classes)[1]){
  
  if(metadata[metadata$V1==classes[i,1],2]=='2BIO1')
    classes[i,2]=1
  else if(metadata[metadata$V1==classes[i,1],2]=='2BIO2')
    classes[i,2]=2
  else
    classes[i,2]=3
}

### ------------------------------------------------------------------------------
### 3. Stan for the posterior
### ------------------------------------------------------------------------------
library(rstan)

X_data <- list(
  n=n, #number of nodes
  k=3, #number of classes 
  X=data, #data 
  C=classes, #node-classes
  rates_std_prior=c(10,10)
)

fit = stan(
  file    = "Stan/NegBin-SBM.stan",   # Stan program
  data    = X_data,         # named list of data
  chains  = 1,              # number of Markov chains
  warmup  = 1000,           # number of warmup iterations per chain
  iter    = 2000,           # total number of iterations per chains
  cores   = 1,              # number of cores (could use one per chain)
  refresh = 1               
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


par(mfrow=c(1,2))
plot(data, asp=T, col=cols2, main='Original Dataset', border=NA)
plot(round(Q), asp=T, col=cols2, main='Rounded-Posterior', border=NA)


##############
# Traceplots #
##############

par(mfrow=c(2,5))
plot(ts(observations$`mu[1]`), ylab='', main='mu_0')
abline(h=mu_0, col='red', lwd=5)
plot(ts(observations$`mu[2]`), ylab='', main='mu_1')
abline(h=mu_1, col='red', lwd=5)
plot(ts(observations$`phi[1]`), ylab='', main='phi_0')
abline(h=phi_0, col='red', lwd=5)
plot(ts(observations$`phi[2]`), ylab='', main='phi_1')
abline(h=phi_1, col='red', lwd=5)
plot(ts(observations$`omega[1,1]`), ylab='', main='omega[1,1]')
abline(h=fit_estimated_omega[1,1], col='red', lwd=5)
plot(ts(observations$`omega[2,2]`), ylab='', main='omega[2,2]')
abline(h=fit_estimated_omega[2,2], col='red', lwd=5)
plot(ts(observations$`omega[3,3]`), ylab='', main='omega[3,3]')
abline(h=fit_estimated_omega[3,3], col='red', lwd=5)
plot(ts(observations$`omega[1,2]`), ylab='', main='omega[1,2]')
abline(h=fit_estimated_omega[1,2], col='red', lwd=5)
plot(ts(observations$`omega[1,3]`), ylab='', main='omega[1,3]')
abline(h=fit_estimated_omega[1,3], col='red', lwd=5)
plot(ts(observations$`omega[2,3]`), ylab='', main='omega[2,3]')
abline(h=fit_estimated_omega[2,3], col='red', lwd=5)


###################
# Autocorrelation #
###################

par(mfrow=c(2,5))
acf(observations$`mu[1]`, ylab='', main='mu_0')
acf(observations$`mu[2]`, ylab='', main='mu_1')
acf(observations$`phi[1]`, ylab='', main='phi_0')
acf(observations$`phi[2]`, ylab='', main='phi_1')
acf(observations$`omega[1,1]`, ylab='', main='omega[1,1]')
acf(observations$`omega[2,2]`, ylab='', main='omega[2,2]')
acf(observations$`omega[3,3]`, ylab='', main='omega[3,3]')
acf(observations$`omega[1,2]`, ylab='', main='omega[1,2]')
acf(observations$`omega[1,3]`, ylab='', main='omega[1,3]')
acf(observations$`omega[2,3]`, ylab='', main='omega[2,3]')

####################
#   Post-Stan: CI  # 
####################

par(mfrow=c(2,5))
densplot(mcmc(observations$`mu[1]`), ylab='', main='mu_0')
abline(v=quantile(observations$`mu[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mu_0, col='red', lwd=5)
densplot(mcmc(observations$`mu[2]`), ylab='', main='mu_1')
abline(v=quantile(observations$`mu[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=mu_1, col='red', lwd=5)
densplot(mcmc(observations$`phi[1]`), ylab='', main='phi_0')
abline(v=quantile(observations$`phi[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=phi_0, col='red', lwd=5)
densplot(mcmc(observations$`phi[2]`), ylab='', main='phi_1')
abline(v=quantile(observations$`phi[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=phi_1, col='red', lwd=5)
densplot(mcmc(observations$`omega[1,1]`), ylab='', main='omega[1,1]')
abline(v=quantile(observations$`omega[1,1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=fit_estimated_omega[1,1], col='red', lwd=5)
densplot(mcmc(observations$`omega[2,2]`), ylab='', main='omega[2,2]')
abline(v=quantile(observations$`omega[2,2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=fit_estimated_omega[2,2], col='red', lwd=5)
densplot(mcmc(observations$`omega[3,3]`), ylab='', main='omega[3,3]')
abline(v=quantile(observations$`omega[3,3]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=fit_estimated_omega[3,3], col='red', lwd=5)
densplot(mcmc(observations$`omega[1,2]`), ylab='', main='omega[1,2]')
abline(v=quantile(observations$`omega[1,2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=fit_estimated_omega[1,2], col='red', lwd=5)
densplot(mcmc(observations$`omega[1,3]`), ylab='', main='omega[1,3]')
abline(v=quantile(observations$`omega[1,3]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=fit_estimated_omega[1,3], col='red', lwd=5)
densplot(mcmc(observations$`omega[2,3]`), ylab='', main='omega[2,3]')
abline(v=quantile(observations$`omega[2,3]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=fit_estimated_omega[2,3], col='red', lwd=5)

######################################################
#  Histogram of the data with fitted distributions   # 
######################################################

dataframeQ=as.data.frame(cbind(as.vector(round(Q)),as.vector(data)))
A0=dataframeQ[dataframeQ$V1==0,]
A1=dataframeQ[dataframeQ$V1==1,]
A11=A1[A1$V2!=0,]
par(mfrow=c(1,2))
hist(A0$V2, freq=FALSE, breaks=10, xlab='X|A=0', main="Histogram of NegativeBinomial(mu_0,phi_0)")
x=seq(from=0, to=12, by=1)
y=dnbinom(x,phi_0, mu=mu_0)
lines(x,y, col='red')
hist(A11$V2, freq=FALSE, breaks=50, xlim=c(10, 250), ylim=c(0,0.06), xlab='X|A=1', main="Histogram of NegativeBinomial(mu_1,phi_1)")
x=seq(from=0, to=250, by=1)
y=dnbinom(x,phi_1, mu=mu_1)
lines(x,y, col='red')

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
plot(d_data[d_data<d_artificial],d_artificial[d_data<d_artificial], col='red', xlab='D(X|theta)', ylab='D(X_tilde|theta)', xlim=c(15000,30000), ylim=c(15000,30000))
points(d_data[d_data>d_artificial],d_artificial[d_data>d_artificial], col='blue')
abline(a=0,b=1)
pp_p_value<-length(d_artificial[d_data<d_artificial])/500
pp_p_value

### ------------------------------------------------------------------------------
### 5. Plot of the Network  
### ------------------------------------------------------------------------------

classi = as.data.frame(indici)
classi[,'nuovo_numero'] = 1:length(indici)
classi[,'classe'] = 0*(1:length(indici))

metadata_bio$V2=as.factor(metadata_bio$V2)

for(i in 1:dim(classi)[1]){
  classi[i,3]=metadata_bio[metadata_bio$V1==classi[i,1],2]
}

classi$classe = as.factor(classi$classe)


num_elem = 111
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
