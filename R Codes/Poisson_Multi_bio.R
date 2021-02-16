library(network)
library(Matrix)
library(igraph)
library(networkD3)
library(plot.matrix)
library(lattice)
library(rstan)
library(coda)
library(igraph)
library(lsa)

#setwd("")


### ------------------------------------------------------------------------------
### 1. Genero la matrice
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

# Build nxn matrix
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
cols3 = colorRampPalette(c("red","blue"))(256)
plot(data, asp=T, col=cols2, border=NA)

### ------------------------------------------------------------------------------
### 3. Stan for the Posterior
### ------------------------------------------------------------------------------

X_data <- list(
  n = n, 
  X = data, 
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
    temp0 <- mean(as.matrix(observations[var1]))
    temp1 <- mean(as.matrix(observations[var2]))
    temp2 <- mean(as.matrix(observations[var3]))
    
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

par(mfrow=c(1,2))
plot(data, asp=T, col=cols2, main='Original Dataset', border=NA)
plot(Q-2*diag(dim(Q)[1]), asp=T, col=cols2, main='Posterior', border=NA)


##############
# Traceplots #
##############
par(mfrow=c(2,3))
plot(ts(observations$`rates[1]`), ylab='', main='lambda_0')
plot(ts(observations$`rates[2]`), ylab='', main='lambda_1')
plot(ts(observations$`rates[3]`), ylab='', main='lambda_2')
plot(ts(observations$`rho[1]`), ylab='', main='rho_0')
plot(ts(observations$`rho[2]`), ylab='', main='rho_1')
plot(ts(observations$`rho[3]`), ylab='', main='rho_2')

###################
# Autocorrelation #
###################

par(mfrow=c(2,3))
acf((observations$`rates[1]`), ylab='', main='lambda_0')
acf(ts(observations$`rates[2]`), ylab='', main='lambda_1')
acf(ts(observations$`rates[3]`), ylab='', main='lambda_2')
acf(ts(observations$`rho[1]`), ylab='', main='rho_0')
acf(ts(observations$`rho[2]`), ylab='', main='rho_1')
acf(ts(observations$`rho[3]`), ylab='', main='rho_2')

#######################
#   Post-Stan: CI     #
#######################
par(mfrow=c(2,3))
densplot(mcmc(observations$`rates[1]`), ylab='', main='lambda_0')
abline(v=quantile(observations$`rates[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=rates_1, col='red', lwd=5)
densplot(mcmc(observations$`rates[2]`), ylab='', main='lambda_1')
abline(v=quantile(observations$`rates[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=rates_2, col='red', lwd=5)
densplot(mcmc(observations$`rates[3]`), ylab='', main='lambda_2')
abline(v=quantile(observations$`rates[3]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=rates_3, col='red', lwd=5)
densplot(mcmc(observations$`rho[1]`), ylab='', main='rho_0')
abline(v=quantile(observations$`rho[1]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=rho_1, col='red', lwd=5)
densplot(mcmc(observations$`rho[2]`), ylab='', main='rho_1')
abline(v=quantile(observations$`rho[2]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=rho_2, col='red', lwd=5)
densplot(mcmc(observations$`rho[3]`), ylab='', main='rho_2')
abline(v=quantile(observations$`rho[3]`, prob=c(0.025, 0.975)), lwd=2, col='grey')
abline(v=rho_3, col='red', lwd=5)


######################################################
#  Histogram of the data with fitted distributions   # 
######################################################

dataframeQ=as.data.frame(cbind(as.vector(Q_012),as.vector(data)))
A0=dataframeQ[dataframeQ$V1==0,]
A1=dataframeQ[dataframeQ$V1==1,]
A2=dataframeQ[dataframeQ$V1==2,]
A11=A1[A1$V2!=0,]
A22=A2[A2$V2!=0,]
par(mfrow=c(1,3))
hist(A0$V2, freq=FALSE, breaks=10, xlab='X|A=0', main="Histogram of Poisson(lambda_0)")
x=seq(from=0, to=12, by=1)
y=dpois(x, rates_1)
lines(x,y, col='red')
hist(A11$V2, freq=FALSE, breaks=10, xlim=c(0, 50), ylim=c(0,0.125), xlab='X|A=1', main="Histogram of Poisson(lambda_1)")
x=seq(from=0, to=250, by=1)
y=dpois(x, rates_2)
lines(x,y, col='red')
hist(A22$V2, freq=FALSE, breaks=20, xlim=c(40, 270), ylim=c(0,0.05), xlab='X|A=2', main="Histogram of Poisson(lambda_2)")
x=seq(from=0, to=250, by=1)
y=dpois(x, rates_3)
lines(x,y, col='red')

### ------------------------------------------------------------------------------
### 4. Discrepancy
### ------------------------------------------------------------------------------


discrepancy<-function(X, Q0, Q1, lambda0, lambda1, lambda2){
  #Build X_star
  X_star=matrix(nrow=n,ncol=n)
  for(i in 1: n){
    for(j in i:n)
      X_star[i,j]=Q0[i,j]*lambda0+Q1[i,j]*lambda1+(1-Q0[i,j]-Q1[i,j])*lambda2
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
  
  Qq0=matrix(nrow=n,ncol=n)
  Qq1=matrix(nrow=n,ncol=n)
  
  A=matrix(nrow=n,ncol=n)
  
  #Build A from the posterior
  for(i in 1:n){
    for(j in i:n){
      var1  = paste('Q[',i,',',j,',',1,']', sep='')
      var2  = paste('Q[',i,',',j,',',2,']', sep='')
      Qq0[i,j] = as.matrix(observations[var1])[param_id,1]
      Qq1[i,j] = as.matrix(observations[var2])[param_id,1]
      rand=runif(1,0,1)
      A[i,j]=as.numeric(rand<Qq1[i,j]*1+(rand>(Qq1[i,j]+Qq0[i,j]))*2)
      A[j,i]=A[i,j]
    }
  }
  #Build new observation from the posterior
  X_tilde_0=matrix(rpois(n*n, observations["rates[1]"][param_id,1]),nrow=n,ncol=n)
  X_tilde_1=matrix(rpois(n*n, observations["rates[2]"][param_id,1]),nrow=n,ncol=n)
  X_tilde_2=matrix(rpois(n*n, observations["rates[3]"][param_id,1]),nrow=n,ncol=n)
  X_tilde=as.numeric(A==0)*X_tilde_0+as.numeric(A==1)*X_tilde_1+as.numeric(A==2)*X_tilde_2
  d_data[param_id]=discrepancy(data,Qq0, Qq1, observations["rates[1]"][param_id,1],observations["rates[2]"][param_id,1],observations["rates[3]"][param_id,1])
  d_artificial[param_id]=discrepancy(X_tilde,Qq0,Qq1, observations["rates[1]"][param_id,1],observations["rates[2]"][param_id,1],observations["rates[3]"][param_id,1])
}

par(mfrow=c(1,1))
plot(d_data[d_data<d_artificial],d_artificial[d_data<d_artificial], col='red', xlab='D(X|theta)', ylab='D(X_tilde|theta)', xlim=c(-2500,5000), ylim=c(-2500,5000))
points(d_data[d_data>d_artificial],d_artificial[d_data>d_artificial], col='blue')
abline(a=0,b=1)
p_value<-length(d_artificial[d_data<d_artificial])/500
p_value


### ------------------------------------------------------------------------------
### 5. Plot del network
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
matrice_plot_posterior = round(Q[1:num_elem, 1:num_elem]) - 2*diag(dim(Q)[1])


colori = as.numeric(classi$classe[1:num_elem])
colori[which(colori==1)] = 'tomato'
colori[which(colori==2)] = 'royalblue3' #'gray50'
colori[which(colori==3)] = 'springgreen3' #'gold'

g = graph.adjacency(matrice_plot_posterior, mode="undirected", weighted=NULL)
coords = layout_with_fr(g)
# plot the graph
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



