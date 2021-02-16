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
  file    = "Stan/Poisson-ER.stan",   # Stan program
  data    = X_data,         # named list of data
  chains  = 1,              # number of Markov chains
  warmup  = 1000,           # number of warmup iterations per chain
  iter    = 2000,           # total number of iterations per chains
  cores   = 1,              # number of cores (could use one per chain)
  refresh = 1               # progress shown
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

par(mfrow=c(1,2))
plot(data, asp=T, col=cols2, main='Original Dataset', border=NA)
plot(Q, asp=T, col=cols2, main='Posterior', border=NA)

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

###########################################
#   Histograms with fitted distributions  #
###########################################

dataframeQ=as.data.frame(cbind(as.vector(round(Q)),as.vector(data)))
A0=dataframeQ[dataframeQ$V1==0,]
A1=dataframeQ[dataframeQ$V1==1,]
A11=A1[A1$V2!=0,]
par(mfrow=c(1,2))
hist(A0$V2, freq=FALSE, xlab='X|A=0', main="Histogram of Binomial(p_0)")
x=seq(from=0, to=15, by=1)
y=dpois(x,rates_1)
lines(x,y, col='red')
hist(A11$V2, freq=FALSE, xlab='X|A=1', breaks=20, ylim=c(0,0.06), main="Histogram of Binomial(p_1)")
x=seq(from=10, to=270, by=1)
y=dpois(x,rates_2)
lines(x,y, col='red')

### ------------------------------------------------------------------------------
### 4. Posterior Predictive Assessment with Discrepancy
### ------------------------------------------------------------------------------


discrepancy<-function(X, Q, lambda0, lambda1){
  #creo la matrice X_tilde mat come dato artificiale
  X_star=matrix(nrow=n,ncol=n)
  for(i in 1: n){
    for(j in i:n)
      X_star[i,j]=Q[i,j]*lambda0+(1-Q[i,j])*lambda1
  }
  
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
  #Creo la matrice Q e A
  for(i in 1:n){
    for(j in i:n){
      var = paste('Q[',i,',',j,']', sep='')
      Qq[i,j] = as.matrix(observations[var])[param_id,1]
      A[i,j]=as.numeric(runif(1,0,1)<Qq[i,j])
      A[j,i]=A[i,j]
    }
  }
  #Creo la matrice X_tilde di osservazioni artificial
  X_tilde_0=matrix(rpois(n*n, observations["rates[1]"][param_id,1]),nrow=n,ncol=n)
  X_tilde_1=matrix(rpois(n*n, observations["rates[2]"][param_id,1]),nrow=n,ncol=n)
  X_tilde=(1-A)*X_tilde_0+A*X_tilde_1
  d_data[param_id]=discrepancy(data,Qq,observations["rates[1]"][param_id,1],observations["rates[2]"][param_id,1])
  d_artificial[param_id]=discrepancy(X_tilde,Qq,observations["rates[1]"][param_id,1],observations["rates[2]"][param_id,1])
}

par(mfrow=c(1,1))
plot(d_data[d_data<d_artificial],d_artificial[d_data<d_artificial], col='red', xlab='D(X|theta)', ylab='D(X_tilde|theta)', xlim=c(30000,50000), ylim=c(30000,50000))
points(d_data[d_data>d_artificial],d_artificial[d_data>d_artificial], col='blue')
abline(a=0,b=1)
post_pred_pval<-length(d_artificial[d_data<d_artificial])/500
post_pred_pval

### ------------------------------------------------------------------------------
### 5. Plot del Network
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

# Plot the graph
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



