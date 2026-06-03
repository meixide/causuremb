set.seed('1023')

# Initial setup
nlambdas = 5
lambdas = seq(0.1, 1, length.out = nlambdas)
source('simulate_tail.R')
source('count_me.R')

# Gaussian Kernel Function
gaussian_kernel <- function(x, y, sigma2) {
  if (length(x) > 1 || length(y) > 1) return('Only scalars are accepted as inputs!')
  return(exp(-((x - y)^2) / sigma2))
}


# Auxiliary function to perform goodness-of-fit
gi <- function(x, ui, deltai) {
  g = numeric(length(x))
  for (k in 1:length(x)) {
    xk = x[k]
    if (xk <= ui) {
      g[k] = 0
    } else {
      g[k] = deltai + (1 - deltai) * (xk - ui) / (1 - ui)
    }
  }
  return(g)
}

# Simulation Parameters
n = 800
m = 800
B = 100
lgrid = 200

# Let us simulate an scenario with confounding
beta0=c(1,0)
beta1=beta0
shift=2
shift_time=0
delta0=1.5
delta1=1.5


# Simulate data
sim = simulate_tail(n, m, beta0, delta0, shift, shift_time, delta1, beta1)
tmax = max(sim$Tstar)
tmin = 0

library(Matching)

tall=sim$Tstar
zall=sim$z
deltall=sim$delta
xall=sim$x

matched_data <- Match(Y = NULL, Tr = zall, X = xall,caliper=0.2,exact=F,ties=FALSE)
summary(matched_data)

lgrid=100
xproba=seq(0,1,length.out=100)

# Cross-validation loop
K = 5
folds = sample(cut(seq(1, length(tall)), breaks = K, labels = FALSE), replace = FALSE)
average_errors = numeric(nlambdas)

for (f in 1:K) {
  test.ind <- which(folds == f)
  matc_subset <- which(matched_data$index.control %in% test.ind)
  print(length(matc_subset))
  placebot=tall[matched_data$index.control[matc_subset]]
  deltam=deltall[matched_data$index.control[matc_subset]]
  ntest=length(placebot)
  t=tall[-test.ind]
  z=zall[-test.ind]
  x=xall[-test.ind,]
  delta=deltall[-test.ind]
  sigma2kernel=median((as.numeric(dist(t[z==0 & delta==1 ]))^2))
  
  
  # Lambda tuning
  for (j in 1:nlambdas) {
    fd = numeric(ntest)
    weights = count_me(t, delta, z, x, lgrid, tmax, lambdas[j], sigma2kernel, tmin)[[2]]

    for(i in 1:ntest) {
      fd[i]=sum(weights[which(t[which(z==0)]<placebot[i])])
    }
    
    
    m=ntest
    
    ghat = Reduce('+', lapply(1:m, function(y) gi(xproba, fd[y], deltam[y]) / m))
    average_errors[j] = mean((ghat - xproba)^2)
  }
  
  average_errors = average_errors + average_errors / K

}

# Plot results
plot(1:nlambdas,average_errors)


t=sim$Tstar
z=sim$z
delta=sim$delta
x=sim$x


km0=survfit(Surv(t[z==0],delta[z==0])~1)
km1=survfit(Surv(t[z==1],delta[z==1])~1)

weights=count_me(t,delta,z,x,lgrid,tmax,lambdas[which.min(average_errors)]
                 ,sigma2kernel,tmin)[[2]]
pt=seq(tmin,tmax,length.out=lgrid)
surv01=numeric(lgrid)
for(i in 1:lgrid) {
  surv01[i]=1-sum(weights[which(t[which(z==0)]<pt[i])])
}

plot(km0,ylim=c(0,1),conf.int=F,ylab='Survival probability',xlab='Log-time')
lines(km1,col='blue',conf.int=F)
lines(pt,surv01,col='red')



    







