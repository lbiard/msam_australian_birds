library(stringr)
library(R2jags)
library(tictoc)

## Need to load RData file containing the input data (it already contains the output as well).

###############
#### Model ####
###############

model_string <-"
    model {
  
  #PRIORS
  mu.b0 ~ dnorm(0,0.1)
  mu.a0 ~ dnorm(0,0.1)
  mu.phi ~ dnorm(0,0.1)
  tau.b0 ~ dgamma(0.1,0.1)
  tau.a0 ~ dgamma(0.1,0.1)
  tau.phi ~ dgamma(0.1,0.1)
  tau.phiSpp ~ dgamma(1,5)
  
  
  for (i in 1:(n.species)) {
    #Create Priors for species i from community level prior distributions
    for (j in 1:n.sites) {
    b0[i,j] ~ dnorm(mu.b0, tau.b0)
    spp.phi[j,i] ~ dnorm(mu.phi, tau.phi)
    }
    a0[i] ~ dnorm(mu.a0, tau.a0)

    
    #Season-specific autologistic term
    for(l in 1:n.seasons) {
        for (j in 1:n.sites) {
      phi[j, l, i] ~ dnorm(spp.phi[j,i], tau.phiSpp)
    }
    }
    
    #Estimation of N matrix (true abundance for species i at point j in season t = 1)
    for (j in 1:n.sites) {
      log.lambda[j, f1[j], i] <- b0[i,j]
      lambda[j, f1[j], i] <- exp(log.lambda[j, f1[j], i])
      N[j, f1[j], i] ~ dpois(lambda[j, f1[j], i])
      
      #Estimate detection for species i at point j during sampling period k
      for (k in 1:n.visits) {
        logit.p[j, k, f1[j], i] <- a0[i]
        p[j, k, f1[j], i] <- 1 / (1 + exp(-logit.p[j, k, f1[j], i]))
        y[j, k, f1[j], i] ~ dbin(p[j, k, f1[j], i], N[j, f1[j], i])
      }
      
      ## Time loop for season t > 1
      for (t in (f1[j]+1):(f2[j])) {
        log.lambda[j, t, i] <- b0[i,j] + phi[j, t-1, i] * N[j, t-1, i]
        lambda[j, t, i] <- exp(log.lambda[j, t, i])
        N[j, t, i] ~ dpois(lambda[j, t, i])
        
        #Estimate detection for species i at point j during sampling period k and season t
        for (k in 1:n.visits) {
          logit.p[j, k, t, i] <- a0[i]
          p[j, k, t, i] <- 1 / (1 + exp(-logit.p[j, k, t, i]))
          y[j, k, t, i] ~ dbin(p[j, k, t, i], N[j, t, i])
        }
      }		
    }
  }
    }
    "

jags.data <- list(
  y = y,
  n.sites = dim(y)[1],
  n.visits = dim(y)[2],
  n.seasons = dim(y)[3],
  n.species = dim(y)[4],
  f1 = f1,
  f2 = f2)

# Setting initial values for latent variable N
N.naive <- apply(y, MARGIN = c(1,3,4), max, na.rm=T) 
N.naive[N.naive == -Inf] <- NA
N.miss <- apply(N.naive, c(1, 3), max, na.rm = T)
N.naive <- array(dim = dim(N.naive))
for(i in 1:dim(N.naive)[3]) {
  for(j in 1:dim(N.naive)[1]) {
    N.naive[j, f1[j]:f2[j] , i] <- N.miss[j, i]
  }
}
inits <- list(N = N.naive)

## MCMC parameters ##
nc <- 3
n.burn <- 75000
n.iter <- 150000
thin <- 375


inits <- function(){list(N = N.naive)}  

# Parameters monitored
# run first saving underlying parameters to ensure convergence, and then another time saving only "N" to get abundance

#parameters <- c("N")   
parameters <- c("a0", "b0", "spp.phi", "phi", "tau.phiSpp", "mu.phi", "tau.phi",
                "mu.b0", "tau.b0", "mu.a0", "tau.a0")

# Call JAGS 
tic()
set.seed(26)
model.Nmat <- jags(jags.data, inits, parameters, model.file = textConnection(model_string), n.chains = nc,
                   n.thin = thin, n.iter = n.iter, n.burnin = n.burn, working.directory = getwd())
toc()
print(model.Nmat, digits = 3, intervals=c(0.025, 0.975))

# r-hat values for all saved parameters, to check convergence of underlying parameters
plot(density(model.Nmat$BUGSoutput$summary[,8]), xlim=c(1,1.5))



#########################
#### Post-processing ####
#########################

## get N.mat array with estimated abundance from posterior draws (when saving parameters "N")

N.1 <- as.mcmc(model.Nmat)
N.1 <- as.matrix(N.1)
N.2 <- as.data.frame(N.1)
N.2 <- N.2[,-1] #remove deviance parameter

store_key_array <- matrix(NA, nrow = length(colnames(N.2)), ncol = 3) # columns are site, season, species
temporary_storage <- NA

# extract site, season, species from the name of each saved parameter
for (i in 1:length(colnames(N.2))) {
  temporary_storage <- str_extract_all(colnames(N.2)[i], "[0-9]+")
  for (j in 1:3) {
    store_key_array[i,j] <- temporary_storage[[1]][j]
  }
}

# create N.mat array to store estimated abundance
N.mat <- array(dim = c(jags.data$n.sites, jags.data$n.seasons, jags.data$n.species, dim(N.1)[1]))
# Fill N.mat with estimated abundance
for (i in 1:dim(store_key_array)[1]) {
  N.mat[as.numeric(store_key_array[i,1]),
        as.numeric(store_key_array[i,2]),
        as.numeric(store_key_array[i,3]), ] <- N.2[,i]
}

# cap N.mat at 300
N.mat[N.mat>300] <- 300  # cap at 300

