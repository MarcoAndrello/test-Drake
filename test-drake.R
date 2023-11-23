library(nimble)
library(sf)
library(tidyverse)
library(units)
library(coda)


flexispom <- nimbleCode({
  #~~~~~~~PRIORS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #PSI1 prior
  psi1 ~ dunif(0,1)
  
  # Detection prior
  p_mu ~ dnorm(0,0.001)
  p_sd ~ dunif(0,10)
  p_tau <- pow(p_sd, -2)
  for(t in 1:(nyear.obs)){
    P_t[t] ~ dnorm(p_mu, p_tau)
    logit(p_t[t]) <- P_t[t]
  }
  # Connectivity model priors
  b1_mu ~ dnorm(0, 0.01)
  b1_sd ~ dunif(0,10)
  b1_tau <- pow(b1_sd, -2)
  alpha_mu ~ dnorm(0, 0.01)
  alpha_sd ~ dunif(0, 10)
  alpha_tau <- pow(alpha_sd, -2)
  for(t in 1:(nyear.sim-1)){
    Alpha[t] ~ dnorm(0, alpha_tau)
    alpha[t] <- alpha_mu + c.dyn*Alpha[t]
    sigterm[t] <- 1/(exp(alpha[t])) # sigterm is mean dispersal distance
    
    B1_t[t] ~ dnorm(0, b1_tau)
    b1_t[t] <- exp(b1_mu + c.dyn*B1_t[t])
  }
  
  # Extinction model priors
  # logit(ext) = g0 + g1 * Area
  g0_mu ~ dnorm(0, 0.01)
  g0_sd ~ dunif(0,10)
  g0_tau <- pow(g0_sd, -2)
  g1_mu ~ dnorm(0, 0.01)
  g1_sd ~ dunif(0,10)
  g1_tau <- pow(g0_sd, -2)
  
  #time specific random transition parameters
  for(t in 1:(nyear.sim-1)){
    G0_t[t] ~ dnorm(0, g0_tau)
    G1_t[t] ~ dnorm(0, g1_tau)
    g0_t[t] <- g0_mu + e.dyn*G0_t[t]
    g1_t[t] <- g1_mu + e.dyn*G1_t[t]
  }
  
  #~~~~~~~Likelihood~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  for(i in 1:nsite){ #initial occupancy t0
    z[i,1] ~ dbern(psi1)
  }
  
  for(k in 2:nyear.sim){ #for occupancy t1 and after
    for(i in 1:nsite){
      for(j in 1:nsite){ 
        con[i,j,k-1] <- exp(-sigterm[k-1] * dmat[i,j]) * #kernel
          (1 - equals(i,j)) * #self 
          max(z[j,k-1], struct) *#functional weight
          Area[j] #area weight contrib
      }
      #transition probs
      conx[i,k-1] <- sum(con[i,1:nsite,k-1])
      
      col[i,k-1] <- 1-exp(-b1_t[k-1]*conx[i,k-1]) # akin to Sutherland et al. 2014 to help with model convergence
      logit(ext[i,k-1]) <- g0_mu + g1_mu * Area[i]
      
      #occupancy
      mu.z[i,k-1] <- z[i,k-1] * max(0.001, min((1-ext[i,k-1]), 0.999)) + (1 - z[i,k-1]) * max(0.001, min(col[i,k-1], 0.999)) # min-max trick to prevent calculation issues
      z[i,k] ~ dbern(mu.z[i,k-1]) 
    } 
  }
  #### observation model
  for(i in 1:nsite){
    for (t in 1:nyear.obs){
      mu.p[i, t] <- z[i,t] * p_t[t]
      Y[i, t] ~ dbin(mu.p[i, t], K[i,t])
    }
  }
  #### Derived parameters
  for(t in 1:nyear.sim){
    m.occ[t] <- sum(z[1:nsite,t]) 
  } 
})


## Make up data
nsite <- 98; nyear.obs <- nyear.sim <- 17
# Distance matrix
sp <- st_read("sampling_points.shp")
dmat <- st_distance(sp)
dmat %>% drop_units -> dmat
set.seed(123)
#Area
Area = runif(nsite,50,3000)
# frequencies of observation
Y <- matrix(sample(c(0,0,0,1,2,3,4),nsite*nyear.obs,replace=T),
            nrow=nsite, ncol=nyear.obs)
K <- matrix(4,nrow=nsite, ncol=nyear.obs)
z <- matrix(as.numeric(as.logical(Y)),
       nrow=nsite, ncol=nyear.obs)

# load("Data.RData")

data <- list(Area=Area,
             Y=Y,
             K=K,
             dmat=dmat,
             z=z)
#1. struct. connectivity (nwork position only) + static effect (beta_t = beta) (model UI)
#2. struct. connectivity (nwork position only) + dynamic effect (beta_t) (model UV)
#3. funct. connectivity (z-weighted) + static effect (beta_t = beta) (model DI)
#4. funct. connectivity (z-weighted) + dynamic effect (beta_t) (model DV)
#1 and 2 are *non*-demographic or demographically naive #3 and 4 are demographic connectivity
#model 1
sta.consts.struct <- list(nyear.obs=nyear.obs, nyear.sim=nyear.sim, nsite=nsite, c.dyn=0, #0=invariant, 1=time-varying
                          e.dyn=0, #0=invariant, 1=time-varying
                          struct=1) #1=structural, 0=functional
#model 2
dyn.consts.struct <- list(nyear.obs=nyear.obs, nyear.sim=nyear.sim, nsite=nsite, c.dyn=1, e.dyn=0, struct=1)#1=structural, 0=functional
#model 3
sta.consts <- list(nyear.obs=nyear.obs, nyear.sim=nyear.sim, nsite=nsite, c.dyn=0, e.dyn=0, struct=0) #1=structural, 0=functional
#model 4
dyn.consts <- list(nyear.obs=nyear.obs, nyear.sim=nyear.sim, nsite=nsite, c.dyn=1, e.dyn=0, struct=0) #1=structural, 0=functional
# Parameters to track
params <- c("alpha","b1_t","m.occ","sigterm", "Alpha", "alpha_mu", "b1_mu", "B1_t")
inits <- function(){
  list( psi1=runif(1,0.1,0.9),
        p_mu=rnorm(1,0,0.1),
        p_sd=runif(1,0.1,1),
        P_t=rnorm(nyear.obs,0,0.1),
        alpha_mu=rnorm(1,0,0.1),
        alpha_sd=runif(1,0.1,1),
        b1_mu=rnorm(1,0,0.1),
        b1_sd=runif(1,0.1,1),
        B1_t=rnorm(nyear.sim-1,0,0.1),
        g0_mu=runif(1,-1,1),
        g1_mu=rnorm(1,-1,0.1),
        G0_t=rnorm(nyear.sim-1,0,0.1),
        G1_t=rnorm(nyear.sim-1,0,0.1))
  }
set.seed(123)
inits_list <- inits()


mDV <- nimbleModel(code=flexispom,
                   constants=dyn.consts,
                   data=data,
                   inits=inits_list,
                   name="DV")

MCMCconf <- configureMCMC(mDV, monitors = params, enableWAIC = T)
mDVmcmc <- buildMCMC(MCMCconf,niter = 80000, nburnin = 30000,summary = TRUE, WAIC =FALSE, check= TRUE, samples = TRUE, samplesAsCodaMCMC=TRUE)
cmDV <- compileNimble(mDV)
cmDVmcmc <- compileNimble(mDVmcmc , project = mDV)

samples <- runMCMC(cmDVmcmc, niter = 80000, nburnin = 30000,summary = TRUE, WAIC =FALSE, samples = TRUE, samplesAsCodaMCMC=TRUE)
save(samples, file="samples_DV.RData")
colnames(samples$samples)
barplot(apply(samples$samples,2,sd),las=2)
plot(samples$samples[, "Alpha[2]"], type = "l")
plot(samples$samples[, "b1_t[1]"], type = "l")
plot(samples$samples[, "m.occ[1]"], type = "l") # non cambia: errore o dovuto ai dati troppo informativi?


#######
save(mp_DV, file=paste("mp_DV",format(Sys.time(), "%Y%m%d"), ".RData", sep=""))
mp_UV <- nimbleMCMC(code=flexispomv, constants=dyn.consts.struct, data=data, inits=inits,
                    monitors = params, nchains=3, niter = 80000, nburnin = 30000,
                    thin = 1, summary = TRUE, WAIC = FALSE, check= TRUE, samples = TRUE, samplesAsCodaMCMC=TRUE)
save(mp_UV, file=paste("mp_UV",format(Sys.time(), "%Y%m%d"), ".RData", sep=""))

mp_DI <- nimbleMCMC(code=flexispom, constants=sta.consts, data=data, inits=inits,
                    monitors = params, nchains=3, niter = 80000, nburnin = 30000,
                    thin = 1, summary = TRUE, WAIC = FALSE, check= TRUE, samples = TRUE, samplesAsCodaMCMC=TRUE)
save(mp_DI, file=paste("mp_DI",format(Sys.time(), "%Y%m%d"), ".RData", sep=""))

mp_UI <- nimbleMCMC(code=flexispom, constants=sta.consts.struct, data=data, inits=inits,
                    monitors = params, nchains=3, niter = 80000, nburnin = 30000,
                    thin = 1, summary = TRUE, WAIC = FALSE, check= TRUE, samples = TRUE, samplesAsCodaMCMC=TRUE)
save(mp_UI, file=paste("mp_UI",format(Sys.time(), "%Y%m%d"), ".RData", sep=""))
