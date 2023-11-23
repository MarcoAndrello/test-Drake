# Test Eryngium
rm(list=ls())

library(tidyverse)
library(nimble)
library(coda)

load("data.surv.RData")
load("xvar_1.RData")

left_join(data.surv, xvar_1, by=c("Site", "Year")) %>% filter(Year %in% c(2014:2021)) -> data
data <- na.omit(data)
# Scale predictors
data[,c("ddays", "num_s_Tmin_15")] <- scale(data[,c("ddays", "num_s_Tmin_15")])
summary(data)

data$Site_ids <- as.numeric(data$Site)
data$State_ids <- as.numeric(data$State)
data$ID_ids <- as.numeric(data$ID)

EAconstants <- list(num_sites = 7,
                    num_states = 3,
                    # num_IDs = 3750,
                    num_obs = 6123,
                    Site = data$Site_ids,
                    State = data$State_ids,
                    # ids = data$ID_ids,
                    ddays = data$ddays,
                    num_s_Tmin_15 = data$num_s_Tmin_15
                    )

EAdata <- list(fateSurv = as.numeric(as.character(data$fateSurv))) # to get back the 0/1

EAcode <- nimbleCode({
  # logit link and Bernoulli data probabilities
  for(i in 1:num_obs) {
    logit(survival_probability[i]) <-
      Site_int[ Site[i] ] +
      State_int[ State[i] ] +
      ddays_coef * ddays[i] +
      num_s_Tmin_15_coef * num_s_Tmin_15[i] #+
      # ID_effect[ ids[i] ]
    fateSurv[i] ~ dbern(survival_probability[i])
  }
  #Priors
  for (i in 1 : num_sites){
    Site_int[i] ~ dnorm(0, sd = 1000)
  }
  for (i in 1 : num_states){
    State_int[i] ~ dnorm(0, sd = 1000)
  }
  # Priors for farm random effects and their standard deviation.
  # ID_sd ~ dunif(0, 20)
  # for(i in 1:num_IDs) {
  #   ID_effect[i] ~ dnorm(0, sd = ID_sd)
  # }
  ddays_coef ~ dnorm(0, sd = 1000)
  num_s_Tmin_15_coef ~ dnorm(0, sd = 1000)
})


EAinitialize <- function() {
  list(Site_int = rnorm(7, 0, 1),
       State_int = rnorm(3, 0, 1),
       ddays_coef = 0,
       num_s_Tmin_15_coef = 0)#,
       # ID_sd = 1,
       # ID_effect = rnorm(3750, 0, 1))
}

set.seed(123)
EAinits <- EAinitialize()

EAmodel <- nimbleModel(EAcode,
                       constants = EAconstants,
                       data = EAdata,
                       inits = EAinits)
# Build the MCMC
EAmcmc <- buildMCMC(EAmodel, enableWAIC = TRUE)

# Compile the model and MCMC
cEAmodel <- compileNimble(EAmodel)
cEAmcmc <- compileNimble(EAmcmc, project = EAmodel)

# Run the MCMC
time_baseline <- system.time(EAresults <- runMCMC(cEAmcmc, niter=50000, nburnin=25000, thin=100, WAIC=TRUE))
cat("Sampling time: ", time_baseline[3], "seconds.\n")

samples <- EAresults$samples
WAIC <- EAresults$WAIC

plot(as.mcmc(samples))
