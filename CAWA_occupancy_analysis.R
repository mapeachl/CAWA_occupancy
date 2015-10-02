# Analyze CAWA data using single-visit multi-season Bayesian aproach
# Set working drive
setwd('C:/Users/Michelle/Desktop/WorkSpace/Research/Occupancy')

# Load the data
load("CAwA_sadoti.saved")

# Create dataframe for analysis
# Occupancy data

y <- as.matrix(cbind(cawax$cawa80,cawax$cawa00))

# Covariates -- create and standardize as necessary
# Effort in both time periods
effort <- as.matrix(cbind(cawax$Effort80,cawax$Effort00))
# Replace NA with averages for each atlas period
effort[,1][which(is.na(effort[,1]))] = mean(effort[,1],na.rm=T)
effort[,2][which(is.na(effort[,2]))] = mean(effort[,2],na.rm=T)

effort[,1] <- (effort[,1]-mean(effort[,1])/sd(effort[,1])
effort[,2] <- (effort[,2]-mean(effort[,2])/sd(effort[,2])

# Detection in the earlier atlas
detect80 <- as.matrix(cbind(cawax$Occ80,cawax$Occ80))
detect80[,1] = 0

# Earlier detection in neighboring blocks during the same atlas
detect3 <- as.matrix(cbind(cawax$Priors80,cawax$Priors00))

# Elevation
ele <- cawax$ele
ele <- (ele-mean(ele))/sd(ele)

ele2 <- cawax$ele2
ele2 <- (ele2-mean(ele2))/sd(ele2)

# Forest cover
pfor <- cawax$pfor

# Spatial auto covariate
acov <- cawax$acov

# Neigh
neigh <- cawax$neigh

# Edge residuals
edge_resid <- cawax$edge_resid
edge_resid <- (edge_resid-mean(edge_resid))/sd(edge_resid)

# Create year variables
years <- c(1980:1985)
for (i in 1:length(years)){
  cawax[,dim(cawax)[2]+1] <- ifelse(cawax$X90year80==i,1,0)
}
years <- c(2000:2005)
for (i in 1:length(years)){
  cawax[,dim(cawax)[2]+1] <- ifelse(cawax$X90year00==i,1,0)
}

names(cawax[(dim(cawax)[2]-11):dim(cawax)[2]]) = c("y1980","y1981","y1982","y1983","y1984","y1985","y2000","y2001","y2002","y2003","y2004","y2005")

# Load library
library(R2WinBUGS)

# Specify model in BUGS language
sink("CAWAocc.txt")
cat("
    model {
    
    # Specify priors
    # Beta parameters for the different explanatory variables
    # Occupancy model
    for (k in 1:5){
    BetaO[k] ~ dunif(-20,20)
    }
    
    # Colonization model
    for (k in 1:2){
    BetaC[k] ~ dunif(-20,20)
    }
  
    # Extinction model
    for(k in 1:3){
    BetaE[k] ~ dunif(-20,20)
    }
    
    # Detection model
    gamma~dunif(-20,20)
    
    # Ecological submodel: Define state conditional on parameters
    for (i in 1:nsite){
    
    z[i,2] ~ dbern(psi2[i])
    psi2[i] <- z[i,1]*extant[i] + (1-z[i,1])*col[i]
    z[i,1] ~ dbern(psi1[i])
    psi1[i] <- 1/(1+exp(-logitpsi[i]))
    logitpsi[i] <- BetaO[1] + BetaO[2]*ele[i] + BetaO[3]*ele2[i] + BetaO[4]*pfor[i] + BetaO[5]*acov[i]
    
    extant[i] <- 1-ext[i]

    col[i] <- 1/(1+exp(-logitcol[i]))
    ext[i] <- 1/(1+exp(-logitexp[i]))
    logitcol[i] <- BetaC[1] + BetaC[2]*neigh[i]
    logitexp[i] <- BetaE[1] + BetaE[2]*edge_resid[i] + BetaE[3]*neigh[i]
    
    
    # Observation model
    
    for (j in 1:nyear){
    
    y[i,j] ~ dbern(muy[i,j])
    muy[i,j] <- z[i,j]*pstar[i,j]
    pstar[i,j] <- 1-(1-logitp[i,j])^effort[i,j]
    logitp[i,j] <- 1/(1+exp(-gamma))
    
    } #j
    } #i
    
    }
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nyear = dim(y)[2], 
                 effort=effort,ele=ele,ele2=ele2,pfor=pfor,
                 acov=acov,neigh=neigh,edge_resid=edge_resid)

# Initial values
#zst <- apply(y, c(1, 3), max)       # Observed occurrence as inits for z

inits <- function()
          list(z = matrix(1,dim(y)[1],dim(y)[2]),BetaO=rnorm(5,0,0.001),
                BetaC=rnorm(2,0,0.001),BetaE=rnorm(3,0,0.001),gamma=runif(1,-10,-3))

# Parameters monitored
params <- c("BetaO","BetaC","BetaE","gamma")

# MCMC settings
ni <- 5000
nt <- 4
nb <- 500
nc <- 3

out1 <- bugs(data=win.data, inits=inits, parameters.to.save=params, model.file="CAWAocc.txt",
        n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni,debug=TRUE,DIC=TRUE,working.directory=getwd())

