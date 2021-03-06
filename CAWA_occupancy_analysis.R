# Analyze CAWA data using single-visit multi-season Bayesian aproach
# Set working drive
setwd('C:/Users/Michelle/Desktop/WorkSpace/Research/Occupancy')

# Load the data
load("CAwA_sadoti.saved")

# Create dataframe for analysis
# Occupancy data

y <- data.frame(cawax$cawa80,cawax$cawa00)

# Covariates
# Effort in both time periods
effort <- data.frame(cawax$Effort80,cawax$Effort00)

# Detection in the earlier atlas
detect80 <- data.frame(cawax$Occ80,cawax$Occ80)
detect80[,1] = 0

# Earlier detection in neighboring blocks during the same atlas
detect3 <- data.frame(cawax$Priors80,cawax$Priors00)

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
library(rjags)

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
    for(k in 1:16){
    BetaD[k] ~ dunif(-20,20)
    }
    
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
    logitp[i,j] <- 1/(1+exp(-a[i,j]))
    a[i,j] <- BetaD[1] + BetaD[2]*pfor[i] + BetaD[3]*detect80[i,j] + BetaD[4]*detect3[i,j] + 
              BetaD[5]*y1980[i] + BetaD[6]*y1981[i] + BetaD[7]*y1982[i] + BetaD[8]*y1983[i] +
              BetaD[9]*y1984[i] + BetaD[10]*y1985[i] + BetaD[11]*y2000[i] + BetaD[12]*y2001[i] +
              BetaD[13]*y2002[i] + BetaD[14]*y2003[i] + BetaD[15]*y2004[i] + BetaD[16]*y2005[i]
    
    } #j
    } #i
    
    # Derived parameters: Sample and population occupancy, growth rate and turnover
    #don't run this for now just makes everything take longer 
    #meanpsi[1] <- mean(psi1[])
    #meanpsi[2] <- mean(psi2[])
    #n.occ[1]<-sum(z[1:nsite,1])
    #pocc[1] <- n.occ[1]/nsite
    #n.occ[2]<-sum(z[1:nsite,2])
    #pocc[2] <- n.occ[2]/nsite
    
    #for(i in 1:nsite){
    #colyn[i] <- step(z[i,2]-z[i,1]-1)
    #extyn[i] <- step(z[i,1]-z[i,2]-1)
    #}
    #scol <- sum(colyn[])
    #sext <- sum(extyn[])
    # Step function and equals funcion are work arounds for if/then  
    # # growthr <- psi1/psi2
    # # turnover <- (1 - psi1) * col/psi2
    
    }
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nyear = dim(y)[2], 
                 effort=effort,ele=cawax$ele,ele2=cawax$ele2,pfor=cawax$pfor,
                 acov=cawax$acov,neigh=cawax$neigh,edge_resid=cawax$edge_resid,
                 detect80=detect80,detect3=detect3,y1980=cawax$y1980,y1981=cawax$y1981,
                 y1982=cawax$y1982,y1983=cawax$y1983,y1984=cawax$y1984,y1985=cawax$y1985,
                 y2000=cawax$y2000,y2001=cawax$y2001,y2002=cawax$y2002,y2003=cawax$y2003,
                 y2004=cawax$y2004,y2005=cawax$y2005)

# Initial values
#zst <- apply(y, c(1, 3), max)       # Observed occurrence as inits for z

inits <- function(){list(z = matrix(1,dim(y)[1],dim(y)[2]),BetaO=rnorm(5,0,0.001),
                         BetaC=rnorm(2,0,0.001),BetaE=rnorm(3,0,0.001),BetaD=rnorm(16,0,0.001))}

# Parameters monitored
params <- c("BetaO","BetaC","BetaE","BetaD")

# MCMC settings
ni <- 5000
nt <- 4
nb <- 500
nc <- 3

test1.glmm <- jags.model("CAWAocc.txt", win.data, inits, n.chain=nc, n.adapt=100)

