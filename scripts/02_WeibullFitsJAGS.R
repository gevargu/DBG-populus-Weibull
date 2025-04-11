#' # Model fit: Xylem Vulnerability Curves
#' 
#' ### Author: German Vargas G.
#' 
# Load Libraries ----
library(tidyverse)
library(udunits2)
library(rjags)
library(runjags)
load.module('dic')
library(mcmcplots)
library(bayestestR)
#devtools::install_github("fellmk/PostJAGS/postjags")
library(postjags)

# Load custom functions ----
set.seed(123)
n_cores = parallel::detectCores()
'%ni%' <- Negate('%in%')
source(file = "scripts/00_GammaModeFunction.R")

# Load Data ----
plc_curves <- read.csv("data_clean/CavitronDataClean.csv",row.names = "X")# load clean cavitron data
plc_curves <- plc_curves[,c("Sample_ref_1","Sample_ref_2","Raw_conductance_kg_Mpa_s","Conductivity_SI_corrT","Pressure_Mpa","kmax","Ksmax","PLCclean")]# select columns of interest

# Extract Tree.No and sample Replicate
plc_curves$Tree.No <- as.numeric(x = gsub(pattern = "\\..*",replacement = "",x = plc_curves$Sample_ref_1))# extract tree number
plc_curves$Rep.No <- as.numeric(x = ifelse(test = grepl(pattern = "\\.", x = plc_curves$Sample_ref_1),
                                           yes = sub("^[^.]*\\.([^.]*)", "\\1", plc_curves$Sample_ref_1),
                                           no = "1"))# extract replicate number
head(plc_curves)

# merge data with genotype information
sample_id <- read.csv("data_clean/dbg_cottonwood_plantID.csv")# load sample id information
plc_curves <- merge(x = plc_curves, y = sample_id, by = "Tree.No", all.x = TRUE)# merge sample id information with plc data
plc_curves$sample.id <- paste(plc_curves$Popn.Geno,plc_curves$Tree.No,plc_curves$Rep,sep = "-")# create sample id
plc_curves <- plc_curves[,c("Population","Genotype","Popn.Geno","Tree.No","Rep.No","sample.id","Elevation.m","Pressure_Mpa","PLCclean")]# select columns of interest
plc_curves <- plc_curves %>% arrange(Population,Genotype,Tree.No,Rep.No,desc(Pressure_Mpa))# arrange data
plc_curves2 <- plc_curves %>% group_by(Population,Genotype,Popn.Geno,Tree.No,Rep.No,sample.id,Elevation.m,Pressure_Mpa) %>%
  summarise(PLCclean = mean(PLCclean,na.rm = T)) %>% # summarise data
  arrange(Population,Genotype,Tree.No,Rep.No,desc(Pressure_Mpa))# arrange data
row.names(plc_curves2) <- NULL
head(plc_curves2)

# Quick plot of the data
ggplot(plc_curves, aes(x=Pressure_Mpa, y=PLCclean, col=factor(sample.id)))+ geom_point() + 
  geom_smooth(se=F)  + facet_wrap(~sample.id) +theme_bw() + theme(legend.position = "none") 
ggplot(plc_curves2, aes(x=Pressure_Mpa, y=PLCclean, col=factor(sample.id)))+ geom_point() + 
  geom_smooth(se=F)  + facet_wrap(~sample.id) +theme_bw() + theme(legend.position = "none") 

# JAGS model ----
# P12 and P50 distributions ----------------------------------------------------
# Bayesian hierarchical model --------------------------------------------------
xylemvc.hmodel <- "model {
  # Define “expected” process according to the reparameterized - line 2
  # Weibull function: - line 3
  for(i in 1:N){
    # Likelihood for observed data: - line 5
    KSrel[i] ~ dnorm(muK[i],nuInv)
    # Model for replicated data: - line 7
    KSrep[i] ~ dnorm(muK[i],nuInv)
    
    # Expected PLC - line 10
    # muK[i,m] <- ksat[i]*p4[i,m], this is an extra step that I don't understand why?
    muK[i] <- pow((1-X/100),p3[i])
    p3[i] <- pow(p1[i],p2[i])
    p2[i] <- (Px[i]*Sx[i])/V
    p1[i] <- MPa[i]/Px[i]
    
    # Individual level parameters - line 17
    Px[i] <- Tree.Px[treeID[i]]
    Sx[i] <- Tree.Sx[treeID[i]]
  }

  # Parameter models: - line 22
  # Tree level random effects: - line 23
  for(i in 1:NTree){
    Tree.Px[i] ~ dnorm(mu.Px[i],tau.Px.Tree)T(0.5,)
    mu.Px[i] <- Px.0 + inprod(beta.Px[],Xmat[i,])
    Tree.Sx[i] ~ dnorm(mu.Sx[i],tau.Sx.Tree)T(0,)
    mu.Sx[i] <- Sx.0 + inprod(beta.Sx[],Xmat[i,])
    # Note: used indicator function, T(0,), to ensure that Tree.Px, and Tree.Sx are positive - line 29
  }
  
  # Constant for Weibull model: - line 32
  X <- 50
  V <- (X-100)*log((1-X/100))
  
  # hyper-priors on population effects - line 36
  for (i in 1:NTRT) {# level 2 (ex: populations) - line 37
    beta.Px[i] ~ dnorm(0.0,tau.Px.TRT)
    beta.Sx[i] ~ dnorm(0.0,5)# this is a compromise to get good parameter estimation given the data. We are assuming the shape should not change much
  }
  
  # Priors for standard deviation parameters, and compute associated variances
  # and precision estimates: - line 43
  nuInv <- 1/varPr
  varPr ~ dunif(0,1)
  Px.0 ~ dnorm(0,0.001)T(0.5,)
  Sx.0 ~ dnorm(0,0.001)T(20,)
  tau.Px.Tree <- 1/pow(tau.Px,2)
  tau.Px ~ dunif(0,2)
  tau.Sx.Tree <- 1/pow(tau.Sx,2)
  tau.Sx ~ dunif(0,2)
  tau.Px.TRT <- 1/(sigma.Px.TRT)^2
  sigma.Px.TRT ~ dgamma(TRT.GammaShRa.Px[1],TRT.GammaShRa.Px[2])# informed by literature values
  #tau.Sx.TRT <- 1/(sigma.Sx.TRT)^2
  #sigma.Sx.TRT ~ dgamma( 1.105125, 0.1051249 )# sigma.Sx.TRT follows a gamma dist based on mode 1 and sd 10
  # Note: used indicator function, T(0,), to ensure that mu.Px, and mu.Sx are positive - line 56
  
  # Derived quantities - line 58
  # Global level means - line 59
  P50.global <- Px.0
  S50.global <- Sx.0
  P12.global <- P50.global*(pow((log(.88)/log(.5)),(V/(P50.global*S50.global))))
  
  # Treatment level means: - line 64
  # P50 - 65
  for ( i in 1:NTRT) { 
    P50.TRT[i] <- Px.0 + beta.Px[i] # Treatment level P50  - line 61
    S50.TRT[i] <- Sx.0 + beta.Sx[i] # Treatment level shape - line 62
  }
  # P12 - line 70
  for (i in 1:NTRT) {
    P12.TRT[i] <- P50.TRT[i]*(pow((log(.88)/log(.5)),(V/(P50.TRT[i]*S50.TRT[i]))))
  }
}"

# Prepare data ----------------------------------------------------------------------
####### a. Get relative conductivity ---------------------------------------------------------
plc_curves2$KSrel <- 1-plc_curves2$PLCclean/100
ggplot(plc_curves2, aes(x = Pressure_Mpa, y = KSrel, color = Population,group = sample.id)) + 
  facet_wrap(~Population) +
  geom_point(color="#000000") + 
  geom_line(stat = "smooth",alpha=0.5,lwd=1)+
  ylim(-0.15,1.15) +
  theme_bw()
row.names(plc_curves2) <- NULL
plc_curves2 <- droplevels(plc_curves2)
head(plc_curves2)

####### b. data list for JAGS --------------------------------------------------
# xylem vc data
# relative conductivity and water potentials
KSrel <- plc_curves2$KSrel
MPa <- plc_curves2$Pressure_Mpa*-1
N <- length(MPa)

# level 1 : sample-level id
plc_curves2$ID <- factor(plc_curves2$sample.id,levels = unique(plc_curves2$sample.id))
ID <- as.numeric(plc_curves2$ID)
NTree <- length(unique(ID))

# level 2: assemble model matrix for elevation band effects
plc_curves2$PopulationF <- factor(plc_curves2$Population,levels = unique(plc_curves2$Population),
                                 labels = unique(plc_curves2$Population))
Xmat <- model.matrix(~PopulationF-1,plyr::ddply(plc_curves2[,c("Pressure_Mpa","ID","PopulationF")],~ID,plyr::catcolwise(unique)))
NTRT <- ncol(x = Xmat)

# hyper prior deflection, used data here from the xylem FDB and the genus populus
read.csv("data_clean/XFT_populus.csv") -> populus.p50
GammaShRa.Px <- unlist(gammaShRaFromModeSD(mode = abs(median(x = populus.p50$P50..MPa.,na.rm = T)),
                                           sd = 2*sd(x = populus.p50$P50..MPa.,na.rm = T)))
# assemble data list
dataList <- list(KSrel = KSrel,
                 MPa = MPa, 
                 N = N, 
                 treeID = ID, 
                 NTree = NTree , 
                 Xmat = Xmat, 
                 NTRT = NTRT, 
                 TRT.GammaShRa.Px = GammaShRa.Px)

####### c. JAGS model specification --------------------------------------------------
# parameters to monitor
parameters <- c("P50.global","S50.global","P12.global","P50.TRT","S50.TRT","P12.TRT","varPr","tau.Px","tau.Sx","Px.0","Sx.0","beta.Px","beta.Sx","Tree.Px","Tree.Sx","KSrep")

# initial values
# check the distribution of P50 for the genus Populus
hist((populus.p50[!is.na(populus.p50$P50..MPa.),"P50..MPa."]))

inits <- function() {
  list(
    varPr = runif(1, 0.001, 1), 
    tau.Px = runif(1, 0.001, 2),
    tau.Sx = runif(1, 0.001, 2),
    sigma.Px.TRT = runif(1, 0.1, 2),
    Px.0 = runif(1, 1, 4), 
    Sx.0 = runif(1, 20, 200), 
    beta.Px = runif(NTRT, -1, 1),
    beta.Sx = runif(NTRT, -1, 1),
    Tree.Px = runif(NTree, 1, 4), 
    Tree.Sx = runif(NTree, 10, 200) 
  )
}
inits_list <- list(inits(), inits(), inits(),inits())

# Or, load previous saved state
load("output/coda/inits.Rdata")
saved_state[[2]]

# chains specification. Note: modify until ESS > 10000
thinSteps <- 150
nIter <- 750000
nChains <- 4
burnInSteps <- floor(nIter/nChains)
adaptSteps <- round(burnInSteps*0.25)

# Initialize model
jm <- jags.model(textConnection(xylemvc.hmodel), data = dataList, inits = inits_list, n.chains = nChains,n.adapt = adaptSteps)
#jm <- jags.model(textConnection(xylemvc.hmodel), data = dataList, inits = saved_state[[2]], n.chains = nChains,n.adapt = adaptSteps)

# Monitor coda samples
jm_coda <- coda.samples(model = jm, variable.names = parameters, n.iter = nIter, thin = thinSteps)

# Save the JAGS output to an .RData file
save(jm_coda, file = "output/coda/jm_coda.Rdata")
load("output/coda/jm_coda.Rdata")

# Visualize chains
mcmcplot(jm_coda,parms = c("P50.global","P50.TRT","S50.global","S50.TRT","P12.global","P12.TRT"))
mcmcplot(jm_coda,parms = c("tau.Px","tau.Sx","R2"))
mcmcplot(jm_coda,parms = c("Px.0","Sx.0","beta.Px","beta.Sx"))
mcmcplot(jm_coda,parms = c("Tree.Px","Tree.Sx"))

# Save state
# newinits <- initfind(jm_coda, OpenBUGS = FALSE)
# newinits[[1]]
# saved_state <- removevars(initsin = newinits, 
#                           variables = c(1:5,7,8))
# save(saved_state, file = "output/coda/inits.Rdata")

# Check convergence
gel <- gelman.diag(jm_coda, multivariate = FALSE)
gel$psrf[gel$psrf[,1] > 1.01,]# still many variables are not converging. Need to increase the number of iterations
ess <- effectiveSize(jm_coda)# check ESS
ess[ess < 1000]
