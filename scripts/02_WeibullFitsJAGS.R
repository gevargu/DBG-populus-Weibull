#' # Model parameters: Xylem Vulnerability Curves
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
library(fitplc)

# Load custom functions ----
set.seed(123)
n_cores = parallel::detectCores()
'%ni%' <- Negate('%in%')
source(file = "scripts/00_GammaModeFunction.R")

# Load Data ----
plc_curves <- read.csv("data_clean/CavitronDataClean.csv",row.names = "X")# load clean cavitron data
plc_curves <- plc_curves[,c("Sample_ref_1","Sample_ref_2","Raw_conductance_kg_Mpa_s","Conductivity_SI_corrT","Pressure_Mpa","kmax","Ksmax","PLCclean")]# select columns of interest

# Extract Tree.No and sample Replicate
plc_curves$Tree.No <- trunc(plc_curves$Sample_ref_1)# extract tree number
plc_curves$Rep <- plc_curves$Sample_ref_1# extract replicate number
plc_curves$Rep <- ifelse(test = (plc_curves$Rep - plc_curves$Tree.No)*10==0, yes = 1, no = (plc_curves$Rep - plc_curves$Tree.No)*10)# fix replicate number

# merge data with genotype information
sample_id <- read.csv("data_clean/dbg_cottonwood_plantID.csv")# load sample id information
plc_curves <- merge(x = plc_curves, y = sample_id, by = "Tree.No", all.x = TRUE)# merge sample id information with plc data
plc_curves$sample.id <- paste(plc_curves$Tree.No,plc_curves$Rep,sep = "-")# create sample id
head(plc_curves)

# Quick plot of the data
ggplot(plc_curves, aes(x = Pressure_Mpa, y = PLCclean, color = Population,group = sample.id)) + 
  facet_wrap(~Population) +
  geom_point(color="#000000") + 
  geom_line(stat = "smooth",alpha=0.5,lwd=1)+
  theme_bw()

# JAGS model ----
# P12 and P50 distributions ----------------------------------------------------
# Bayesian hierarchical model --------------------------------------------------
xylemvc.hmodel <- "model {
  # Define “expected” process according to the reparameterized - line 2
  # Weibull function: - line 3
  for(i in 1:N){
    # Likelihood for observed data: - line 5
    KSrel[i] ~ dnorm(muK[i],tau)
    # Model for replicated data: - line 7
    KSrep[i] ~ dnorm(muK[i],tau)
    
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
  
  # Compute Bayesian R2 value
  var.pred <- pow(sd(muK[]), 2)
  var.resid <- 1/tau
  R2 <- 1 - var.pred/(var.pred+var.resid)
  
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
  tau ~ dgamma(0.01,0.01)
  sig <- pow(tau,-0.5)
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
plc_curves$KSrel <- 1-plc_curves$PLCclean/100
ggplot(plc_curves, aes(x = Pressure_Mpa*-1, y = KSrel, color = Population,group = sample.id)) + 
  facet_wrap(~Population) +
  geom_point(color="#000000") + 
  geom_line(stat = "smooth",alpha=0.5,lwd=1)+
  theme_bw()
row.names(plc_curves) <- NULL
plc_curves <- droplevels(plc_curves)
head(plc_curves)

####### b. data list for JAGS --------------------------------------------------
# xylem vc data
# relative conductivity and water potentials
KSrel <- plc_curves$KSrel
plc_curves$WP.MPa <- plc_curves$Pressure_Mpa
MPa <- plc_curves$Pressure_Mpa*-1
N <- length(MPa)

# level 1 : sample-level id
plc_curves$ID <- factor(plc_curves$sample.id)
ID <- as.numeric(plc_curves$ID)
NTree <- length(unique(ID))

# level 2: assemble model matrix for elevation band effects
plc_curves$PopulationF <- factor(plc_curves$Population,levels = c("CCR-COL","JLA-JAK","NRV-NEW","TSZ-SAN"),
                                 labels = c("CCR-COL","JLA-JAK","NRV-NEW","TSZ-SAN"))
Xmat <- model.matrix(~PopulationF-1,plyr::ddply(plc_curves[,c("WP.MPa","ID","PopulationF")],~ID,plyr::catcolwise(unique)))
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
parameters <- c("P50.global","S50.global","P12.global","P50.TRT","S50.TRT","P12.TRT","KSrep")

# initial values
mean(sample(x = abs(populus.p50[!is.na(populus.p50$P50..MPa.),"P50..MPa."]),
            size = length(abs(populus.p50[!is.na(populus.p50$P50..MPa.),"P50..MPa."])),replace = T))
inits <- list()
inits[[1]] <- list(Px.0 = 1.56, Sx.0 = 25)  
inits[[2]] <- list(Px.0 = 1.65, Sx.0 = 28) 
inits[[3]] <- list(Px.0 = 1.49, Sx.0 = 31) 

# chains specification. Note: modify until ESS > 10000
thinSteps <- 10
numSavedSteps <- 7500
nChains <- 3
nIter <- ceiling( ( numSavedSteps * thinSteps ) / nChains ) # estimate desired number of iterations
burnInSteps <- floor(nIter/nChains)
adaptSteps <- round(burnInSteps*0.25)

# run model
output <- run.jags(model = xylemvc.hmodel,
                        monitor = parameters,
                        data = dataList,
                        n.chains = nChains,
                        burnin = burnInSteps,
                        sample = numSavedSteps,
                        adapt = adaptSteps,
                        thin = thinSteps,
                        method = "rjags",inits = inits)

# make sure rhat ~ 1.0001 and ESS > 10000
summary(output)
vc.mcmc <- output$mcmc
dim(vc.mcmc[[1]])
bayesplot::mcmc_trace(vc.mcmc,pars = c("P50.global","P50.TRT[1]","P50.TRT[2]","P50.TRT[3]","P50.TRT[4]"))+theme_bw()
bayesplot::mcmc_trace(vc.mcmc,pars = c("P12.global","P12.TRT[1]","P12.TRT[2]","P12.TRT[3]","P12.TRT[4]"))+theme_bw()
bayesplot::mcmc_acf(vc.mcmc,pars = c("P50.global","P50.TRT[1]","P50.TRT[2]","P50.TRT[3]","P50.TRT[4]"))+theme_bw()
bayesplot::mcmc_acf(vc.mcmc,pars =c("P12.global","P12.TRT[1]","P12.TRT[2]","P12.TRT[3]","P12.TRT[4]"))+theme_bw()


# Plot global xylem v-curves and predicted data --------------------------------
estimates <- data.frame(bayestestR::point_estimate(vc.mcmc,c("Median","Mean"),dispersion=T),bayestestR::ci(x=vc.mcmc,ci=0.90,method = "ETI"))
head(estimates)

mpo_pipo <- lm(plc_curves$KSrel~estimates[16:length(rownames(estimates)),"Median"])
summary(mpo_pipo)# 0.97

global_vcurve <- ggplot()+
  geom_point(data = plc_curves,mapping = aes(y = KSrel,x = WP.MPa*-1),alpha=0.45,size=2.5)+
  geom_line(mapping = aes(x = seq(from=0.5,to=5,length.out=1000),
                          y = fweibull(P = seq(from=0.5,to=5,length.out=1000),
                                       SX = estimates[estimates$Parameter == "S50.global","Median"],
                                       PX = estimates[estimates$Parameter == "P50.global","Median"],X = 50)),
            lwd=1.25,color="#3366FF")+
  geom_rect(mapping = aes(xmin=estimates[estimates$Parameter == "P50.global","Median"]-estimates[estimates$Parameter == "P50.global","MAD"],
                          xmax=estimates[estimates$Parameter == "P50.global","Median"]+estimates[estimates$Parameter == "P50.global","MAD"],
                          ymin=-Inf, ymax=Inf),fill="red2",alpha=0.15)+
  geom_vline(mapping = aes(xintercept = estimates[estimates$Parameter == "P50.global","Median"]),lwd=1.25,col="red2")+
  annotate(geom = "text",x = 7.5,y = 0.85,label = paste("P12 =",round(estimates[estimates$Parameter == "P12.global","Median"],digits = 2),sep=" "),size=5)+
  annotate(geom = "text",x = 7.5,y = 0.75,label = paste("P50 =",round(estimates[estimates$Parameter == "P50.global","Median"],digits = 2),sep=" "),size=5)+
  coord_cartesian(xlim = c(-0.05,10.75),ylim = c(-0.05,1.05),expand = F)+
  ylab(label = "Relative conductivity (0 - 1)")+
  xlab(label = "Water potential (-MPa)")+
  labs(title = "Populus fremontii")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size=12),
        plot.title = element_text(size = 16,hjust = 0.5,face = "bold"),
        panel.border = element_rect(linewidth = 1.2),
        panel.grid = element_blank())

obs_pred_pipo <- ggplot(mapping = aes(y=plc_curves$KSrel,x=estimates[16:length(rownames(estimates)),"Median"]))+
  geom_point(size=2.5,alpha=0.45)+
  geom_abline(slope = 1,linetype=2)+
  stat_poly_line(lwd=1.25)+
  stat_poly_eq()+
  ylab(label = expression(Obs.~K[rel]~("0 - 1")))+
  xlab(label = expression(Mod.~K[rel]~("0 - 1")))+
  ylim(c(-0.05,1.05))+
  xlim(c(-0.05,1.05))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size=12),
        plot.title = element_text(size = 16,hjust = 0.5,face = "bold"),
        panel.border = element_rect(linewidth = 1.2),
        panel.grid = element_blank())