#
# Plot model fit and parameters

#' # Model parameters: Xylem Vulnerability Curves
#' 
#' ### Author: German Vargas G.
#' 
# Load Libraries ----
library(coda)
library(broom.mixed)
library(udunits2)
library(tidyverse)

# Custom functions ----
calculate_P12 <- function(px, sx, X=50) {
  # Define V
  V <- (X - 100) * log((1 - X / 100))
  px * (log(0.88) / log(0.5))^(V / (px * sx))
}

weibull.params <- function(x1,x2,px1,px2){
  x1 = x1/100# proportion plc
  x2 = x2/100# proportion plc
  px1 = px1*-1
  px2 = px2*-1
  
  # calculate weibull parameter c
  c = (log(x = log(x = 1-x1)/log(x = 1-x2)))/(log(x = px1)-log(x = px2))
  # calculate weibull parameter b
  b = px1/(-log(x = 1-x1))^(1/c)
  return(data.frame(c=c,b=b))
}

Weibull.PLC <- function(WP,b,c){
  PLC <- (100*(1-exp(-((-WP/b)^c))))
  return(PLC)
}

# Load data ----
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
plc_curves2$KSrel <- 1-plc_curves2$PLCclean/100
row.names(plc_curves2) <- NULL
plc_curves2 <- droplevels(plc_curves2)
head(plc_curves2)
# level 1 : sample-level id
plc_curves2$ID <- factor(plc_curves2$sample.id,levels = unique(plc_curves2$sample.id))
# level 2: assemble model matrix for elevation band effects
plc_curves2$PopulationF <- factor(plc_curves2$Population,levels = unique(plc_curves2$Population),
                                  labels = unique(plc_curves2$Population))
Xmat <- model.matrix(~PopulationF-1,plyr::ddply(plc_curves2[,c("Pressure_Mpa","ID","PopulationF")],~ID,plyr::catcolwise(unique)))

# Quick plot of the data
ggplot(plc_curves, aes(x=Pressure_Mpa, y=PLCclean, col=factor(sample.id)))+ geom_point() + 
  geom_smooth(se=F)  + facet_wrap(~sample.id) +theme_bw() + theme(legend.position = "none") 
ggplot(plc_curves2, aes(x=Pressure_Mpa, y=PLCclean, col=factor(sample.id)))+ geom_point() + 
  geom_smooth(se=F)  + facet_wrap(~sample.id) +theme_bw() + theme(legend.position = "none") 

# Load codas ----
load(file = "output/coda/jm_coda.Rdata")

# Tidy parameters
sum_param <- tidyMCMC(jm_coda, 
                      conf.int =  TRUE, 
                      conf.method = "HPDinterval",
                      conf.level = 0.95)

# only predictions
kpred <- sum_param[grep("KSrep",sum_param$term),]
kpred <- cbind(plc_curves2, kpred)
head(kpred)

plot(kpred$Pressure_Mpa, kpred$estimate)

# Calculate model fit ----
##### All data -----
m1 <- lm(KSrel ~ estimate, data = kpred)
summary(m1) # R2 = 0.9778!!!!! great fit!

#### Calculate R2, coverage, and bias for each Population ----
pred2 <- kpred %>%
  mutate(cover = if_else(KSrel >= conf.low & KSrel <= conf.high, 1, 0))

# total data frame
tot.df <- data.frame(Population = "all",
                     R2 = summary(m1)$adj.r.squared,
                     coverage = mean(pred2$cover),
                     bias = summary(m1)$coef[2,1])

# for each TRT
trt.df <- data.frame(Population = levels(as.factor(pred2$Population)),
                     R2 = pred2 %>%
                       group_by(Population) %>%
                       group_map(~ summary(lm(estimate ~ KSrel, data = .x))$adj.r.squared) %>%
                       purrr:::simplify(),
                     coverage = pred2 %>%
                       group_by(Population) %>%
                       summarize(coverage = mean(cover)) %>%
                       pull(coverage),
                     bias = pred2 %>%
                       group_by(Population) %>%
                       group_map(~ summary(lm(estimate ~ KSrel, data = .x))$coef[2,1]) %>%
                       purrr:::simplify())

all.df <- bind_rows(trt.df, tot.df)

# write out model fits
write_csv(all.df, "output/model_fit.csv")

#### Plot model fit -----

# PLC curve: observed and fitted
#png(filename = "output/figures/plc_obs_pred.png", width = 8, height = 7, units = "in", res = 600)
kpred %>% 
  ggplot(aes(x = Pressure_Mpa)) +
  geom_pointrange(aes(y = estimate, 
                      ymin = conf.low,
                      ymax = conf.high,
                      color = "Predicted"),alpha=0.5) +
  geom_point(aes(y = KSrel, 
                 color = "Observed"),alpha=0.65) +
  facet_wrap(~Population) +
  ylab(label="Relative Conductivity") +
  xlab(label="Xylem Pressure (MPa)") +
  scale_color_manual("",values = c("#ED4D40","#1E3550"))+
  ggpubr::theme_pubr()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11,face="bold"))
#dev.off()

# Plot Observed vs. fitted
trt.df.labels <- trt.df |> 
  mutate(lab = paste0("R^2==", round(R2, 3)))

#png(filename = "output/figures/obs_pred.png", width = 8, height = 7, units = "in", res = 600)
kpred %>%
  ggplot(aes(x = KSrel, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low, 
                      ymax = conf.high,
                      color = ""),alpha=0.45) +
  geom_abline(slope = 1, intercept = 0, lty = 1,
              color = "#ED4D40",lwd=0.9) +
  geom_text(data = trt.df.labels, 
            aes(x = 0.22, y = 1,
                label = lab),
            parse = TRUE,
            hjust = 0,
            vjust = 0) +
  scale_color_manual("",values = c("#1E3550"))+
  ylab(label="Estimated Relative Conductivity") +
  xlab(label="Observed Relative Conductivity") +
  facet_wrap(~Population) +
  ggpubr::theme_pubr()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12,face="bold"),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
#dev.off()


# Weibull parameters ------
species <- sum_param %>%
  filter(grepl("P12\\.g", term) |
           grepl("P50\\.g", term)) %>%
  mutate(parameter = sub("\\.global", "", term),
         Population = "all") %>%
  pivot_wider(id_cols = c(Population), # Columns to keep as identifiers
              names_from = parameter,
              values_from = c(estimate, conf.low, conf.high)) %>%
  mutate(weibull_params = purrr::pmap(list(x1 = 12,
                                           x2 = 50,
                                           px1 = (estimate_P12) * -1,
                                           px2 = (estimate_P50) * -1),
                                      ~ weibull.params(x1 = .x, x2 = .y, px1 = ..3, px2 = ..4)),
         weibull_params.low = purrr::pmap(list(x1 = 12,
                                               x2 = 50,
                                               px1 = (conf.low_P12) * -1,
                                               px2 = (conf.low_P50) * -1),
                                          ~ weibull.params(x1 = .x, x2 = .y, px1 = ..3, px2 = ..4)),
         weibull_params.high = purrr::pmap(list(x1 = 12,
                                                x2 = 50,
                                                px1 = (conf.high_P12) * -1,
                                                px2 = (conf.high_P50) * -1),
                                           ~ weibull.params(x1 = .x, x2 = .y, px1 = ..3, px2 = ..4)),
         c = purrr::map_dbl(weibull_params, "c"),
         b = purrr::map_dbl(weibull_params, "b"),
         conf.low_c = purrr::map_dbl(weibull_params.low, "c"),
         conf.low_b = purrr::map_dbl(weibull_params.low, "b"),
         conf.high_c = purrr::map_dbl(weibull_params.high, "c"),
         conf.high_b = purrr::map_dbl(weibull_params.high, "b")) %>% 
  select(-weibull_params,-weibull_params.low,-weibull_params.high)


population <- sum_param %>%
  filter(grepl("P12\\.T", term) | grepl("P50\\.T", term)) %>%
  mutate(
    parameter = sub("\\[[0-9]\\]", "", term),
    parameter = sub("\\.TRT", "", parameter),
    Population = case_when(
      grepl("P12\\.T[^0-9]*1($|[^0-9])", term) ~ "CCR-COL",
      grepl("P12\\.T[^0-9]*2($|[^0-9])", term) ~ "JLA-JAK",
      grepl("P12\\.T[^0-9]*3($|[^0-9])", term) ~ "NRV-NEW",
      grepl("P12\\.T[^0-9]*4($|[^0-9])", term) ~ "TSZ-SAN",
      grepl("P50\\.T[^0-9]*1($|[^0-9])", term) ~ "CCR-COL",
      grepl("P50\\.T[^0-9]*2($|[^0-9])", term) ~ "JLA-JAK",
      grepl("P50\\.T[^0-9]*3($|[^0-9])", term) ~ "NRV-NEW",
      grepl("P50\\.T[^0-9]*4($|[^0-9])", term) ~ "TSZ-SAN",
      TRUE ~ NA_character_),
    parameter = factor(parameter, levels = c("P12", "P50"))) %>%
  pivot_wider(id_cols = c(Population), # Columns to keep as identifiers
              names_from = parameter,
              values_from = c(estimate, conf.low, conf.high)) %>%
  mutate(weibull_params = purrr::pmap(list(x1 = 12,
                                           x2 = 50,
                                           px1 = (estimate_P12) * -1,
                                           px2 = (estimate_P50) * -1),
                                      ~ weibull.params(x1 = .x, x2 = .y, px1 = ..3, px2 = ..4)),
         weibull_params.low = purrr::pmap(list(x1 = 12,
                                               x2 = 50,
                                               px1 = (conf.low_P12) * -1,
                                               px2 = (conf.low_P50) * -1),
                                          ~ weibull.params(x1 = .x, x2 = .y, px1 = ..3, px2 = ..4)),
         weibull_params.high = purrr::pmap(list(x1 = 12,
                                                x2 = 50,
                                                px1 = (conf.high_P12) * -1,
                                                px2 = (conf.high_P50) * -1),
                                           ~ weibull.params(x1 = .x, x2 = .y, px1 = ..3, px2 = ..4)),
         c = purrr::map_dbl(weibull_params, "c"),
         b = purrr::map_dbl(weibull_params, "b"),
         conf.low_c = purrr::map_dbl(weibull_params.low, "c"),
         conf.low_b = purrr::map_dbl(weibull_params.low, "b"),
         conf.high_c = purrr::map_dbl(weibull_params.high, "c"),
         conf.high_b = purrr::map_dbl(weibull_params.high, "b")) %>% 
  select(-weibull_params,-weibull_params.low,-weibull_params.high) 

tree <- sum_param %>%
  filter(grepl("Tree\\.Sx", term) |
           grepl("Tree\\.Px", term)) %>%
  mutate(tree.no = as.numeric(gsub("Tree\\.(Px|Sx)\\[(\\d+)\\]", "\\2", term)),
         type= gsub("Tree\\.(Px|Sx)\\[(\\d+)\\]", "\\1", term),
         type= gsub(pattern = "Px",replacement = "P50",x = type)) %>%
  pivot_wider(id_cols = c(tree.no), # Columns to keep as identifiers
              names_from = type,
              values_from = c(estimate, std.error, conf.low, conf.high)) %>%
  mutate(estimate_P12= calculate_P12(px=estimate_P50,sx = estimate_Sx,X = 50),
         estimate_P12.max = calculate_P12(px=estimate_P50+std.error_P50,sx = estimate_Sx+std.error_Sx,X = 50),
         std.error_P12 = estimate_P12.max-estimate_P12,
         conf.low_P12 = calculate_P12(px=conf.low_P50,sx = conf.low_Sx,X = 50),
         conf.high_P12 = calculate_P12(px=conf.high_P50,sx = conf.high_Sx,X = 50)) %>%
  select(tree.no,estimate_P12,std.error_P12,conf.low_P12,conf.high_P12,estimate_P50,std.error_P50,conf.low_P50,conf.high_P50) %>%
  mutate(weibull_params = purrr::pmap(list(x1 = 12,
                                           x2 = 50,
                                           px1 = (estimate_P12) * -1,
                                           px2 = (estimate_P50) * -1),
                                      ~ weibull.params(x1 = .x, x2 = .y, px1 = ..3, px2 = ..4)),
         weibull_params.low = purrr::pmap(list(x1 = 12,
                                               x2 = 50,
                                               px1 = (conf.low_P12) * -1,
                                               px2 = (conf.low_P50) * -1),
                                          ~ weibull.params(x1 = .x, x2 = .y, px1 = ..3, px2 = ..4)),
         weibull_params.high = purrr::pmap(list(x1 = 12,
                                               x2 = 50,
                                               px1 = (conf.high_P12) * -1,
                                               px2 = (conf.high_P50) * -1),
                                          ~ weibull.params(x1 = .x, x2 = .y, px1 = ..3, px2 = ..4)),
    c = purrr::map_dbl(weibull_params, "c"),
    b = purrr::map_dbl(weibull_params, "b"),
    conf.low_c = purrr::map_dbl(weibull_params.low, "c"),
    conf.low_b = purrr::map_dbl(weibull_params.low, "b"),
    conf.high_c = purrr::map_dbl(weibull_params.high, "c"),
    conf.high_b = purrr::map_dbl(weibull_params.high, "b")) %>% 
  select(-weibull_params,-weibull_params.low,-weibull_params.high) %>%
  mutate(ID = levels(plc_curves2$ID)) 

# Determine the number of columns based on your populations
num_cols <- 4

# Reshape Xmat into a matrix with the correct number of columns
Xmat_matrix <- matrix(Xmat, ncol = num_cols, byrow = FALSE) # Adjust byrow if needed

# Convert the matrix to a tibble with appropriate column names
Xmat_df <- as_tibble(Xmat, .name_repair = "minimal")

tree <- tree %>%
  mutate(
    Population = case_when(
      Xmat_df$`PopulationFCCR-COL` == 1 ~ "CCR-COL",
      Xmat_df$`PopulationFJLA-JAK` == 1 ~ "JLA-JAK",
      Xmat_df$`PopulationFNRV-NEW` == 1 ~ "NRV-NEW",
      Xmat_df$`PopulationFTSZ-SAN` == 1 ~ "TSZ-SAN",
      TRUE ~ NA_character_ # Default if no match (shouldn't happen if all rows are covered)
    ),
    Popn.Geno = strsplit(ID, "-") %>% sapply("[", 1),
    Tree.No = strsplit(ID, "-") %>% sapply("[", 2),
    Rep.No = strsplit(ID, "-") %>% sapply("[", 3)) %>%
  select(Population,Popn.Geno,Tree.No,Rep.No,tree.no,ID,estimate_P12,estimate_P50,c,b,
         conf.low_P12,conf.low_P50,conf.low_c,conf.low_b,
         conf.high_P12,conf.high_P50,conf.high_c,conf.high_b) 

species
population
colnames(tree)[5] <- "bayes.tree.no"
write.csv(tree, "output/weibull_params_tree_level.csv", row.names = F)
write.csv(population, "output/weibull_params_population_level.csv", row.names = F)
write.csv(species, "output/weibull_params_species_level.csv", row.names = F)

# Plot curves global and per population -----
sum_param %>%
  filter(grepl("Sx\\.", term) |
           grepl("Px\\.", term)) %>%
  mutate(parameter = sub("\\.global", "", term),
         Population = "all")

ggplot()+
  geom_point(data = plc_curves2,
             aes(x = Pressure_Mpa, y = PLCclean),alpha=0.35,size=2.5,color="#1E3550")+
  geom_line(mapping = aes(x = seq(from=max(plc_curves2$Pressure_Mpa),to=min(plc_curves2$Pressure_Mpa),length.out=1000),
                          y = Weibull.PLC(WP = seq(from=max(plc_curves2$Pressure_Mpa),to=min(plc_curves2$Pressure_Mpa),length.out=1000),
                                          b = pop.weibull$b,
                                          c = pop.weibull$c)),lwd=1.25,color="#ED4D40")+
  geom_ribbon(mapping = aes(x = seq(from=max(plc_curves2$Pressure_Mpa),to=min(plc_curves2$Pressure_Mpa),length.out=1000),
                            ymin = Weibull.PLC(WP = seq(from=max(plc_curves2$Pressure_Mpa),to=min(plc_curves2$Pressure_Mpa),length.out=1000),
                                            b = 2.357248,
                                            c = 10.1471),
                            ymax = Weibull.PLC(WP = seq(from=max(plc_curves2$Pressure_Mpa),to=min(plc_curves2$Pressure_Mpa),length.out=1000),
                                               b = 3.104619,
                                               c = 13.46862)),lwd=1.25,fill="#ED4D40",alpha=0.2)+
  ylab(label = "PLC (%)")+
  xlab(label = "Xylem pressure (MPa)")+
  ggpubr::theme_pubr()+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size=12),
        plot.title = element_text(size = 16,hjust = 0.5,face = "bold"))


#### Plot model parameters by population ----
png(filename = "output/figures/Px_pars.png", width = 8, height = 7, units = "in", res = 600)
ggplot() +
  geom_pointrange(data = population,
                  aes(x = parameter,
                      y = estimate,
                      ymin = conf.low, 
                      ymax = conf.high,
                      color = Population),
                  position = position_dodge(width = 1)) +
  geom_pointrange(data = pop,
                  aes(x = parameter,
                      y = estimate,
                      ymin = conf.low, 
                      ymax = conf.high),color="#1E3550") +
  facet_wrap(~parameter, scales = "free") +
  scale_color_manual(values = c("#27355E","#3A8B91","#FD7C42","#FA4D31")) +
  ylab(label="Xylem Pressure (MPa)") +
  ggpubr::theme_pubr()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12,face="bold"),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()


ggplot() +
  geom_pointrange(data = tree[,],
                  aes(x = ID,
                      y = c,
                      ymin = conf.low_c, 
                      ymax = conf.high_c,
                      color = Population),
                  position = position_dodge(width = 1)) +
  facet_wrap(~Population,scales="free_x") +
  ggpubr::theme_pubr()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12,face="bold"),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

  geom_pointrange(data = pop,
                  aes(x = parameter,
                      y = estimate,
                      ymin = conf.low, 
                      ymax = conf.high),color="#1E3550") +
  facet_wrap(~parameter, scales = "free") +
  scale_color_manual(values = c("#27355E","#3A8B91","#FD7C42","#FA4D31")) +
  ylab(label="Xylem Pressure (MPa)") +
  ggpubr::theme_pubr()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12,face="bold"),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

