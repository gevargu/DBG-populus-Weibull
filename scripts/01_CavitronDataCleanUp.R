# Load libraries
require(ggplot2)
require(tidyverse)
require(cowplot)
library(fitplc)
require(RColorBrewer)
library(readxl)

mypal <- brewer.pal(n=8, "Dark2")
palette(mypal)
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} 


######### Load raw Cavitron Data ####################
datadir <- "data_raw/RawCavitronData" # Directory where raw data is stored
data.files <- dir(path=datadir)# List of files in the directory
data.files <- data.files[grep("cavitron", data.files)]# Select only files with cavitron in the name
fulldata <- data.frame()# Initialize empty data frame to store all data

for(i in 1:length(data.files)) {# Loop through all files
  if(i==1){# If it is the first file, read it in and store it in fulldata
    fulldata0 <- read.csv(paste(datadir, data.files[i],sep="/"),skip = 1)# Read in the data
    finalpoint <- tail(fulldata0)[1:2,] # add in a -5.5 point for curve fitting
    finalpoint$Pressure_Mpa <- -5.5# Set the pressure to -5.5
    finalpoint$Raw_conductance_kg_Mpa_s <- 0# Set the conductance to 0
    finalpoint$Conductivity_SI_corrT <- 0# Set the conductivity to 0
    fulldata0 <- rbind(fulldata0, finalpoint)# Add the final point to the data
    fulldata <- fulldata0[-which(is.na(fulldata0$Pressure_Mpa) | fulldata0$Note=="bad" ),c(1:43)]# Remove any rows with NA or bad data
  }else{# If it is not the first file, read it in and append it to fulldata
    dataz0 <- read.csv(paste(datadir, data.files[i],sep="/"),skip = 1)# Read in the data
    finalpoint <- tail(dataz0)[1:2,] # add in a -5.5 point for curve fitting
    finalpoint$Pressure_Mpa <- -5.5# Set the pressure to -5.5
    finalpoint$Raw_conductance_kg_Mpa_s <- 0# Set the conductance to 0
    finalpoint$Conductivity_SI_corrT <- 0# Set the conductivity to 0
    dataz0 <- rbind(dataz0, finalpoint)# Add the final point to the data
    dataz <- dataz0[-which(is.na(dataz0$Pressure_Mpa) | dataz0$Note=="bad" ),c(1:43)]# Remove any rows with NA or bad data
    fulldata <- rbind(fulldata, dataz)# Append the data to fulldata
  }# End if
}# End for

# remove duplicates
fulldata <- fulldata[!duplicated(fulldata$Conductivity_SI_corrT),]
# somehow there were a bunch of dups in there.


## first let's clean Kmax values
# get rid of 0 values at low tensions
fulldata$Conductivity_SI_corrT[which(fulldata$Raw_conductance_kg_Mpa_s<= 0 & fulldata$Pressure_Mpa> -2)] <- NA
fulldata$Raw_conductance_kg_Mpa_s[which(fulldata$Raw_conductance_kg_Mpa_s<= 0 & fulldata$Pressure_Mpa> -2)] <- NA


# currently let's just turn high tension negative ks into 0
fulldata$Conductivity_SI_corrT[which(fulldata$Conductivity_SI_corrT<= 0 & fulldata$Pressure_Mpa < -2)] <- 0
fulldata$Raw_conductance_kg_Mpa_s[which(fulldata$Raw_conductance_kg_Mpa_s<= 0 & fulldata$Pressure_Mpa < -2)] <- 0

#### Code for pulling out specific curves that you're interested in
# just curves with low conductance values that otherwise don't visualize well
lowvals <- c(2,3,12,21,30,32,44,45,59,60,62,68,72,83) # [which(fulldata$Sample_ref_1 %in% lowvals),]
# just curves processed in final batch (3rd batch)
newvals <- c(7.2, 7.3,
             7.4,
             8.2,
             44.2,
             45.3,
             49.2,
             49.3,
             61.2,
             63.2,
             63.3)

# some visualizations of things I was worried about
ggplot(fulldata[which(fulldata$Sample_ref_1 %in% newvals),], aes(x=Pressure_Mpa, y=Raw_conductance_kg_Mpa_s, col=factor(Sample_ref_1))) + geom_point() + geom_smooth(se=F) + 
  facet_wrap(~Sample_ref_1) + theme(legend.position = "none")
# actually, all the new curves look pretty damn good. except maybe 63.3

# compare conductance and conductivity
p1 <- ggplot(fulldata, aes(x=Pressure_Mpa, y=Raw_conductance_kg_Mpa_s, col=factor(Sample_ref_1))) + geom_point() + geom_smooth(se=F) + theme(legend.position = "none")

p2 <- ggplot(fulldata, aes(x=Pressure_Mpa, y=Conductivity_SI_corrT, col=factor(Sample_ref_1))) + geom_point() + geom_smooth(se=F)  + theme(legend.position = "none")


multiplot(p1, p2)
# weirdly more normally distributed variation in raw conductance than conductivity. But conductivity is probably driven by 7 high outliers


# plot all curves. 
# in one plotting device
ggplot(fulldata, aes(x=Pressure_Mpa, y=Raw_conductance_kg_Mpa_s, col=factor(Sample_ref_1)))+ geom_point() + 
  geom_smooth(se=F) + xlim(-2.5,0) + facet_wrap(~Sample_ref_1) + theme(legend.position = "none")

for(j in unique(fulldata$Sample_ref_1)){
  print(ggplot(fulldata[which(fulldata$Sample_ref_1==j),], aes(x=Pressure_Mpa, y=Raw_conductance_kg_Mpa_s, col=factor(Sample_ref_1)))+ geom_point()+theme(legend.position = "none")+ggtitle(j))
}

############# Curve Cleaning #################
####### Flag bad curves and bad points ########
fulldata$flag <- 0

# bad (never moved water)
# 9, 50, 7.2 sample_ref_2 == 1 (removed data file)
fulldata$flag[which(fulldata$Sample_ref_1=="9")] <- 1
fulldata$flag[which(fulldata$Sample_ref_1=="50")] <- 1

# 85: increased until 2
fulldata$flag[which(fulldata$Sample_ref_1=="85" & fulldata$Pressure_Mpa > -1.75 )] <- 1
# 76: increased until 2.25, removed above -1.5
fulldata$flag[which(fulldata$Sample_ref_1=="76" & fulldata$Pressure_Mpa > -1.55 )] <- 1
# 75.2 decrased smoothly until 1.5, platueas then decreases normally after 2.25. Removed -1 and -1.25
fulldata$flag[which(fulldata$Sample_ref_1=="75.2" & fulldata$Pressure_Mpa > -1.35 )] <- 1
# 75 bit funky above 2, leaving as is
# 70 decrased smoothly until 2.25 then decreases normally after 2.25. left as is
# 72: points above -1.5 bad
fulldata$flag[which(fulldata$Sample_ref_1=="72"& fulldata$Pressure_Mpa> -1.5 )] <- 1
# 7: -1 point looks a bit high
fulldata$flag[which(fulldata$Sample_ref_1=="7"& fulldata$Pressure_Mpa> -1.1 )] <- 1
# 69: above -1.5 high
fulldata$flag[which(fulldata$Sample_ref_1=="69"& fulldata$Pressure_Mpa> -1.5 )] <- 1
# 68: -1 points too high
fulldata$flag[which(fulldata$Sample_ref_1=="68"& fulldata$Pressure_Mpa> -1.1 )] <- 1
#67 -1.5 is bad
fulldata$flag[which(fulldata$Sample_ref_1=="67"& fulldata$Pressure_Mpa == -1.505 )] <- 1
# 66: above -1.5 too high
fulldata$flag[which(fulldata$Sample_ref_1=="66"& fulldata$Pressure_Mpa> -1.5 )] <- 1
# 63: slow decline from 1?
# ***60: slow incline to 2.25? Maybe remove >1.5?
fulldata$flag[which(fulldata$Sample_ref_1=="60"& fulldata$Pressure_Mpa> -1.5 )] <- 1
# 58: -1.25 bad
fulldata$flag[which(fulldata$Sample_ref_1=="58"& fulldata$Pressure_Mpa< -1.1 & fulldata$Pressure_Mpa>-1.4)] <- 1
# 51: 1.25 and 1 noisy
fulldata$flag[which(fulldata$Sample_ref_1=="51"& fulldata$Pressure_Mpa> -1.5 )] <- 1
# 47: -1 bad
fulldata$flag[which(fulldata$Sample_ref_1=="47"& fulldata$Pressure_Mpa> -1.1 )] <- 1
# 46: -1 bad
fulldata$flag[which(fulldata$Sample_ref_1=="46"& fulldata$Pressure_Mpa> -1.1 )] <- 1
# ***45: a total cluster. 3.25 has a bad one, and everything below 2.5 looks weird # left in
#fulldata$Raw_conductance_kg_Mpa_s[which(fulldata$Sample_ref_1=="46"& fulldata$Pressure_Mpa< -3 )] #<- 1
# ***44: total cluster, increases sharply until -3 then dies # left in
# 41.2: -1.5 and above totally different from -1.76 up
fulldata$flag[which(fulldata$Sample_ref_1=="41.2"& fulldata$Pressure_Mpa> -1.6 )] <- 1
# 41: similar to 41.2
fulldata$flag[which(fulldata$Sample_ref_1=="41"& fulldata$Pressure_Mpa> -1.6 )] <- 1
# 38: -1 bad
fulldata$flag[which(fulldata$Sample_ref_1=="38"& fulldata$Pressure_Mpa> -1.1 )] <- 1
# 36: -1 and one high -1.25 bad
fulldata$flag[which(fulldata$Sample_ref_1=="36"& fulldata$Pressure_Mpa> -1.3 & fulldata$Raw_conductance_kg_Mpa_s>1.64e-05)] <- 1
# 37: -1 bad
fulldata$flag[which(fulldata$Sample_ref_1=="41.2"& fulldata$Pressure_Mpa> -1 )] <- 1
# 34: -1 and -1.25 real high
fulldata$flag[which(fulldata$Sample_ref_1=="34"& fulldata$Pressure_Mpa> -1.3 )] <- 1
#32.2: -1 is bad
fulldata$flag[which(fulldata$Sample_ref_1=="32.2"& fulldata$Pressure_Mpa> -1 )] <- 1
# ***33: BATSHIT. increases until 2.75 then crashes. # left in
# ***32: BATSHIT. increases until 2.75 then crashes. # left in, but gets pulled out after PLC
# ***31.2. REAL bad, increases until 2.5, and all over the place # currently left in
#fulldata$flag[which(fulldata$Sample_ref_1=="31.2")] <- 1
# 31.1 - < -1.5 prob bad (could be in 'constant decline cat)
fulldata$flag[which(fulldata$Sample_ref_1=="31.1"& fulldata$Pressure_Mpa> -1.3 )] <- 1
# ***3: BATSHIT. Increases until 2.6 and then crashes-ish # left in currently, but gets pulled after PLC
#fulldata$flag[which(fulldata$Sample_ref_1=="3")] <- 1
# 22: 1, maybe 1.25 bad
fulldata$flag[which(fulldata$Sample_ref_1=="22"& fulldata$Pressure_Mpa> -1.1 )] <- 1
# 21: 1 bad
fulldata$flag[which(fulldata$Sample_ref_1=="21"& fulldata$Pressure_Mpa> -1.1 )] <- 1
# 71: 1 and 1.25 bad
fulldata$flag[which(fulldata$Sample_ref_1=="71"& fulldata$Pressure_Mpa> -1.3 )] <- 1
# 62: 1 and 1.25 bad
fulldata$flag[which(fulldata$Sample_ref_1=="62"& fulldata$Pressure_Mpa> -1.3 )] <- 1
# 83: 1 and 1.25 bad
fulldata$flag[which(fulldata$Sample_ref_1=="83"& fulldata$Pressure_Mpa> -1.3 )] <- 1
# 2: wow, 1.25,1.5, prob 1.75 bad
fulldata$flag[which(fulldata$Sample_ref_1=="2"& fulldata$Pressure_Mpa< -1.1 & fulldata$Pressure_Mpa>-1.7)] <- 1
# 19: 1 bad
fulldata$flag[which(fulldata$Sample_ref_1=="19"& fulldata$Pressure_Mpa> -1.1 )] <- 1
# 17: -1.25 and up prob bad
fulldata$flag[which(fulldata$Sample_ref_1=="17"& fulldata$Pressure_Mpa> -1.3 )] <- 1
# 12: 1,1.25 and most of 1.5 bad, plus a -3 with way too high k
fulldata$flag[which(fulldata$Sample_ref_1=="12"& fulldata$Pressure_Mpa> -1.6 & fulldata$Raw_conductance_kg_Mpa_s<4e-6 )] <- 1
fulldata$flag[which(fulldata$Sample_ref_1=="12"& fulldata$Pressure_Mpa< -2.9 & fulldata$Raw_conductance_kg_Mpa_s>0 )] <- 1
# 10: -1.25 probs bad
fulldata$flag[which(fulldata$Sample_ref_1=="10"& fulldata$Pressure_Mpa == -1.221 )] <- 1
# 1.2: -1.5 and up clearly bad
fulldata$flag[which(fulldata$Sample_ref_1=="1.2"& fulldata$Pressure_Mpa> -1.6 )] <- 1
# 9.2: -1 bad
fulldata$flag[which(fulldata$Sample_ref_1=="9.2"& fulldata$Pressure_Mpa> -1 )] <- 1
# 7.2: original didn't move water
fulldata$flag[which(fulldata$Sample_ref_1=="7.2"& fulldata$Sample_ref_2==1)] <- 1
# 63.3: Rises until -1.75, then -1.75 and -2 and -2.25 look similar
fulldata$flag[which(fulldata$Sample_ref_1=="63.3"& fulldata$Pressure_Mpa> -1.7 )] <- 1


# new samples from last batch:
# 7.2again
# 7.3 - Knat
# 7.4 - Knat
# 8.2
# 44.2
# 45.3
# 49.2
# 49.3 - Knat
# 61.2
# 63.2
# 63.3 - Knat


# remove data that was flagged
dataclean <- fulldata[which(fulldata$flag==0 & !is.na(fulldata$Raw_conductance_kg_Mpa_s)),]

ggplot(dataclean, aes(x=Pressure_Mpa, y=Raw_conductance_kg_Mpa_s, col=factor(Sample_ref_1)))+ geom_point() + 
  geom_smooth(se=F) + xlim(-2,0) + facet_wrap(~Sample_ref_1) + theme(legend.position = "none")

ggplot(dataclean, aes(x=Pressure_Mpa, y=Raw_conductance_kg_Mpa_s, col=factor(Sample_ref_1)))+ geom_point() + 
  geom_smooth(se=F)  + facet_wrap(~Sample_ref_1) + theme(legend.position = "none")

for(j in unique(dataclean$Sample_ref_1)){
  print(ggplot(dataclean[which(dataclean$Sample_ref_1==j),], aes(x=Pressure_Mpa, y=Raw_conductance_kg_Mpa_s, col=factor(Sample_ref_1)))+ geom_point()+theme(legend.position = "none")+ggtitle(j))
}

###### calculate Kmax for each branch #############
# then recalculate PLC by hand (don't trust Cavisoft version)
# we'll use everything above -2, and try upper quartile to get Kmax. 
# Note, for some branches Kmax actually seemed to be reached at 2 or 2.5, but some branches were losing by 2
# 7,64,38,49, 63.3 - 1.75 for kmax,
# 63, 61.2 slow decline from 1
# 60,71 - increase to 2.25,2
# 32,3,2 - increase to 2.75


ksmax.prob <- 0.75 # taking the 75th percentile
kmax1 <- dataclean[which(dataclean$Sample_ref_1 %in% c(7,64,38,49,63)),] %>% filter(Pressure_Mpa>-1.8) %>% group_by(Sample_ref_1) %>% summarise(kmax = quantile(Raw_conductance_kg_Mpa_s, probs = ksmax.prob, na.rm=T), Ksmax = quantile(Conductivity_SI_corrT, probs = ksmax.prob, na.rm=T))
kmax2 <- dataclean[which(dataclean$Sample_ref_1 %in% c(60,71)),] %>% filter(Pressure_Mpa< -1.9 & Pressure_Mpa> -2.3) %>% group_by(Sample_ref_1) %>% summarise(kmax = quantile(Raw_conductance_kg_Mpa_s, probs = ksmax.prob, na.rm=T), Ksmax = quantile(Conductivity_SI_corrT, probs = ksmax.prob, na.rm=T))
kmax3 <- dataclean[which(dataclean$Sample_ref_1 %in% c(32,3,2)),] %>% filter(Pressure_Mpa< -2.2 & Pressure_Mpa> -2.8) %>% group_by(Sample_ref_1) %>% summarise(kmax = quantile(Raw_conductance_kg_Mpa_s, probs = ksmax.prob, na.rm=T), Ksmax = quantile(Conductivity_SI_corrT, probs = ksmax.prob, na.rm=T))
kmax4 <- dataclean[which(!dataclean$Sample_ref_1 %in% c(7,64,38,49,63,60,71,32,3,2)),] %>% filter(Pressure_Mpa>-1.8) %>% group_by(Sample_ref_1) %>% summarise(kmax = quantile(Raw_conductance_kg_Mpa_s, probs = ksmax.prob, na.rm=T), Ksmax = quantile(Conductivity_SI_corrT, probs = ksmax.prob, na.rm=T))
kmax <- rbind(kmax1, kmax2, kmax3, kmax4)


datac <- left_join(dataclean, kmax)
datac$PLCclean <- (1 - (datac$Raw_conductance_kg_Mpa_s /datac$kmax))*100
datac$PLCclean <- (1 - (datac$Conductivity_SI_corrT/datac$Ksmax))*100
# let's force things to be between 0 and 100
datac$PLCclean[which(datac$PLCclean<0)] <- 0
# and get rid of the bad early values for many plants
datac$PLCclean[which(datac$PLCclean > 20 & datac$Pressure_Mpa> -1.75)] <- 0

ggplot(datac, aes(x=Pressure_Mpa, y=PLCclean, col=factor(Sample_ref_1))) + geom_point() + geom_smooth(se=F) + theme(legend.position = "none")

ggplot(datac, aes(x=Pressure_Mpa, y=PLCclean, col=factor(Sample_ref_1)))+ geom_point() + 
  geom_smooth(se=F)  + facet_wrap(~Sample_ref_1) + theme(legend.position = "none")

#### pulling out bad curves based on PLC
badbranches <- c("3","31.2","32","33","85") # 44 is the only questionable one. But doesn't fit with fitplc
bad <- datac %>% filter(Sample_ref_1 %in% badbranches)
good <- datac %>% filter(!Sample_ref_1 %in% badbranches)

ggplot(good, aes(x=Pressure_Mpa, y=PLCclean, col=factor(Sample_ref_1)))+ geom_point() + 
  geom_smooth(se=F)  + facet_wrap(~Sample_ref_1) + theme(legend.position = "none")

# save the data--------------------------------
head(good)
write.csv(good, "data_clean/CavitronDataClean.csv")
