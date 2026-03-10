## BLL Model Selection Process ##

setwd("~/Research Materials/DepredationPotential_Overlay/DepredationCode")


## load all necessary libraries 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sdmTMB)
library(spdep)
library(pROC)


## read in BLLCoord data to use for model selection
BLLCoordSelect <- read.csv("data/processed/BLL_model_selection.csv")
BLLCoordSelect$X <- NULL

## Add UTM points
BLLCoordSelect <-sdmTMB::add_utm_columns(BLLCoordSelect,ll_names=c("Longitude","Latitude"),units="km", 
                                   utm_crs = 32617)
  # UTM Zone 17N

## create a mesh object 
mesh <- make_mesh(BLLCoordSelect, xy_cols = c("X", "Y"), cutoff = 20)
mesh$mesh$n # 263
plot(mesh)

### Shark Abundance Model selection 
## Comparing Intercept (Null), Depth linear, and depth smooth
SharkDensityNull <- sdmTMB(
  formula = TotalOver150 ~ 1, 
  data = BLLCoordSelect,
  mesh = mesh, 
  family = nbinom2(link = "log"),
  spatial = "off"
)

SharkDensitySmooth <- sdmTMB(
  formula = TotalOver150 ~ s(Depth) + as.factor(Year), 
  data = BLLCoordSelect,
  mesh = mesh, 
  family = nbinom2(link = "log"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "IID"
)

SharkDensityLinear <- sdmTMB(
  formula = TotalOver150 ~ Depth + as.factor(Year), 
  data = BLLCoordSelect,
  mesh = mesh, 
  family = nbinom2(link = "log"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "IID"
)

### AIC table to compare depth 
aic_tableDensity <- AIC(SharkDensityNull, SharkDensitySmooth, SharkDensityLinear)
aic_tableDensity$delta_aic <- round(aic_tableDensity$AIC-min(aic_tableDensity$AIC),2)
aic_tableDensity <- aic_tableDensity %>% arrange(delta_aic)
aic_tableDensity
  # Model favors depth as a linear predictor (by 1.39)
  # use as final model
  # not considering different random effects 
    # IID and Year predictor variable will not constrain year to year abundance estimates


### Probability of presence  models 
## Using hierarchical approach, first determining if linear or smooth depth is favored
# Null
Null <- sdmTMB(
  formula = Over150Present ~ 1,
  data = BLLCoordSelect,
  mesh = mesh, 
  family = binomial(link = "logit") 
)

# Depth Smooth, no spatial random effect
SmoothDepth <- sdmTMB(
  formula = Over150Present ~ s(Depth),
  data = BLLCoordSelect,
  mesh = mesh,  
  family = binomial(link = "logit") 
)

# Depth Linear, no spatial random effect
LinearDepth <- sdmTMB(
  formula = Over150Present ~ Depth,
  data = BLLCoordSelect,
  mesh = mesh,  
  family = binomial(link = "logit") 
)

# AIC to determine if depth should be linear or smooth
aic_tableDepth <- AIC(Null, SmoothDepth, LinearDepth)
aic_tableDepth$delta_aic <- round(aic_tableDepth$AIC-min(aic_tableDepth$AIC),2)
aic_tableDepth <- aic_tableDepth %>% arrange(delta_aic)
aic_tableDepth
  # linear depth is favored, include as linear effect. 

### Next decision: Spatiotemporal effect
RW <- sdmTMB(
  formula = Over150Present ~ Depth,
  data = BLLCoordSelect,
  mesh = mesh, 
  family = binomial(link = "logit"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "RW"
)
sanity(RW) 

AR1 <- sdmTMB(
  formula = Over150Present ~ Depth,
  data = BLLCoordSelect,
  mesh = mesh, 
  family = binomial(link = "logit"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "AR1"
)
sanity(AR1) 

IID <- sdmTMB(
  formula = Over150Present ~ Depth,
  data = BLLCoordSelect,
  mesh = mesh, 
  family = binomial(link = "logit"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "IID"
)


aic_tableST <- AIC(RW, AR1, IID)
aic_tableST$delta_aic <- round(aic_tableST$AIC-min(aic_tableST$AIC),2)
aic_tableST <- aic_tableST %>% arrange(delta_aic)
aic_tableST

aic_tableTotal <- AIC(RW, AR1, IID, LinearDepth, SmoothDepth, Null)
aic_tableTotal$delta_aic <- round(aic_tableTotal$AIC-min(aic_tableTotal$AIC),2)
aic_tableTotal <- aic_tableTotal %>% arrange(delta_aic)
aic_tableTotal

### Rw is favored


