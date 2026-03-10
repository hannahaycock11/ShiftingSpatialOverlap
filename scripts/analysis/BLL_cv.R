#### cross validation

## Load all necessary libraries 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sdmTMB)
library(future)


# Set working directory 
setwd("~/Research Materials/DepredationPotential_Overlay/DepredationCode")

## read in BLL data 
BLL <- read_csv("data/processed/BLLPresence.csv")

## Filter BLL data to only coast of Florida
# Lon  Min: -87.94  Max: -80.00 Lat Min: 23.73 Max: 30.46
BLLCoord <- BLL %>% filter(between(Latitude, 23, 31), between(Longitude, -87, -79.00),
                           YEAR >= 2000) %>% 
  rename(Year = YEAR)

## Identify East and West coast and filter out the points due to the shape of Florida. 
BLLCoord <- BLLCoord %>% 
  mutate(Coast = ifelse((Longitude < -80.4 & Latitude < 28) | 
                          (Longitude <= -82 & Latitude >= 28),"W","E")) %>%
  filter( Coast == "W")

## Plot points to check filter
ggplot(BLLCoord, aes(x=Longitude, y=Latitude, color= Coast))+
  geom_point()

## Remove depths that are NA
BLLCoord <- filter(BLLCoord, !is.na(Depth))
  # Make depth negative to match bathymetry data 
BLLCoord$Depth <- BLLCoord$Depth * -1
summary(BLLCoord)


## After model selection
  ## Add UTM points
BLLCoord <-sdmTMB::add_utm_columns(BLLCoord,ll_names=c("Longitude","Latitude"),units="km", 
                                   utm_crs = 32617)
  # UTM Zone 17N

## Create a mesh object
BLLmesh <- make_mesh(BLLCoord, xy_cols = c("X", "Y"), cutoff = 20)
BLLmesh$mesh$n #263
plot(BLLmesh)


### random k folds
library(future)
plan(multisession, workers = 2)


#### Depredating Shark Abundance Model ####

### first, training model (one used in analysis)
SharkDensity_training <- sdmTMB(
  formula = TotalOver150 ~ Depth + as.factor(Year), 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = nbinom2(link = "log"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "IID"
)

# model performance 
logLik(SharkDensity_training) 
# -1346.045  
train_pred_resp <- exp(predict(SharkDensity_training)$est)

train_rmse <- sqrt(mean((BLLCoord$TotalOver150 - train_pred_resp)^2))
train_rmse
# 0.8762951
train_mae <- mean(abs(BLLCoord$TotalOver150 - train_pred_resp))
train_mae
# 0.5026338


## Random k folds, k = 4, 8 
SharkDensity_cv <- sdmTMB_cv(
  formula = TotalOver150 ~ Depth + as.factor(Year), 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = nbinom2(link = "log"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "IID",
  k_folds = 8
)

SharkDensity_cv$fold_loglik # fold log-likelihood

SharkDensity_cv$sum_loglik # total log-likelihood
# -1378.747, -1373.525

# RMSE across entire dataset:
sqrt(mean((SharkDensity_cv$data$TotalOver150- SharkDensity_cv$data$cv_predicted)^2))
# 0.995756, 0.9970072
# MAE across entire dataset:
mean(abs(SharkDensity_cv$data$TotalOver150- SharkDensity_cv$data$cv_predicted))
#  0.5653651, 0.5638832

### Spatial k folds ####
# folds by cluster kmeans = 10, 20, 50
BLLCoord$clust <- kmeans(BLLCoord[, c("X", "Y")], 120)$cluster

ggplot(BLLCoord, aes(Longitude, Latitude, color = clust)) +
  geom_point() +
  facet_wrap(~clust)

ggplot(BLLCoord, aes(Longitude, Latitude, color = clust)) +
  geom_point()

SharkDensity_spatial <- sdmTMB_cv(
  formula = TotalOver150 ~ Depth + as.factor(Year), 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = nbinom2(link = "log"),
  spatial = "off",
  time = "Year",
  spatiotemporal = "off",
  fold_ids = clust
)

SharkDensity_spatial$fold_loglik 
SharkDensity_spatial$sum_loglik
# -1411.969, -1408.238, -1409.653

# RMSE across entire dataset:
sqrt(mean((SharkDensity_spatial$data$TotalOver150- SharkDensity_spatial$data$cv_predicted)^2))
# 1.03434,  1.026975, 1.026117
# MAE across entire dataset:
mean(abs(SharkDensity_spatial$data$TotalOver150- SharkDensity_spatial$data$cv_predicted))
# 0.6015424, 0.5954956, 0.5951232



### Year k folds #### 
SharkDensity_Spatialtraining <- sdmTMB(
  formula = TotalOver150 ~ Depth, 
  data = BLLCoord,
  mesh = BLLmesh, 
  time= "Year",
  family = nbinom2(link = "log")
)

# model performance 
logLik(SharkDensity_Spatialtraining) 
# -1397.043   
train_pred_resp <- exp(predict(SharkDensity_Spatialtraining)$est)

train_rmse <- sqrt(mean((BLLCoord$TotalOver150 - train_pred_resp)^2))
train_rmse
# 0.9393269
train_mae <- mean(abs(BLLCoord$TotalOver150 - train_pred_resp))
train_mae
# 0.5409381



### Year k folds #### 
SharkDensity_Temptraining <- sdmTMB(
  formula = TotalOver150 ~ Depth, 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = nbinom2(link = "log")
  )

# model performance 
logLik(SharkDensity_Temptraining) 
# -1460.521  
train_pred_resp <- exp(predict(SharkDensity_Temptraining)$est)

train_rmse <- sqrt(mean((BLLCoord$TotalOver150 - train_pred_resp)^2))
train_rmse
# 1.004424
train_mae <- mean(abs(BLLCoord$TotalOver150 - train_pred_resp))
train_mae
# 0.6043167

clust <- as.numeric(as.factor(BLLCoord$Year))
SharkDensity_time <- sdmTMB_cv(
  formula = TotalOver150 ~ Depth, 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = nbinom2(link = "log"),
  fold_ids = clust
)

SharkDensity_time$fold_loglik 
SharkDensity_time$sum_loglik
# -1467.09

# RMSE across entire dataset:
sqrt(mean((SharkDensity_time$data$TotalOver150- SharkDensity_time$data$cv_predicted)^2))
  # 1.046071
# MAE across entire dataset:
mean(abs(SharkDensity_time$data$TotalOver150- SharkDensity_time$data$cv_predicted))
  # 0.6328816

#### Shark Presence Model ####

### start with training model
SharkPresence_training <- sdmTMB(
  formula = Over150Present ~ Depth, 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = binomial(link = "logit"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "RW"
)

# model performance 
logLik(SharkPresence_training) 
  #   -793.9336
lp <- predict(SharkPresence_training)$est
# Convert to probability using inverse logit
train_pred_resp <- 1 / (1 + exp(-lp))

train_rmse <- sqrt(mean((BLLCoord$Over150Present - train_pred_resp)^2))
train_rmse
  #  0.3669089
train_mae <- mean(abs(BLLCoord$Over150Present - train_pred_resp))
train_mae
  #  0.2831805

# random folds k = 4, 8 
SharkPresence_cv <- sdmTMB_cv(
  formula = Over150Present ~ Depth, 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = binomial(link = "logit"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "RW",
  k_folds = 4
)

SharkPresence_cv$fold_loglik # fold log-likelihood
SharkPresence_cv$sum_loglik # total log-likelihood
  # -799.4414,-784.0106

# RMSE across entire dataset:
sqrt(mean((SharkPresence_cv$data$Over150Present- SharkPresence_cv$data$cv_predicted)^2))
  # 0.39236, 0.3879907
# MAE across entire dataset:
mean(abs(SharkPresence_cv$data$Over150Present- SharkPresence_cv$data$cv_predicted))
  # 0.3036977, 0.2993451

#### spatial Presence ####

SharkPresence_SpatialTraining <- sdmTMB(
  formula = Over150Present ~ Depth, 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = binomial(link = "logit"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "off"
)

# model performance 
logLik(SharkPresence_SpatialTraining) 
# -859.9367 
lp <- predict(SharkPresence_SpatialTraining)$est
# Convert to probability using inverse logit
train_pred_resp <- 1 / (1 + exp(-lp))

train_rmse <- sqrt(mean((BLLCoord$Over150Present - train_pred_resp)^2))
train_rmse
# 0.4015221
train_mae <- mean(abs(BLLCoord$Over150Present - train_pred_resp))
train_mae
#   0.3281001

clust <- kmeans(BLLCoord[, c("X", "Y")], 50)$cluster
# kmeans = 10, 20, 50 

SharkPresence_Spatial <- sdmTMB_cv(
  formula = Over150Present ~ Depth, 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = binomial(link = "logit"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "off",
  fold_ids = clust
)

SharkPresence_Spatial$fold_loglik 
SharkPresence_Spatial$sum_loglik
#-869.5747, -870.2004

# RMSE across entire dataset:
sqrt(mean((SharkPresence_Spatial$data$Over150Present- SharkPresence_Spatial$data$cv_predicted)^2))
#   0.4133958, 0.4143246
# MAE across entire dataset:
mean(abs(SharkPresence_Spatial$data$Over150Present- SharkPresence_Spatial$data$cv_predicted))
#  0.3363918, 0.3401492 




#### Time Presence ####
#fold = Year 
clust <- as.numeric(as.factor(BLLCoord$Year))

SharkPresence_time <- sdmTMB_cv(
  formula = Over150Present ~ Depth, 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = binomial(link = "logit"),
  spatial = "on",
  fold_ids = clust
)

SharkPresence_time$fold_loglik 
SharkPresence_time$sum_loglik
#  -866.944
# RMSE across entire dataset:
sqrt(mean((SharkPresence_time$data$Over150Present- SharkPresence_time$data$cv_predicted)^2))
# 0.4132313
# MAE across entire dataset:
mean(abs(SharkPresence_time$data$Over150Present- SharkPresence_time$data$cv_predicted))
#  0.3383248

SharkPresence_timeNull <- sdmTMB(
  formula = Over150Present ~ Depth, 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = binomial(link = "logit"),
  spatial = "on"
)

# model performance 
logLik(SharkPresence_timeNull) 
#   -859.9367
lp <- predict(SharkPresence_timeNull)$est
# Convert to probability using inverse logit
train_pred_resp <- 1 / (1 + exp(-lp))

train_rmse <- sqrt(mean((BLLCoord$Over150Present - train_pred_resp)^2))
train_rmse 
# 0.4015221
train_mae <- mean(abs(BLLCoord$Over150Present - train_pred_resp))
train_mae
# 0.3281001

