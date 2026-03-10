## Load all necessary libraries 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sdmTMB)
library(future)


# Set working directory 
setwd("~/Research Materials/DepredationPotential_Overlay/DepredationCode")

## Read in weighted ASOS data (see ASOS_Model.R)
ASOSCoord_25 <- read.csv("data/processed/ASOSCoord2_5.csv")
ASOSCoord_25$X.1 <- NULL

# small values at presence locations
ASOSCoord_25$wt <- 1e-6
# pseudo-absences: area per quadrature point
n_zeros <- length(which(ASOSCoord_25$Present == 0))
ASOSCoord_25$wt <- ifelse(ASOSCoord_25$Present == 1,
                          1e-6, 269679.1 / n_zeros
)

ASOSMesh <- make_mesh(ASOSCoord_25, xy_cols = c("X", "Y"), cutoff = 30)
  
trainingASOS <- sdmTMB(
  formula = Present / wt ~ s(Depth),
  data = ASOSCoord_25,
  mesh = ASOSMesh,
  family = poisson(link = "log"),
  spatial = "on",
  weights = ASOSCoord_25$wt,
  time = "Year",
  spatiotemporal = "RW"
)

# model performance 
logLik(trainingASOS) 
# -193709.3 

train_pred_resp <- exp(predict(trainingASOS)$est) * ASOSCoord_25$wt
train_rmse_NF <- sqrt(mean((ASOSCoord_25$Present - train_pred_resp)^2))
train_rmse_NF
# 0.4439637 
train_mae_NF <- mean(abs(ASOSCoord_25$Present - train_pred_resp))
train_mae_NF
#  0.1510937


#----------------------------------------
### Run the loops 
# random k-fold validation
set.seed(123)
# Assign each row to one of 4 folds
ASOSCoord_25$fold <- sample(1:8, nrow(ASOSCoord_25), replace = TRUE)

results4 <- map_dfr(1:8, function(k) {
  
  # ---- Train = all rows not in fold k ----
  train <- ASOSCoord_25[ASOSCoord_25$fold != k, ]
  test  <- ASOSCoord_25[ASOSCoord_25$fold == k, ]
  
  # ---- Build mesh from training data only ----
  mesh_k <- make_mesh(train, xy_cols = c("X", "Y"), cutoff = 30)
  
  # ---- Fit model ----
  fit_k <- sdmTMB(
    formula = Present / wt ~ s(Depth),
    data = train,
    mesh = mesh_k,
    family = poisson(link = "log"),
    spatial = "on",
    weights = train$wt,
    time = "Year",
    spatiotemporal = "RW"
  )
  
  # ---- Predict on TEST (out-of-fold) ----
  pred_test <- predict(fit_k, newdata = test)
  test_pred_count <- exp(pred_test$est) *test$wt
  
  test %>% 
    summarize(
      Fold = k,
      Test_RMSE = sqrt(mean((Present - test_pred_count)^2)),
      Test_MAE  = mean(abs(Present - test_pred_count)),
      Test_LL   = sum(dpois(Present, lambda = test_pred_count, log = TRUE)), 
      test_mse = sum((Present - test_pred_count)^2),
      test_maeNew = sum(abs(Present - test_pred_count)),
      n = nrow(test))
})


results_summary8 <- results4 %>% 
  summarize(
    Test_RMSE_mean  = mean(Test_RMSE),
    Test_MAE_mean   = mean(Test_MAE),
    Test_LL_mean    = mean(Test_LL),
    Test_RMSE_sum  = sum(test_mse),
    Test_MAE_sum   = sum(test_maeNew),
    Test_LL_sum    = sum(Test_LL),
    ntotal = sum(n)
  ) %>% 
  mutate(
    RMSE = sqrt(sum(Test_RMSE_sum)/ntotal),
    MAE = sum(Test_MAE_sum)/ntotal
  )
results_summary4
results_summary8



### Year folds
timeblocks <- sort(unique(ASOSCoord_25$Year))
resultsTB <- map_df(timeblocks, function(tb) {
  
  message("Running fold leaving out Year = ", tb)
  
  # ---- Split train / test ----
  train <- ASOSCoord_25 %>% 
    filter(Year != tb) 
  
  test  <- ASOSCoord_25 %>% 
    filter(Year == tb) 
  
  # ----------------------
  # MESH ON TRAINING DATA ONLY
  # ----------------------
  mesh_tb <- make_mesh(train, xy_cols = c("X", "Y"), cutoff = 30)
  
  # ----------------------
  # FIT MODEL ON TRAINING SET
  # ----------------------
  fit_tb <- sdmTMB(
    formula = Present / wt ~ s(Depth),
    data = train,
    mesh = mesh_tb,
    family = poisson(link = "log"),
    spatial = "on",
    weights = train$wt
  )
  
  # predict on held-out Year
  test_pred <- predict(fit_tb, newdata = test)
  test_pred_rate  <- exp(test_pred$est)
  test_pred_count <- test_pred_rate * test$wt
  
  # ---- Predict on TEST (out-of-fold) ----
  test %>% 
    summarize(
      Fold = tb,
      Test_RMSE = sqrt(mean((Present - test_pred_count)^2)),
      Test_MAE  = mean(abs(Present - test_pred_count)),
      Test_LL   = sum(dpois(Present, lambda = test_pred_count, log = TRUE)), 
      test_mse = sum((Present - test_pred_count)^2),
      test_maeNew = sum(abs(Present - test_pred_count)),
      n = nrow(test)
    )
})

results_summary_time <- resultsTB %>% 
  summarize(
    Test_RMSE_mean  = mean(Test_RMSE),
    Test_MAE_mean   = mean(Test_MAE),
    Test_LL_mean    = mean(Test_LL),
    Test_RMSE_sum  = sum(test_mse),
    Test_MAE_sum   = sum(test_maeNew),
    Test_LL_sum    = sum(Test_LL),
    ntotal = sum(n)
  ) %>% 
  mutate(
    RMSE = sqrt(sum(Test_RMSE_sum)/ntotal),
    MAE = sum(Test_MAE_sum)/ntotal
  )
results_summary_time

trainingASOS_Time <- sdmTMB(
  formula = Present / wt ~ s(Depth),
  data = ASOSCoord_25,
  mesh = ASOSMesh,
  family = poisson(link = "log"),
  spatial = "on",
  weights = ASOSCoord_25$wt
)

# model performance 
logLik(trainingASOS_Time) 
#  -194551.8 
train_pred_resp <- exp(predict(trainingASOS_Time)$est) * ASOSCoord_25$wt
train_rmse_NF <- sqrt(mean((ASOSCoord_25$Present - train_pred_resp)^2))
train_rmse_NF
# 0.4271811 
train_mae_NF <- mean(abs(ASOSCoord_25$Present - train_pred_resp))
train_mae_NF
#  0.151099


# Spatial Loops 
remove(ASOSCoord_25)
ASOSCoord_25 <- read.csv("data/processed/ASOSCoord2_5.csv")
ASOSCoord_25$X.1 <- NULL

# small values at presence locations
ASOSCoord_25$wt <- 1e-6

# pseudo-absences: area per quadrature point
n_zeros <- length(which(ASOSCoord_25$Present == 0))

ASOSCoord_25$wt <- ifelse(
  ASOSCoord_25$Present == 1,
  1e-6,
  269679.1 / n_zeros
)

### Calculating grid size needed for 500 grid cells 
sqrt(269679.1/500)
  ## 23.93064 
  ## use 52 as ceiling and floor 

ASOSCoord_25 <- ASOSCoord_25 %>% 
  mutate(
    grid_x = floor(X / 23),
    grid_y = floor(Y / 23),
    grid_id = paste(grid_x, grid_y, sep = "_")
  )


cluster_order <- ASOSCoord_25 %>% 
  group_by(grid_id) %>% 
  summarise(Total = sum(Present), .groups = "drop") %>% 
  arrange(desc(Total)) %>% 
  mutate(cluster_rank = row_number())

cluster_order <- cluster_order %>%
  mutate(
    fold_10 = ((cluster_rank - 1) %% 10) + 1
  )

## Create random grid selection
set.seed(123)
set.seed(456)
set.seed(789)


# create dataframe of unique grids
grid_folds <- ASOSCoord_25 %>%
  distinct(grid_id) %>%
  mutate(fold_rand10 = sample(1:10, n(), replace = TRUE))

ASOSCoord_25 <- ASOSCoord_25 %>%
  left_join(cluster_order, by = "grid_id")


# join back to original data
ASOSCoord_25Rand3 <- ASOSCoord_25 %>%
  left_join(grid_folds, by = "grid_id")

ggplot(ASOSCoord_25Rand3, aes(Longitude, Latitude, color = fold_rand10)) +
  geom_point()
# ------------------------
# Step 2: loop through folds
# ------------------------
folds <- sort(unique(ASOSCoord_25Rand3$fold_rand10))

results_spatial <- map_df(folds, function(f) {
  
  message("Running spatial fold = ", f)
  
  # ---- Split train / test ----
  train <- ASOSCoord_25Rand3 %>% filter(fold_rand10 != f)
  test  <- ASOSCoord_25Rand3 %>% filter(fold_rand10 == f)
  
  # ---- Build mesh on training data only ----
  mesh_f <- make_mesh(train, xy_cols = c("X", "Y"), cutoff = 30)
  # ---- Fit model ----
  fit_f <- sdmTMB(
    formula = Present / wt ~ s(Depth),
    data = train,
    mesh = mesh_f,
    family = poisson(link = "log"),
    spatial = "off",
    time = "Year",
    weights = train$wt
  )
  
  # ---- Predict on TEST (out-of-fold) ----
  test_pred <- predict(fit_f, newdata = test)
  test_pred_rate  <- exp(test_pred$est)
  test_pred_count <- test_pred_rate * test$wt
  
  test %>% 
    summarize(
      Fold = f,
      Test_RMSE = sqrt(mean((Present - test_pred_count)^2)),
      Test_MAE  = mean(abs(Present - test_pred_count)),
      Test_LL   = sum(dpois(Present, lambda = test_pred_count, log = TRUE)), 
      test_mse = sum((Present - test_pred_count)^2),
      test_maeNew = sum(abs(Present - test_pred_count)),
      n = nrow(test)
    )
})


Spatial_Summary10YR_456 <- results_spatial %>% 
  summarize(
    Test_RMSE_mean  = mean(Test_RMSE),
    Test_MAE_mean   = mean(Test_MAE),
    Test_LL_mean    = mean(Test_LL),
    Test_RMSE_sum   = sum(test_mse),
    Test_MAE_sum    = sum(test_maeNew),
    Test_LL_sum     = sum(Test_LL),
    ntotal          = sum(n)
  ) %>% 
  mutate(
    RMSE = sqrt(Test_RMSE_sum / ntotal),
    MAE  = Test_MAE_sum / ntotal
  )

Spatial_Summary10YR


Spatial_Summary10Y


ggplot(ASOSCoord_25, aes(Longitude, Latitude, color = factor(fold_5))) +
  geom_point()

Clust10 <- ggplot(ASOSCoord_25, aes(Longitude, Latitude, color = fold_10)) +
  geom_point()
  # 10 per fold

Clust2
Clust3
Clust5
Clust6


### Spatial Training
Spatial_Training <- sdmTMB(
  formula = Present / wt ~ s(Depth),
  data = ASOSCoord_25,
  mesh = ASOSMesh,
  family = poisson(link = "log"),
  spatial = "off",
  time = "Year",
  weights = ASOSCoord_25$wt
)

# model performance 
logLik(Spatial_Training) 
# -214642.6  (Y: -194861)
 
train_pred_resp <- exp(predict(Spatial_Training)$est) * ASOSCoord_25$wt
train_rmse_NF <- sqrt(mean((ASOSCoord_25$Present - train_pred_resp)^2))
train_rmse_NF
# 0.2925064   (Y: 0.4456516)
train_mae_NF <- mean(abs(ASOSCoord_25$Present - train_pred_resp))
train_mae_NF
#  0.1510764  (Y: 0.1509799)

