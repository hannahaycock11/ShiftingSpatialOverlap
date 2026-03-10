## Code used Depredating shark abundnace and presence models ##

# Set working directory 
setwd("~/Research Materials/DepredationPotential_Overlay/DepredationCode")

# Function to save all plots to BLL figures folder# 
save_plot <- function(obj, filename, folder = "figures/BLL", type = "ggplot", device = "png", width = 800, height = 600, ...) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  filepath <- file.path(folder, filename)
  
  if (type == "ggplot") {
    ggsave(filepath, plot = obj, width = width / 100, height = height / 100, dpi = 300, ...)
  } else if (type == "base") {
    do.call(device, c(list(filename = filepath, width = width, height = height), list(...)))
    eval(obj)
    dev.off()
  } else {
    stop("Unknown plot type. Use 'ggplot' or 'base'.")
  }
}


## Load all necessary libraries 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sdmTMB)
library(spdep)
library(pROC) 
library(av)
library(mgcv)
library(marmap)
library(gridExtra)
library(gridGraphics)
library(patchwork)
theme_set(theme_bw())

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
summary(BLLCoord$Depth)
BLLCoord <- filter(BLLCoord, !is.na(Depth))
# Make depth negative to match bathymetry data 
BLLCoord$Depth <- BLLCoord$Depth * -1
summary(BLLCoord)

## Create .csv file used for BLL model selection
  # See Rscript BLL_ModelSelection.R for process
  # Favored models will be used to create final models in this RScript
write.csv(BLLCoord,"data/processed/BLL_model_selection.csv")

## After model selection, continue to creating model here 
## Add UTM points
BLLCoord <-sdmTMB::add_utm_columns(BLLCoord,ll_names=c("Longitude","Latitude"),units="km", 
                                   utm_crs = 32617)
  # UTM Zone 17N

## Create a mesh object
BLLmesh <- make_mesh(BLLCoord, xy_cols = c("X", "Y"), cutoff = 20)
BLLmesh$mesh$n #263
plot(BLLmesh)

nrow(BLLCoord)

#### Depredating Shark Abundance Model 
SharkDensity <- sdmTMB(
  formula = TotalOver150 ~ Depth + as.factor(Year), 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = nbinom2(link = "log"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "IID"
)

sanity(SharkDensity) # passes sanity check 
summary(SharkDensity)

## summary of fixed and random effects
tidy(SharkDensity, conf.int = TRUE)
tidy(SharkDensity, "ran_pars", confint = TRUE)

print(n=30, tidy(SharkDensity, conf.int = TRUE))
logLik(SharkDensity) # total log-likelihood
# -1362.237

### Model Diagnostics 
## Normal Q-Q Plot 
BLLCoord$resids <- residuals(SharkDensity) 
qqnorm(BLLCoord$resids) 
qqline(BLLCoord$resids) 
  # Save figure 
save_plot(quote({
  qqnorm(BLLCoord$resids)
  qqline(BLLCoord$resids)
}), "Figure2A.png", type = "base", width = 800, height = 600)
## simulate data 
sharksim <- simulate(SharkDensity, nsim = 500, type = "mle-mvn")
  # Compare fraction of zeroes 
sum(BLLCoord$TotalOver150 == 0) / length(BLLCoord$TotalOver150)
  # 0.735226
sum(sharksim == 0)/length(sharksim)
  # 0.7222028



### Model Performance
predictions <- predict(SharkDensity, type = "response")$est
## Compute ROC curve
roc_curve <- roc(BLLCoord$TotalOver150, predictions)
## Plot ROC curve
plot(roc_curve, main = "ROC Curve", col = "blue")
  # Save figure 
save_plot(quote({
  plot(roc_curve, main = "ROC Curve", col = "blue")
}), "Figure2c.png", type = "base", width = 800, height = 600)
## Print AUC
auc(roc_curve) 
  # 0.8044



qq_plot <- ggplot(BLLCoord, aes(sample = resids)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "A) Residual QQ Plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_bw()


roc_df <- data.frame(
  specificity = roc_curve$specificities,
  sensitivity = roc_curve$sensitivities
)

roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line() +
  geom_abline(linetype = "dashed") +
  annotate("text", x = 0.7, y = 0.2,
           label = paste("AUC =", round(auc(roc_curve), 3))) +
  labs(title = "B) ROC Curve",
       x = "Specifity",
       y = "Sensitivity") +
  theme_bw()
combined_plot <- qq_plot + roc_plot 
combined_plot


### Predicting Index of Abundance 
## Make prediction grid
  # Define the area: lat and lon range
lon_range <- c(-87.5, -80.5)  # longitude range
lat_range <- c(24, 31)   
  # Get bathymetry data
bathy_data <- getNOAA.bathy(lon1 = lon_range[1], lon2 = lon_range[2], 
                            lat1 = lat_range[1], lat2 = lat_range[2], 
                            resolution = 3) # ~5.55 km
bathy_new <- as.xyz(bathy_data)
names(bathy_new) <- c("Longitude","Latitude", "Depth")
  # Assign UTM values
bathy_new <-sdmTMB::add_utm_columns(bathy_new,ll_names=c("Longitude","Latitude"),units="km", 
                                    utm_crs = 32617)
  # Filter lat, lon, and depth to match BLLCoord data
predgrid <- bathy_new
predgrid <- predgrid %>% 
  mutate(Coast = ifelse((Longitude < -80.4 & Latitude < 28) | 
                          (Longitude <= -82 & Latitude >= 28),"W","E")) %>%
  filter( Coast == "W") %>% 
  filter( between(Depth, -376, 0)) %>% 
  filter(!(Latitude > 26.5 & Longitude > -81.5)) %>% 
  filter(!(Latitude < 25 & Longitude < -86))
predgrid$Coast<-NULL
## Plot to check grid
plot(predgrid$Longitude, predgrid$Latitude) 
# Replicate TimeBlock
predgrid <- replicate_df(predgrid, "Year", unique(BLLCoord$Year))

## Predict with sdmTMBobjectreturn = FALSE (default) to estimate per region and write csv file 
DensityYB2 <- predict(SharkDensity, newdata = predgrid)
summary(DensityYB2)
DensityYB2 <- DensityYB2 %>%
  mutate(MeanNumShark = exp(est))

## Model prediction
DensityYB <- predict(SharkDensity, newdata = predgrid, return_tmb_object = TRUE)

## Predict Index of Abundance
indYB <- get_index(DensityYB, bias_correct = TRUE)

time_labels <- c(
  "1" = "2000–2002",
  "2" = "2003–2005",
  "3" = "2006–2008",
  "4" = "2009–2011",
  "5" = "2012–2014",
  "6" = "2015–2017",
  "7" = "2018–2020",
  "8" = "2021–2023"
)

DepredatingIndex <- ggplot(indYB, aes(x = Year, y = est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#69b3a2", alpha = 0.3) +
  geom_line(color = "#1f78b4", linewidth = 1.2) +
  scale_x_continuous(
    #breaks = as.numeric(names(time_labels)),
    #labels = time_labels,
    expand = c(0, 0) 
  ) +
  scale_y_continuous(
    expand = c(0, 0) 
  ) +
  labs(
    x = "Year",
    y = "Number of Depredating Sharks"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) 
DepredatingIndex
save_plot(DepredatingIndex, "Figure3.png", type = "ggplot", width = 800, height = 600)

ggplot() + 
  geom_tile(data = DensityYB2, 
            aes(x = Longitude, y = Latitude, fill = MeanNumShark)) +  
  geom_sf(data = florida_shape, fill = "lightgray", color = "black") + 
  labs(
    x = "Longitude", 
    y = "Latitude") +
  scale_fill_viridis_c() +  
  theme_minimal() + 
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_sf(expand = FALSE) +  # Keep correct spatial aspect
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))


### Depredating Shark Presence Model 
SharkPresence <- sdmTMB(
  formula = Over150Present ~ Depth, 
  data = BLLCoord,
  mesh = BLLmesh, 
  family = binomial(link = "logit"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "RW"
)

sanity(SharkPresence) # passes sanity check
summary(SharkPresence)

## Summary of fixed and random effects
tidy(SharkPresence, conf.int = TRUE)
tidy(SharkPresence, "ran_pars", confint = TRUE)

### Model Diagnostics 
## Normal Q-Q Plot 
BLLCoord$resids <- residuals(SharkPresence) 
qqnorm(BLLCoord$resids) 
qqline(BLLCoord$resids) 
# Save figure 

## simulate data 
sharksim <- simulate(SharkPresence, nsim = 500, type = "mle-mvn")
# Compare fraction of zeroes 
sum(BLLCoord$Over150Present == 0) / length(BLLCoord$Over150Present)
# 0.735226
sum(sharksim == 0)/length(sharksim)
# 0.7127045

PresenceResid <- ggplot(BLLCoord, aes(x = X, y = Y, color = resids)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradient2(
    name = "Residuals"
  ) +
  facet_wrap(~Year, ncol = 6) +
  coord_fixed() +
  labs(
    x = "X (UTM)",
    y = "Y (UTM)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(size = 10),
    strip.text = element_text( size = 10),
    legend.title = element_text(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank()
  )
PresenceResid

save_plot(PresenceResid, "FigureS10.png", type = "ggplot", width = 1000, height = 1000)
### Model Performance
predictions <- predict(SharkPresence, type = "response")$est
## Compute ROC curve
roc_curve <- roc(BLLCoord$Over150Present, predictions)
## Plot ROC curve
plot(roc_curve, main = "ROC Curve", col = "blue")
# Save figure 
auc(roc_curve) 
  # 0.8531

qq_plot <- ggplot(BLLCoord, aes(sample = resids)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "A) Residual QQ Plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_bw()


roc_df <- data.frame(
  specificity = roc_curve$specificities,
  sensitivity = roc_curve$sensitivities
)

roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line() +
  geom_abline(linetype = "dashed") +
  annotate("text", x = 0.7, y = 0.2,
           label = paste("AUC =", round(auc(roc_curve), 3))) +
  labs(title = "B) ROC Curve",
       x = "Specifity",
       y = "Sensitivity") +
  theme_bw()
combined_plot <- qq_plot + roc_plot 
combined_plot


### Predict Probability of Presence using same prediction grid
SharkPredPres <- predict(SharkPresence, newdata = predgrid)
### logit transform to get probability of presence 
SharkPredPres <- SharkPredPres %>% 
  mutate(ProbabilityPresence = exp(est)/(1+exp(est)))

  # write .csv file of output to use in depredation risk
write.csv(SharkPredPres, "data/output/BLL_shark_presence_output.csv")

## Assess random and fixed effects
Fixed <- ggplot() +
  geom_tile(data = SharkPredPres, 
            aes(X, Y, fill = exp(est_non_rf)), 
            width = 6, height = 6) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Shark Probability of Presence Model - Prediction with Fixed Effects Only")+
  coord_fixed()
save_plot(Fixed, "FigureS11a.png", type = "ggplot", width = 800, height = 600)

Random <- ggplot() +
  geom_tile(data = SharkPredPres, 
            aes(X, Y, fill = omega_s), 
            width = 6, height = 6) +
  scale_fill_gradient2() +
  ggtitle("Shark Probability of Presence Model - Spatial Random Effects Only")+
  coord_fixed()
save_plot(Random, "FigureS11b.png", type = "ggplot", width = 800, height = 600)

Spatiotemporal <- ggplot() +
  geom_tile(data = SharkPredPres, 
            aes(X, Y, fill = epsilon_st), 
            width = 6, height = 6) +
  scale_fill_gradient2() +
  facet_wrap(~Year, ncol = 6) + 
  ggtitle("Shark Probability of Presence Model - Spatiotemporal Random Effects Only")+
  coord_fixed()
save_plot(Spatiotemporal, "FigureS11c.png", type = "ggplot", width = 800, height = 600)

### Make Probability of presence plot
## Add Florida 
florida_shape <- st_read("scripts/analysis/Detailed_Florida_State_Boundary.shp")
# Reproject shapefile to UTM Zone 17N to match the bathymetry data
if (st_crs(florida_shape)$epsg != 4326) {
  florida_shape <- st_transform(florida_shape, crs = 4326)
}

### For visualizations sake - take the average every three years
shark_spatial_blocks <- SharkPredPres %>%
  mutate(timeblock = cut(Year, 
                         breaks = seq(2000, 2024, by = 3), 
                         right = FALSE, 
                         labels = c("2000-2002", "2003-2005", "2006-2008", 
                                    "2009-2011", "2012-2014", "2015-2017", 
                                    "2018-2020", "2021-2023"))) %>%
  # Group by time AND space (X and Y coordinates)
  group_by(timeblock, Latitude, Longitude) %>% 
  summarize(
    mean_prob = mean(ProbabilityPresence, na.rm = TRUE),
    .groups = "drop"
  )

SharkPresMap <- ggplot() + 
  geom_tile(data = shark_spatial_blocks, 
            aes(x = Longitude, y = Latitude, fill = mean_prob)) +  
  geom_sf(data = florida_shape, fill = "lightgray", color = "black") + 
  labs(
       x = "Longitude", 
       y = "Latitude") +
  scale_fill_viridis_c(limits = c(0, 1),option = "D", 
                       name = "Probability") +  
  theme_minimal() + 
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+
  facet_wrap(~timeblock, ncol=4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_sf(expand = FALSE) +  # Keep correct spatial aspect
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
SharkPresMap
save_plot(SharkPresMap, "Figure5.png", type = "ggplot", width = 800, height = 600)


### Regional Uncertainty predictions for presence ###
# Simulate from joint posterior
sim_Pres <- predict(SharkPresence, newdata = predgrid, nsim = 1000)

sim_df_Pres_Region <- as.data.frame(sim_Pres)

# Add Year and Region
sim_df_Pres_Region$Year <- predgrid$Year
predgrid <- predgrid  %>%
  mutate(Region = case_when(
    between(Latitude, 28.69, 31) & Longitude < -84.5 ~ "NW",
    Latitude > 28.69 & Longitude >= -84.5 ~ "BB",
    between(Latitude, 27, 28.69) ~ "TB",
    between(Latitude, 25.2, 27) ~ "SW",
    Latitude <= 25.2 ~ "KY",
    TRUE ~ "Other"  # fallback for rows that don't meet the condition
  )) 

sim_df_Pres_Region$Region <- predgrid$Region

sim_long_Pres_Region <- pivot_longer(sim_df_Pres_Region, 
                         cols = starts_with("V"), 
                         names_to = "sim", 
                         values_to = "value")

agg_sim_Pres_Region <- sim_long_Pres_Region %>%
  group_by(Year, Region, sim) %>%
  summarise(sum_pred = mean(plogis(value)), .groups = "drop")

region_summary_df_pres <- agg_sim_Pres_Region %>%
  group_by(Year, Region) %>%
  summarise(
    Presence = mean(sum_pred),
    Lower = quantile(sum_pred, 0.025),
    Upper = quantile(sum_pred, 0.975),
    SE_Total = sd(sum_pred)
  )

ggplot(region_summary_df_pres, aes(x = as.numeric(Year), y = Presence)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "lightblue", alpha = 0.4) +
  labs(
    title = "Change in Presence Over Time",
    x = "Time Block",
    y = "Total Predicted Abundance"
  ) +
  theme_minimal() +
  facet_wrap(~Region)

### Write csv to use in regional depredation risk predictions
write.csv(region_summary_df_pres,"data/output/PresRegion.csv")

## Making nice plot 
df_total2 <- region_summary_df_pres %>%
  group_by(Year) %>%
  summarise(
    Presence = mean(Presence, na.rm = TRUE),
    Lower          = mean(Lower, na.rm = TRUE),
    Upper          = mean(Upper, na.rm = TRUE),
    SE_Total       = mean(SE_Total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Region = "Total")

# bind back into main dataset
df_with_total2 <- bind_rows(region_summary_df_pres, df_total2)
summary(df_with_total2)

df_with_total2 <- df_with_total2 %>%
  mutate(Year = as.numeric(as.character(Year)))

region_labels <- c(
  "NW" = "Northwest",
  "BB" = "Big Bend",
  "TB" = "Tampa",
  "SW" = "Southwest",
  "KY" = "Keys",
  "Total" = "West Florida"
)

df_with_total2$RegionFull <- region_labels[df_with_total2$Region]

five_colors_okabe_ito <- c(
  "Northwest" = "#0072B2",     # Blue
  "Big Bend"  = "#009E73",     # Bluish Green
  "Tampa"     = "#E69F00",     # Orange
  "Southwest" = "#56B4E9",     # Sky Blue
  "Keys"      = "#CC79A7",      # Reddish Purple
  "West Florida"= "black"
)  

# Create lighter versions of the colors for the ribbons
five_colors_ribbon <- sapply(five_colors_okabe_ito, function(x) adjustcolor(x, alpha.f = 0.3))

RegionPresence <- ggplot(df_with_total2, aes(x = Year, y = Presence, color = RegionFull, group = RegionFull)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = RegionFull), color = NA) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = five_colors_okabe_ito, name = "Region") +
  scale_fill_manual(values = five_colors_ribbon, guide = "none") +
  scale_x_continuous(
    #breaks = as.numeric(names(time_labels)),
    #labels = time_labels,
    expand = c(0, 0) 
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Time Period",
    y = "Regional Probability of Presence"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, margin = margin(t = 2)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  ) +
  facet_wrap(~RegionFull)
RegionPresence
save_plot(RegionPresence, "Figure6.png", type = "ggplot", width = 800, height = 600) 



