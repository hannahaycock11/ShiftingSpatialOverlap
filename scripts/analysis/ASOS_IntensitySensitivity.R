### Hooked Fish Intensity Sensitivity Analysis and Model Selection

## Set working directory
setwd("~/Research Materials/DepredationPotential_Overlay/DepredationCode")

## Function to save all plots to ASOS figures folder
save_plot <- function(obj, filename, folder = "figures/ASOS", type = "ggplot", device = "png", width = 800, height = 600, ...) {
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
library(ggplot2)
library(sdmTMB)
library(marmap)
library(sf)
library(sp)
library(patchwork) 

## Read in weighted ASOS data (see ASOS_Model.R)
SensitivityData <- read.csv("data/processed/ASOS_SensitivySelectivity.csv")
nrow(SensitivityData)
summary(SensitivityData)

plot(SensitivityData$Longitude, SensitivityData$Latitude)

SensitivityData$X <- NULL
SensitivityData$Region <- NULL

## Add UTM columns
SensitivityData <- sdmTMB::add_utm_columns(SensitivityData,ll_names=c("Longitude","Latitude"),units="km", 
                                    utm_crs = 32617)

### Testing quadrature point grids 
## Define the area: lat and lon range
lon_range <- c(-87.5, -80.5)  # longitude range
lat_range <- c(24, 31)   
## Get bathymetry data
  # Sensitivity analysis tested eleven different resolutions 
  # To run the different resolutions, change resolution = to one of the following resolutions that were tested: 
    # 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8
bathy_data <- getNOAA.bathy(lon1 = lon_range[1], lon2 = lon_range[2], 
                            lat1 = lat_range[1], lat2 = lat_range[2], 
                            resolution = 2.5) 
bathy_new <- as.xyz(bathy_data)
names(bathy_new) <- c("Longitude","Latitude", "Depth")
## Filtering to match ASOS data lon, lat, and depth
bathy_new <-bathy_new %>% 
  mutate(Coast = ifelse((Longitude < -80.4 & Latitude < 28) | 
                          (Longitude <= -82 & Latitude >= 28),"W","E")) %>%
  filter( Coast == "W") %>% 
  filter( between(Depth, -376, 0)) %>% 
  filter(!(Latitude > 26.5 & Longitude > -81.5)) %>% 
  filter(!(Latitude < 25 & Longitude < -86))
bathy_new$Coast<-NULL
plot(bathy_new$Longitude, bathy_new$Latitude)

## Determine Area 
  # Need to to determine area in order to calculate intensity 
  # Since marmap uses lat/lon, need to convert to UTM points to get evenly spaced grids
  # Then use convex hull to turn grid points into a polygon, then calculate the area
UTMBathy <- sdmTMB::add_utm_columns(bathy_new,ll_names=c("Longitude","Latitude"),units="km", 
                        utm_crs = 32617)
  # Get convex hull indices
hpts <- chull(UTMBathy$X, UTMBathy$Y)
  # Close the polygon (repeat the first point at the end)
hpts <- c(hpts, hpts[1])
  # Extract convex hull coordinates
chull_coords <- UTMBathy[hpts, c("X","Y")]
  # Build polygon
chull.poly <- Polygon(chull_coords, hole = FALSE)
  # Extract area (in km²) 
chull_area <- chull.poly@area
chull_area # 269679.1

## Assign a quadrature point for each XY coordinate
QuadPoints <- UTMBathy
QuadPoints <- QuadPoints %>%
  mutate(Present = 0)
  # Replicate quadrature points for each year
QuadPoints <- replicate_df(QuadPoints, "Year", time_values = 2009:2023)

table(QuadPoints$Year)

nrow(QuadPoints)

summary(QuadPoints)
summary(SensitivityData)

##  Bind with ASOS data 
ASOS_Sensitivity <- rbind(SensitivityData, QuadPoints)

nrow(ASOS_Sensitivity)
summary(ASOS_Sensitivity)
plot(ASOS_Sensitivity$Longitude, ASOS_Sensitivity$Latitude)


write.csv(ASOS_Sensitivity, "data/processed/ASOSCoord2_5.csv")
### Downweighted Poisson Regression (DWPR)
    # An intercept only DWPR was run for each of the eleven different resolutions

## Make mesh
FHmesh <- make_mesh(ASOS_Sensitivity,xy_cols = c("X", "Y"), cutoff = 30)
FHmesh$mesh$n 
plot(FHmesh)

## Small values at presence locations
ASOS_Sensitivity$wt <- 1e-6

## Quadrature point weights
n_zeros <- length(which(ASOS_Sensitivity$Present == 0))

ASOS_Sensitivity$wt <- ifelse(ASOS_Sensitivity$Present == 1,
                     1e-6, chull_area / n_zeros
)

## Intercept only DWPR 
fit <- sdmTMB(
  Present / wt ~ 1,
  data = ASOS_Sensitivity,
  mesh = FHmesh,
  family = poisson(link = "log"),
  weights = ASOS_Sensitivity$wt
)

summary(fit)
sanity(fit)

### Compare Model Mapped Predictions
## Make prediction grid 
lon_range <- c(-87.5, -80.5) 
lat_range <- c(24, 31)   
bathy_data <- getNOAA.bathy(lon1 = lon_range[1], lon2 = lon_range[2], 
                            lat1 = lat_range[1], lat2 = lat_range[2], 
                            resolution = 3) # ~5.55 km
bathy_new <- as.xyz(bathy_data)
names(bathy_new) <- c("Longitude","Latitude", "Depth")
## Assign UTM values
bathy_new <-sdmTMB::add_utm_columns(bathy_new,ll_names=c("Longitude","Latitude"),units="km", 
                                    utm_crs = 32617)
predgrid <- bathy_new
predgrid <- predgrid %>% 
  mutate(Coast = ifelse((Longitude < -80.4 & Latitude < 28) | 
                          (Longitude <= -82 & Latitude >= 28),"W","E")) %>%
  filter( Coast == "W") %>% 
  filter( between(Depth, -376, 0)) %>% 
  filter(!(Latitude > 26.5 & Longitude > -81.5)) %>% 
  filter(!(Latitude < 25 & Longitude < -86))
predgrid$Coast<-NULL

FitPred <- predict(fit, newdata = predgrid)
FitPred <- FitPred %>%
  mutate(LogIntensity = est)

## Look at regional predictions
FitPred <- FitPred %>%
  mutate(Region = case_when(
    between(Latitude, 28.69, 31) & Longitude < -84.5 ~ "NW",
    Latitude > 28.69 & Longitude >= -84.5 ~ "BB",
    between(Latitude, 27, 28.69) ~ "TB",
    between(Latitude, 25.2, 27) ~ "SW",
    Latitude <= 25.2 ~ "KY",
    TRUE ~ "Other"  # fallback for rows that don't meet the condition
  )) 
RegionPred <- FitPred %>%
  group_by(Region) %>% 
  summarise(AverageIntensity = mean(est))
RegionPred 

## Look at total predictions
SumPred <- FitPred %>% 
  summarise(MeanIntensity = mean(est))
SumPred

### .csv summaries of model outputs for eleven resolutions assessed 
Output <- read.csv("data/output/ASOS_Sensitivity_Model_Outputs.csv")
Regional <- read.csv("data/output/ASOS_Sensitivity_Regional_Estimates.csv")
Total <- read.csv("data/output/ASOS_Sensitivity_Total_Estimates.csv")

## Look into different outputs
ggplot(Output, aes(Resolution, coef.est)) +
  geom_point()
ggplot(Output, aes(Resolution, coef.se)) +
  geom_point()
ggplot(Output, aes(Resolution, ML.Criterion)) +
  geom_point()
ggplot(Output, aes(Resolution, Spatial.SD)) +
  geom_point()
ggplot(Output, aes(Resolution, Area.km2)) +
  geom_point()
ggplot(Output, aes(Resolution, TotalQuad)) +
  geom_point()

ggplot(Total, aes(Resolution,AvIntensity)) +
  geom_point()
ggplot(Regional, aes(Region, AvIntensity, group = 1))+
  geom_point()+
  geom_line()+
  facet_wrap(~Res)
ggplot(Regional, aes(Res, AvIntensity))+
  geom_point()+
  geom_line()+
  facet_wrap(~Region)

## Make pretty plots 
# Define a clean base theme for all plots
base_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# 1. Model Coefficient Estimate
p1 <- ggplot(Output, aes(Resolution, coef.est)) +
  geom_point(size = 2, color = "#0072B2") +
  geom_smooth(method = "loess", se = FALSE, color = "#0072B2", linetype = "dashed") +
  labs(
    title = "Model Coefficient Estimate",
    x = "Resolution",
    y = "Estimate"
  ) +
  base_theme

# 2. Standard Error
p2 <- ggplot(Output, aes(Resolution, coef.se)) +
  geom_point(size = 2, color = "#D55E00") +
  geom_smooth(method = "loess", se = FALSE, color = "#D55E00", linetype = "dashed") +
  labs(
    title = "Standard Error of Estimate",
    x = "Resolution",
    y = "Standard Error"
  ) +
  base_theme

# 3. Maximum Likelihood Criterion
p3 <- ggplot(Output, aes(Resolution, ML.Criterion)) +
  geom_point(size = 2, color = "#009E73") +
  geom_smooth(method = "loess", se = FALSE, color = "#009E73", linetype = "dashed") +
  labs(
    title = "Maximum Likelihood Criterion",
    x = "Resolution",
    y = "ML Criterion"
  ) +
  base_theme

# 4. Spatial SD
p4 <- ggplot(Output, aes(Resolution, Spatial.SD)) +
  geom_point(size = 2, color = "#CC79A7") +
  geom_smooth(method = "loess", se = FALSE, color = "#CC79A7", linetype = "dashed") +
  labs(
    title = "Spatial Random Effect SD",
    x = "Resolution",
    y = "Spatial SD"
  ) +
  base_theme

# 5. Area of Grid Cells
p5 <- ggplot(Output, aes(Resolution, Area.km2)) +
  geom_point(size = 2, color = "#F0E442") +
  geom_smooth(method = "loess", se = FALSE, color = "#F0E442", linetype = "dashed") +
  labs(
    title = "Area of Grid Cells",
    x = "Resolution",
    y = "Area (km²)"
  ) +
  base_theme

# 6. Number of Quadrature Points
p6 <- ggplot(Output, aes(Resolution, TotalQuad)) +
  geom_point(size = 2, color = "#56B4E9") +
  geom_smooth(method = "loess", se = FALSE, color = "#56B4E9", linetype = "dashed") +
  labs(
    title = "Number of Quadrature Points",
    x = "Resolution",
    y = "Count"
  ) +
  base_theme

combined_plot <- (p1 | p2 | p3) /
  (p4 | p5 | p6)

save_plot(combined_plot, "FigureS8.png", type = "ggplot", width = 900, height = 600)


### Use resolution 2.5 for model selection
  ## read csv file with correct # quadrature points for resolution 2.5 if not already loaded
ASOSCoord_25 <- read.csv("data/processed/ASOSCoord2_5.csv")
view(ASOSCoord_25)

## Create a mesh object 
mesh <- make_mesh(ASOSCoord_25, xy_cols = c("X", "Y"), cutoff = 30)
mesh$mesh$n 
  # 254
plot(mesh)

# small values at presence locations
ASOSCoord_25$wt <- 1e-6
# pseudo-absences: area per quadrature point
n_zeros <- length(which(ASOSCoord_25$Present == 0))
ASOSCoord_25$wt <- ifelse(ASOSCoord_25$Present == 1,
                         1e-6, 269679.1 / n_zeros
)

### Creating Candidate Models 
  ## Using hierarchical approach, first see how depth is favored in model
  ## intercept only, linear, smooth
# Null
Null <- sdmTMB(
  Present / wt ~ 1,
  data = ASOSCoord_25,
  mesh = mesh,
  family = poisson(link = "log"),
  weights = ASOSCoord_25$wt
)


sanity(Null) #Passes sanity check 

# Depth Smooth, no spatial random effect
SmoothDepth <- sdmTMB(
  formula = Present / wt ~ s(Depth),
  data = ASOSCoord_25,
  mesh = mesh,
  family = poisson(link = "log"),
  weights = ASOSCoord_25$wt
)

sanity(SmoothDepth) 

# Depth Linear, no spatial random effect
LinearDepth <- sdmTMB(
  formula = Present / wt ~ Depth,
  data = ASOSCoord_25,
  mesh = mesh,
  family = poisson(link = "log"),
  weights = ASOSCoord_25$wt
)
sanity(LinearDepth) #Passes sanity check

## AIC to determine if depth should be linear or smooth
aic_tableDepth <- AIC(Null, SmoothDepth, LinearDepth)
aic_tableDepth$delta_aic <- round(aic_tableDepth$AIC-min(aic_tableDepth$AIC),2)
aic_tableDepth <- aic_tableDepth %>% arrange(delta_aic)
aic_tableDepth
  # Smooth depth is favored, use smooth


# What timeblock, linear or factor?
TimeNull <- sdmTMB(
  formula = Present / wt ~ s(Depth),
  data = ASOSCoord_25,
  mesh = mesh,
  family = poisson(link = "log"),
  weights = ASOSCoord_25$wt,
  spatial = "on",
  time = "Year"
)
sanity(TimeNull) #Passes sanity check

TimeLinear <- sdmTMB(
  formula = Present / wt ~ s(Depth) + Year,
  data = ASOSCoord_25,
  mesh = mesh,
  family = poisson(link = "log"),
  weights = ASOSCoord_25$wt,
  spatial = "on",
  time = "Year"
)
sanity(TimeLinear) #Passes sanity check

TimeFactor <- sdmTMB(
  formula = Present / wt ~ s(Depth) + as.factor(Year),
  data = ASOSCoord_25,
  mesh = mesh,
  family = poisson(link = "log"),
  weights = ASOSCoord_25$wt,
  spatial = "on",
  time = "Year"
)

sanity(TimeFactor) #Passes sanity check


# AIC to determine if time should be linear or factor
aic_tableTime <- AIC(TimeNull, TimeLinear, TimeFactor)

aic_tableTime$delta_aic <- round(aic_tableTime$AIC-min(aic_tableTime$AIC),2)
aic_tableTime <- aic_tableTime %>% arrange(delta_aic)
aic_tableTime
  # Including time as a linear variable is favored 


### Compare Spatiotemporal Random effects
IID <- sdmTMB(
  formula = Present / wt ~ s(Depth),
  data = ASOSCoord_25,
  mesh = mesh,
  family = poisson(link = "log"),
  weights = ASOSCoord_25$wt,
  spatial = "on",
  time = "Year",
  spatiotemporal = "IID"
)
sanity(IID) # 


RW <- sdmTMB(
  formula = Present / wt ~ s(Depth),
  data = ASOSCoord_25,
  mesh = mesh,
  family = poisson(link = "log"),
  weights = ASOSCoord_25$wt,
  spatial = "on",
  time = "Year",
  spatiotemporal = "RW"
)
sanity(RW) # 


AR1 <- sdmTMB(
  formula = Present / wt ~ s(Depth, k =6),
  data = ASOSCoord_25,
  mesh = mesh,
  family = poisson(link = "log"),
  weights = ASOSCoord_25$wt,
  spatial = "on",
  time = "Year",
  spatiotemporal = "AR1"
)
sanity(AR1) # Does not pass sanity check
summary(AR1)

aic_tableST <- AIC(IID, AR1, RW)
aic_tableST$delta_aic <- round(aic_tableST$AIC-min(aic_tableST$AIC),2)
aic_tableST <- aic_tableST %>% arrange(delta_aic)
aic_tableST
# Random walk is favored

aic_tableTotal <- AIC(IID, AR1, RW, TimeLinear, TimeFactor, TimeNull)
aic_tableTotal$delta_aic <- round(aic_tableTotal$AIC-min(aic_tableTotal$AIC),2)
aic_tableTotal <- aic_tableTotal %>% arrange(delta_aic)
aic_tableTotal
# Model selection favors spatiotemporal random walk, smoothed depth, and linear timeblock fixed effect