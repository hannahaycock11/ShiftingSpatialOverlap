 ## Fishing event intensity models ##

# Set working directory 
setwd("~/Research Materials/DepredationPotential_Overlay/DepredationCode")

## Load all necessary libraries 
library(tidyverse)
library(ggplot2)
library(sdmTMB)
library(marmap)
library(pROC)

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

## Read in ASOS data 
ASOSCoord <- read.csv("data/processed/ASOSPresence.csv")

## Making depth negative to match bathymetry data from marmap 
ASOSCoord$Depth <- ASOSCoord$Depth*-1
nrow(ASOSCoord) 
  # 13302

ASOSGroup <- ASOSCoord %>%   
  group_by(Year) %>% 
  summarise(
    TotalPresence = sum(Present)
  )

### Accounting for sampling regime
## weighting the number of sites per region per year 
# first counting number of surveyed fishing sites per region in each year
Rescale <- ASOSCoord %>% 
  group_by(Year, Region) %>% 
  summarise(NumFS = sum(Present), .groups = 'drop') %>%
  complete(Year, Region, fill = list(NumFS = 0))

sum(ASOSCoord$Present) # 13302
sum(Rescale$NumFS) # 13302
  # summarized correctly 

##  Next, calculate proportion of trips taken in each region per year ASOS and MRIP 
ASOSCoverage <- ASOSCoord %>%
  distinct(SeriesID, .keep_all = TRUE) %>% 
  group_by(Year, Region) %>% 
  summarise(RegionCountASOS = n(), .groups = "drop") %>% 
  complete(Year, Region, fill = list(RegionCountASOS = 0))

ASOSCoverageYear <- ASOSCoverage %>% 
  group_by(Year) %>% 
  summarise(ASOSYear = sum(RegionCountASOS))

ASOSProp <- left_join(ASOSCoverage, ASOSCoverageYear, by = "Year")

## Merge regional number of fishing sites and ASOS proportion data sets 
Join <- left_join(Rescale, ASOSProp, by = c("Year", "Region"))
  # replace zeros with a small number (weight cant use zero in calculation)
Join$NumFS <- replace(Join$NumFS, Join$NumFS == 0, 0.1)
Join$RegionCountASOS <- replace(Join$RegionCountASOS, Join$RegionCountASOS == 0, 0.1)
summary(Join)

## Read in MRIP region trip estimates 
  # sum of trips in each county per region / annual West florida shelf trips
MRIPCountyNum <- read.csv("data/processed/MRIP_CountyNum_2023.csv")
head(MRIPCountyNum)
MRIPCountyNum$X <- NULL

## Count number of MRIP trips per year
MRIPCountyNumYear <- MRIPCountyNum %>% 
  group_by(Year) %>% 
  summarize(YearMRIP = sum(CountyMRIP))

## Merge MRIP data together
Merge_MRIP <- left_join(MRIPCountyNum, MRIPCountyNumYear, by = "Year")

## Merge MRIP, fishing sites, and ASOS region count data 
C.Factor <- left_join(Join, Merge_MRIP, by = c("Year", "Region"))
C.Factor <- C.Factor %>% 
  select(Year, Region, NumFS, RegionCountASOS, CountyMRIP) 
summary(C.Factor)
head(C.Factor)

## Create weight to have ASOS coverage equal proportion of reported MRIP trips
  # Current: Number of fishing sites/Total ASOS observed trips per year 
  # Target: Number of MRIP Trips per REgion/Total MRIP Trips per year
  # Weight: Target/Current
sum(C.Factor$RegionCountASOS) # 5154.2
sum(C.Factor$CountyMRIP) # 11766053
C.Factor <- C.Factor %>% 
  mutate(Current = RegionCountASOS/sum(RegionCountASOS), 
         Target = CountyMRIP/sum(CountyMRIP),
         Weight = Target/Current) 

## Multiply the weight by the total number of fishing sites 
  # Calculates the weighted number of fishing sites 
  # This is the # of fishing sites needed in each region to match the proportion of MRIP trips per region 
C.Factor <- C.Factor %>% 
  group_by(Region) %>% 
  mutate(RegionSite = sum(NumFS),
         TotalTrip = sum(RegionCountASOS),
         SiteTrip = RegionSite/TotalTrip) %>% 
  ungroup() %>% 
  mutate(PointWeight = RegionCountASOS*Weight,
         WeightedSites = PointWeight*SiteTrip) 
summary(C.Factor$SiteTrip)


Table <- C.Factor %>% group_by(Year) %>% 
  summarise(SumTrip = sum(PointWeight),
            Sum = sum(WeightedSites),
            TotalMRIP = sum(CountyMRIP))

## Plot of number of calculated weighted sites per year and MRIP trips per year
ggplot(Table, aes(Year, Sum)) +
  geom_line()
ggplot(Table, aes(Year, TotalMRIP)) +
  geom_line()
  # both follow the same pattern

## Creating a dataframe that has the year, region, site, and weight information
weight.df <- C.Factor %>% 
  select(Year, Region, NumFS, WeightedSites)

## Determining how each region-year will be resampled to achieve proper scaling: 
  # For this model to be scaled properly, the # of fishing sites needs to = the weighted # of fishing sites
  # There are 3 scenarios: 
    # 1. # Fishing sites = Weighted Fishing Sites
      # Will use the same Lon/lat coords (no resampling)
    # 2. # Fishing Sites > Weighted Fishing Sites
      # Resample Lon/Lat coords from that region-year's fishing site locations (Num.FS) 
    # 3. # Fishing Sites < Weighted Fishing Sites
      # Resample Lon/Lat coords from the region's total fishing site locations (Loc)


## An aside, technically, I could resample each year-region using the regional reservoir of fishing locations (loc)
  # BUT, I am trying to keep the fishing locations as accurate to the year-region as possible (just in case it's important)
  # So I'm trying to resample from points only from that year and region when it's possible 
  # some region-years do not have many points (or any at all), which is when resampling the location reservoir is handy

wdf <- weight.df %>% 
  mutate(NumFS = round(NumFS),
         WeightedSites = round(WeightedSites),
         ID = case_when(
           NumFS == WeightedSites ~ "AC",  # 1
           NumFS > WeightedSites  ~ "FS",  # 2
           NumFS < WeightedSites  ~ "LOC"  # 3
         ))

IDCount <- wdf %>% 
  count(ID) 
IDCount
  # No scenario 1's

## Data set that contains fishing site locations by year and region
FS <- ASOSCoord %>% 
  select(Year,Region, Latitude, Longitude, Present, Depth) 
## Data set that contains fishing locations ONLY by region 
LOC <- ASOSCoord  %>% 
  select(Year,Region, Latitude, Longitude, Present, Depth) 

## Filter region and years by how they need to be resampled
FS.wdf <- wdf %>% 
  filter(ID == 'FS')

LOC.wdf <- wdf %>% 
  filter(ID == 'LOC')

## Check to see that filter worked
nrow(FS.wdf) 
  # 35
nrow(LOC.wdf) 
  # 40
nrow(wdf) 
  # 75 
  # filtered correctly 


## Loop to Resample Region-Year Fishing Sites (scenario 2's)
  # For each Year and Region, use the # of Point Weights to randomly select that number of FS locations
FS.list <- data.frame(Longitude = numeric(), Latitude = numeric(), Year = numeric(), 
                      Region = character(), Depth = numeric(), TotalHooked = numeric())

for (i in 1:nrow(FS.wdf)){
  year <- FS.wdf$Year[i]
  region <- FS.wdf$Region[i]
  nsamp <- FS.wdf$WeightedSites[i]
  
  sub.fs = FS %>% filter(Year == year, Region == region)
  
  rand = sample(1:nrow(sub.fs),size=nsamp, replace = FALSE)
  
  FS.list = bind_rows(FS.list, sub.fs[rand,])
  
}

nrow(FS.list) 
  # 6717
sum(FS.wdf$WeightedSites)
  # 6717

## Check to make sure loop worked
checkFS <- FS.list %>% 
  group_by(Year, Region) %>% 
  summarise(count = sum(Present))
c1 <- FS.wdf %>% 
  select(Year,Region,WeightedSites)
cfs <- left_join(checkFS, c1, by = c("Year", "Region"))
cfs <- cfs %>% 
  mutate(check = count - WeightedSites) 
cfs %>% filter(check == 0) %>%  nrow() 
  # 35
nrow(cfs) 
  # 35
  # Loop worked!

## Loop to Resample Region Reservoir Fishing Sites (scenario 3's)
# For each Year and Region, use the # of Point Weights to randomly select that number of LOC locations
LOC.list <- data.frame(Longitude = numeric(), Latitude = numeric(), Year = numeric(), 
                       Region = character(), Depth = numeric(), TotalEffort = numeric())

for (i in 1:nrow(LOC.wdf)){
  year = LOC.wdf$Year[i]
  region = LOC.wdf$Region[i]
  nsamp = LOC.wdf$WeightedSites[i]
  sub.loc = LOC
  
  sub.loc = LOC %>% filter(Region == region)
  
  rand = sample(1:nrow(sub.loc),size=nsamp, replace = FALSE)
  sub.loc$Year = year
  
  LOC.list = bind_rows(LOC.list, sub.loc[rand,])
  
}

nrow(LOC.list)
  # 6303
sum(LOC.wdf$WeightedSites)
  # 6303

## check to make sure loop worked
checkLOC <- LOC.list %>% 
  group_by(Year, Region) %>% 
  summarise(count = sum(Present))
c2 <- LOC.wdf %>% 
  select(Year,Region,WeightedSites)
cloc <- left_join(checkLOC, c2, by = c("Year", "Region"))
cloc <- cloc %>% 
  mutate(check = count - WeightedSites) 
cloc %>% filter(check == 0) %>%  nrow() # 40
nrow(cloc) # 40

combined <- bind_rows(FS.list, LOC.list)
nrow(combined) 
  # 13020
sum(wdf$WeightedSites) 
  # 13020
  # loop worked!

## check to see how representative my sample of the data is 
ggplot() +
  geom_point(data = ASOSCoord, aes(Longitude, Latitude),
  size = 1, alpha = 0.5,  color = "blue") +
  geom_point(data = combined, aes(Longitude, Latitude),
             size = 1, alpha = 0.5,  color = "red") +
  facet_wrap(~Year)

## Select columns of interest
PresObsForHire <- combined %>% 
  select(Year, Region, Longitude, Latitude, Present, Depth)

## write csv file to use for sensitivity analysis
write.csv(PresObsForHire, "data/processed/ASOS_SensitivySelectivity.csv")
  # Go to IntensitySensitivity.R to see sensitivity analysis and model selection 

### Fishing event intensity model 
## Read in data with quadrature points (number selected during sensitivity analysis)
ForHireData <- read.csv("data/processed/ASOSCoord2_5.csv")
ForHireData$X.1 <- NULL

## create a mesh object 
mesh <- make_mesh(ForHireData, xy_cols = c("X", "Y"), cutoff = 30)
mesh$mesh$n # 249
plot(mesh)

## small values at presence locations
ForHireData$wt <- 1e-6
## Add quadrature weights
n_zeros <- length(which(ForHireData$Present == 0))
ForHireData$wt <- ifelse(ForHireData$Present == 1,
                         1e-6, 269679.1 / n_zeros
)
  # 269679.1 is area of W. Florida shelf (see sensitivity analysis for calculation)

## Down Weighted Poisson Regression Model
HookedFish <- sdmTMB(
  formula = Present / wt ~ s(Depth),
  data = ForHireData,
  mesh = mesh,
  family = poisson(link = "log"),
  weights = ForHireData$wt,
  spatial = "on",
  time = "Year",
  spatiotemporal = "RW"
)

sanity(HookedFish) #Passes sanity check
summary(HookedFish)

tidy(HookedFish, conf.int = TRUE)
tidy(HookedFish, "ran_pars", confint = TRUE)

### Model Diagnostics 
## simulate data 
fishsim <- simulate(HookedFish, nsim = 500, type = "mle-mvn")
# Compare fraction of zeroes 
sum(ForHireData$Present == 0) / length(ForHireData$Present)
# 0.9244495
sum(fishsim == 0)/length(fishsim)
# 0.9344329


### Predictions
## Using the same prediction grid ### 
  # Define the area: lat and lon range
lon_range <- c(-87.5, -80.5)  # longitude range
lat_range <- c(24, 31)   
  
# Get bathymetry data
bathy_data <- getNOAA.bathy(lon1 = lon_range[1], lon2 = lon_range[2], 
                            lat1 = lat_range[1], lat2 = lat_range[2], 
                            resolution = 3) #~5.55 km
 bathy_new <- as.xyz(bathy_data)
names(bathy_new) <- c("Longitude","Latitude", "Depth")
## Need to assign UTM values
bathy_new <-sdmTMB::add_utm_columns(bathy_new,ll_names=c("Longitude","Latitude"),units="km", 
                                    utm_crs = 32617)
predgrid <- bathy_new
## Filter Lat, Lon, and Depth to match data 
predgrid <- predgrid %>% 
  mutate(Coast = ifelse((Longitude < -80.4 & Latitude < 28) | 
                          (Longitude <= -82 & Latitude >= 28),"W","E")) %>%
  filter( Coast == "W") %>% 
  filter( between(Depth, -376, 0)) %>% 
  filter(!(Latitude > 26.5 & Longitude > -81.5)) %>% 
  filter(!(Latitude < 25 & Longitude < -86))
predgrid$Coast<-NULL
plot(predgrid$Longitude, predgrid$Latitude) 
predgrid <- replicate_df(predgrid, "Year", unique(ForHireData$Year))

## Predict
HookedPred <- predict(HookedFish, newdata = predgrid)
HookedPred <- HookedPred %>% 
  mutate(Intensity = exp(est))

## write csv of fishinb event intensity model
write.csv(HookedPred, "data/output/ASOS_FishingIntensity_Output.csv")


## Look at Fixed, spatial, and spatiotemporal random effects
Fixed <- ggplot() +
  geom_tile(data = HookedPred, 
            aes(X, Y, fill = exp(est_non_rf)), 
            width = 6, height = 6) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Fishsing Event Intensity Model - Prediction with Fixed Effects Only")+
  coord_fixed()
save_plot(Fixed, "FigureS12a.png", type = "ggplot", width = 800, height = 600)

Random <- ggplot() +
  geom_tile(data = HookedPred, 
            aes(X, Y, fill = omega_s), 
            width = 6, height = 6) +
  scale_fill_gradient2() +
  ggtitle("Fishsing Event Intensity Model  - Spatial Random Effects Only")+
  coord_fixed()
save_plot(Random, "FigureS12b.png", type = "ggplot", width = 800, height = 600)

Spatiotemporal <- ggplot() +
  geom_tile(data = HookedPred, 
            aes(X, Y, fill = epsilon_st), 
            width = 6, height = 6) +
  scale_fill_gradient2() +
  facet_wrap(~Year, ncol = 5) + 
  ggtitle("Fishsing Event Intensity Model  - Spatiotemporal Random Effects Only")+
  coord_fixed()
save_plot(Spatiotemporal, "FigureS12c.png", type = "ggplot", width = 800, height = 600)


### Add Florida 
library(sf)
florida_shape <- st_read("scripts/analysis/Detailed_Florida_State_Boundary.shp")
# Reproject shapefile to UTM Zone 17N to match the bathymetry data
if (st_crs(florida_shape)$epsg != 4326) {
  florida_shape <- st_transform(florida_shape, crs = 4326)
}

time_labels2 <- c(
  "1" = "2009–2011",
  "2" = "2012–2014",
  "3" = "2015–2017",
  "4" = "2018–2020",
  "5" = "2021–2023"
)

intensity_spatial_blocks <- HookedPred %>%
  mutate(timeblock = cut(Year, 
                         breaks = seq(2009, 2024, by = 3), 
                         right = FALSE, 
                         labels = c(
                                    "2009-2011", "2012-2014", "2015-2017", 
                                    "2018-2020", "2021-2023"))) %>%
  group_by(timeblock, Latitude, Longitude) %>% 
  summarize(
    mean_intensity = mean(Intensity, na.rm = TRUE),
    .groups = "drop"
  )


intensity_spatial_blocks$Intensity_bin <- cut(intensity_spatial_blocks$mean_intensity, 
                                breaks = c(- Inf, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, Inf), 
                                include.lowest = TRUE)

IntensityMap <- ggplot() + 
  geom_tile(data = intensity_spatial_blocks, 
            aes(x = Longitude, y = Latitude, fill = Intensity_bin)) +  
  geom_sf(data = florida_shape, fill = "lightgray", color = "black") + 
  labs(
       x = "Longitude", 
       y = "Latitude") +
scale_fill_viridis_d(name = "Intensity (binned)")+
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))+
  facet_wrap(~timeblock) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_sf(expand = FALSE) +  # Keep correct spatial aspect
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

IntensityMap
save_plot(IntensityMap, "Figure7.png", type = "ggplot", width = 800, height = 600)

### Calculating Regional predictions uncertainty
 
# Simulate from joint posterior
sim_Hook <- predict(HookedFish, newdata = predgrid, nsim = 1000)

### add regions and year
sim_df_Hook_Region <- as.data.frame(sim_Hook)
sim_df_Hook_Region$Year <- predgrid$Year

predgrid <- predgrid  %>%
  mutate(Region = case_when(
    between(Latitude, 28.69, 31) & Longitude < -84.5 ~ "NW",
    Latitude > 28.69 & Longitude >= -84.5 ~ "BB",
    between(Latitude, 27, 28.69) ~ "TB",
    between(Latitude, 25.2, 27) ~ "SW",
    Latitude <= 25.2 ~ "KY",
    TRUE ~ "Other"  # fallback for rows that don't meet the condition
  )) 

sim_df_Hook_Region$Region <- predgrid$Region 

sim_long_Hook_Region <- pivot_longer(sim_df_Hook_Region, 
                                     cols = starts_with("V"), 
                                     names_to = "sim", 
                                     values_to = "value")

agg_sim_Hook_Region <- sim_long_Hook_Region %>%
  group_by(Year, Region, sim) %>%
  summarise(sum_pred = mean(exp(value)), .groups = "drop")

region_summary_df_Hook <- agg_sim_Hook_Region %>%
  group_by(Year, Region) %>%
  summarise(
    Intensity = mean(sum_pred),
    Lower = quantile(sum_pred, 0.025),
    Upper = quantile(sum_pred, 0.975),
    SE_Total = sd(sum_pred)
  )
 
write.csv(region_summary_df_Hook,"data/output/HookRegion.csv")

# Making nice plot 
df_total3 <- region_summary_df_Hook %>%
  group_by(Year) %>%
  summarise(
    Intensity = mean(Intensity, na.rm = TRUE),
    Lower          = mean(Lower, na.rm = TRUE),
    Upper          = mean(Upper, na.rm = TRUE),
    SE_Total       = mean(SE_Total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Region = "Total")

# bind back into main dataset
df_with_total3 <- bind_rows(region_summary_df_Hook, df_total3)
summary(df_with_total3)

df_with_total3 <- df_with_total3 %>%
  mutate(Year = as.numeric(as.character(Year)))

region_labels <- c(
  "NW" = "Northwest",
  "BB" = "Big Bend",
  "TB" = "Tampa",
  "SW" = "Southwest",
  "KY" = "Keys",
  "Total" = "West Florida"
)

df_with_total3$RegionFull <- region_labels[df_with_total3$Region]

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

RegionIntensity <- ggplot(df_with_total3, aes(x = Year, y = Intensity, color = RegionFull, group = RegionFull)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = RegionFull), color = NA) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = five_colors_okabe_ito, name = "Region") +
  scale_fill_manual(values = five_colors_ribbon, guide = "none") +  
  scale_x_continuous() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Year",
    y = "Fishing Event Intensity"
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
  )+
  facet_wrap(~RegionFull)
RegionIntensity
save_plot(RegionIntensity, "Figure8.png", type = "ggplot", width = 800, height = 600) 


