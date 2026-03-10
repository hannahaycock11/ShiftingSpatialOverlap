## Overlay Depredating Shark Presence and Fishing Event Intensity to calculate depredation risk ##

# Set working directory 
setwd("~/Research Materials/DepredationPotential_Overlay/DepredationCode")

## Function to save all plots to overlay figures folder
save_plot <- function(obj, filename, folder = "figures/Overlay", type = "ggplot", device = "png", width = 800, height = 600, ...) {
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

## load all necessary libraries 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)

# upload files
SharkProb <- read.csv("data/output/BLL_shark_presence_output.csv")
SharkProb$X.1 <- NULL
HookedIntensity <- read.csv("data/output/ASOS_FishingIntensity_Output.csv")
HookedIntensity$X.1 <- NULL

summary(SharkProb)
summary(HookedIntensity)

# check number of rows
nrow(SharkProb) # 176928
nrow(HookedIntensity) # 110580
# Need to filter SharkPresence  year to match fishing event intensity

SharkProb <- SharkProb %>% 
  filter(between(Year, 2009, 2023))
summary(SharkProb$Year)
nrow(SharkProb) # 110580, perfect!

# Filter and rename data sets to merge
SharkProb <- SharkProb %>% 
  select(Longitude, Latitude, X, Y, Year, ProbabilityPresence) 
HookedIntensity <- HookedIntensity %>% 
  select(Longitude, Latitude, X, Y, Year, Intensity)

### merge data sets 
MergedData <- merge(SharkProb, HookedIntensity, by = c("Longitude", "Latitude", "X", "Y", "Year"))
summary(MergedData)

## Calculate Depredation Risk Index
MergedData <- MergedData %>% 
  mutate(DepIndex = ProbabilityPresence * Intensity
         )
summary(MergedData)

time_labels5 <- c(
  "1" = "2009–2011",
  "2" = "2012–2014",
  "3" = "2015–2017",
  "4" = "2018–2020",
  "5" = "2021–2023"
)

###Add Florida 
florida_shape <- st_read("scripts/analysis/Detailed_Florida_State_Boundary.shp")
# Reproject shapefile to UTM Zone 17N to match the bathymetry data
if (st_crs(florida_shape)$epsg != 4326) {
  florida_shape <- st_transform(florida_shape, crs = 4326)
}

risk_spatial_blocks <- MergedData %>%
  mutate(timeblock = cut(Year, 
                         breaks = seq(2009, 2024, by = 3), 
                         right = FALSE, 
                         labels = c(
                           "2009-2011", "2012-2014", "2015-2017", 
                           "2018-2020", "2021-2023"))) %>%
  group_by(timeblock, Latitude, Longitude) %>% 
  summarize(
    mean_risk = mean(DepIndex, na.rm = TRUE),
    .groups = "drop"
  )

risk_spatial_blocks$DepBin <- cut(risk_spatial_blocks$mean_risk,
                         breaks = c(- Inf,0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, Inf), 
                                include.lowest = TRUE)
DepIndexMap <- ggplot() + 
  geom_tile(data = risk_spatial_blocks, 
            aes(x = Longitude, y = Latitude, fill = DepBin)) +  
  geom_sf(data = florida_shape, fill = "lightgray", color = "black") + 
  labs(
       x = "Longitude", 
       y = "Latitude") +
  scale_fill_viridis_d(name = "Depredation Risk \n\ Index (binned)")+
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  facet_wrap(~timeblock) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_sf(expand = FALSE) +  # Keep correct spatial aspect
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
DepIndexMap

save_plot(DepIndexMap, "Figure9.png", type = "ggplot", width = 800, height = 600) 

MergedData <- MergedData %>%
  mutate(Region = case_when(
    between(Latitude, 28.69, 31) & Longitude < -84.5 ~ "NW",
    Latitude > 28.69 & Longitude >= -84.5 ~ "BB",
    between(Latitude, 27, 28.69) ~ "TB",
    between(Latitude, 25.2, 27) ~ "SW",
    Latitude <= 25.2 ~ "KY",
    TRUE ~ "Other"  # fallback for rows that don't meet the condition
  )) 

## Region Averages
RegionAv <- MergedData %>% 
  group_by(Region) %>% 
  summarise(
    MeanProb = mean(ProbabilityPresence),
    MeanInt = mean(Intensity),
    MeanIndex = mean(DepIndex)
  )
RegionAv

mean(MergedData$DepIndex)
mean(MergedData$ProbabilityPresence)
mean(MergedData$Intensity)


## Year Averages
YearAv <- MergedData %>% 
  group_by(Year) %>% 
  summarise(
    MeanProb = mean(ProbabilityPresence),
    MeanInt = mean(Intensity),
    MeanIndex = mean(DepIndex)
  )
YearAv

### total change in probability of presence 
((0.441-0.301)/0.301)*100 
  # 46.5 %, 2009-2023 
### total change in fishing event intensity
((0.0582-0.0336)/0.0336)*100 
  # 73.2 %, 2009-2023 
### total change in depredation index 
((0.0267-0.0114)/0.0114)*100 
  # 134.2105 %, 2009-2023 

#### Calculating Areas of heightened risk, intensity, presence ###  
YearlySum <- MergedData %>% 
  mutate(NewDepIndex = ifelse(DepIndex > 0.01893554, 1, 0),
         NewPresIndex = ifelse(ProbabilityPresence > 0.3688725, 1, 0),
         NewIntensityIndex = ifelse(Intensity > 0.01893554, 1, 0),
         Total = 1) %>% 
  group_by(Year) %>% 
  summarise(
    Total = sum(Total),
    OverDep = sum(NewDepIndex),
    AverageDep = mean(DepIndex),
    OverPres = sum(NewPresIndex),
    AveragePres = mean(NewPresIndex),
    OverInt= sum(NewIntensityIndex),
    AverageInt = mean(NewIntensityIndex)
  ) %>% 
  mutate(
    FractionDep = OverDep/Total,
    PercentDep = FractionDep*100,
    FractionPres = OverPres/Total,
    PercentPres = FractionPres*100,
    FractionInt = OverInt/Total,
    PercentInt = FractionInt*100
  )
YearlySum %>% select(Year, PercentPres, PercentDep, PercentInt)

### total change in probability of presence 
70.5-37.4
  # 33.1 %, 2009-2023 
### total change in fishing event intensity
28.1-21.1
  # 7 %, 2009-2023 
### total change in depredation index 
21.0-11.1
  # 9.9 %, 2009-2023 


# Reshape data for plotting
PlotData <- YearlySum %>%
  select(Year, PercentDep, PercentPres, PercentInt) %>%
  pivot_longer(cols = -Year, names_to = "Metric", values_to = "Percentage")

AreaHeight <- ggplot(PlotData, aes(x = Year, y = Percentage, color = Metric, shape = Metric)) +
  geom_line(linewidth = 1) + 
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 2)) +
  scale_x_continuous(breaks = unique(PlotData$Year)) +
  scale_color_brewer(palette = "Set1", 
                     labels = c("PercentDep" = "Deprivation Risk", 
                                "PercentInt" = "Fishing Event Intensity", 
                                "PercentPres" = "Probability of Presence")) +
  scale_shape_manual(values = c(16, 17, 15), 
                     labels = c("PercentDep" = "Deprivation Risk", 
                                "PercentInt" = "Fishing Event Intensity", 
                                "PercentPres" = "Probability of Presence")) +
  theme_minimal(base_size = 14) +
  labs(x = "Year", y = "Percentage (%)", color = NULL, shape = NULL) +
  theme(
    axis.title = element_text(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, margin = margin(t = 2)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  ) 

save_plot(AreaHeight, "AreaHeight.png", type = "ggplot", width = 800, height = 600)  



### Regional Risk Plot with SE bars
SharkReg<- read.csv("data/output/PresRegion.csv")
SharkReg$X <- NULL

SharkReg <- SharkReg %>% 
  filter(between(Year, 2009, 2023)) %>% 
  rename(
    Lower_Shark = Lower,
    Upper_Shark = Upper,
    SE_Total_Shark = SE_Total
  )
nrow(SharkReg)
summary(SharkReg)

HookedReg <- read.csv("data/output/HookRegion.csv")
HookedReg$X <- NULL

HookedReg <- HookedReg %>% 
  rename(
    Lower_Hook = Lower,
    Upper_Hook = Upper,
    SE_Total_Hook = SE_Total
  )
nrow(HookedReg)
MergeInterval <- merge(HookedReg, SharkReg, by = c("Year", "Region"))

MergeInterval_Summary <- MergeInterval %>%
  mutate(
    DepRisk = Intensity * Presence,
    SE_DepRisk = sqrt((Intensity^2 * SE_Total_Shark^2) + (Presence^2 * SE_Total_Hook^2)),
    Lower_DepRisk = DepRisk - 1.96 * SE_DepRisk,
    Upper_DepRisk = DepRisk + 1.96 * SE_DepRisk
  )

summary(MergeInterval_Summary)

df_total_risk <- MergeInterval_Summary %>%
  group_by(Year) %>%
  summarise(
    DepRisk = mean(DepRisk, na.rm = TRUE),
    Lower_DepRisk = mean(Lower_DepRisk, na.rm = TRUE),
    Upper_DepRisk = mean(Upper_DepRisk, na.rm = TRUE),
    SE_DepRisk = mean(SE_DepRisk, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Region = "Total")


df_with_total_risk <- bind_rows(MergeInterval_Summary, df_total_risk) %>%
  mutate(Year = as.numeric(as.character(Year)))

# --- 4. Region labels and color palette ---
region_labels <- c(
  "NW" = "Northwest",
  "BB" = "Big Bend",
  "TB" = "Tampa",
  "SW" = "Southwest",
  "KY" = "Keys",
  "Total" = "West Florida"
)

five_colors_okabe_ito <- c(
  "Northwest" = "#0072B2",     # Blue
  "Big Bend"  = "#009E73",     # Green
  "Tampa"     = "#E69F00",     # Gold
  "Southwest" = "#56B4E9",     # Sky Blue
  "Keys"      = "#CC79A7",     # Pink
  "West Florida" = "black"
)

df_with_total_risk$RegionFull <- region_labels[df_with_total_risk$Region]

# Create lighter versions of the colors for the ribbons
five_colors_ribbon <- sapply(five_colors_okabe_ito, function(x) adjustcolor(x, alpha.f = 0.3))

RegionRisk <- ggplot(df_with_total_risk, aes(x = Year, y = DepRisk, color = RegionFull, group = RegionFull)) +
  geom_ribbon(aes(ymin = Lower_DepRisk, ymax = Upper_DepRisk, fill = RegionFull), color = NA) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = five_colors_okabe_ito, name = "Region") +
  scale_fill_manual(values = five_colors_ribbon, guide = "none") +  
  scale_x_continuous(
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Year",
    y = "Average Depredation Risk"
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  facet_wrap(~RegionFull)


RegionRisk

save_plot(RegionRisk, "Figure10.png", type = "ggplot", width = 800, height = 600) 
