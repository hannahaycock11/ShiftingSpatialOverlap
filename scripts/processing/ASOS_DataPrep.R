## Code used to clean FWRI At Sea Observer Survey data for fishing event intensity models ##

# Set working directory 
setwd("~/Research Materials/DepredationPotential_Overlay/DepredationCode")

# Load necessary libraries 
library(reshape2)
library(tidyverse)


# Load in ASOS data set 
Effort <- read.csv("data/raw/ASOS_2023_data_raw.csv")
head(Effort)
nrow(Effort) #438779
summary(Effort)

## Data cleaning protocol 
  # Remove NA values 
      # Use beginLat and beginLon because endLat and endLon are mostly NA
Effort <- Effort %>%  
  filter(!(is.na(Depth))) %>% 
  filter(!(is.na(BeginLat))) %>% 
  filter(!(is.na(BeginLon)))
nrow(Effort) #437449, 1330 NA row's removed 
1330/437449

### Filter depth, latitude, and longitude
EffortCoord <-  Effort %>% 
  filter(between(BeginLat, 24, 31), between(BeginLon, -87.5, -80.5)) %>% 
  mutate(Coast = ifelse((BeginLon < -80.4 & BeginLat < 28) | 
                      (BeginLon <= -82 & BeginLat >= 28),"W","E")) %>%
  filter( Coast == "W") %>% 
  filter(between(Depth, 0, 376)) %>% 
  filter(!(BeginLat > 26.5 & BeginLon > -81.5)) %>% 
  filter(!(BeginLat < 25 & BeginLon < -86))

12920/424529

nrow(EffortCoord)  #424529, 12920 points outside of range removed 
summary(EffortCoord)

### Filter out any outliers/incorrectly recorded data points 
plot(EffortCoord$BeginLon, EffortCoord$BeginLat)

  # Look for outliers 
outlier <- EffortCoord %>% 
  filter(
    (BeginLat < 25 & BeginLon < -84) |
      (between(BeginLat, 29.8, 30) & between(BeginLon, -85.2, -85)) |
      (between(BeginLat, 28.9, 29.1) & between(BeginLon, -86.7, -86.5)) |
      (between(BeginLat, 26, 26.3) & between(BeginLon, -81.7, -81.4)) |
      (between(BeginLat, 26.5, 27) & between(BeginLon, -82.2, -81.9)) |
      (between(BeginLat, 25.5, 26) & between(BeginLon, -81, -80))|
      Depth == 0
  )

nrow(outlier) # 665 rows

outlier %>% count(BeginLat)
  # 13 different sites/latlong points

EffortCoord <- EffortCoord %>%
  mutate(
    outlier =
      (BeginLat < 25 & BeginLon < -84) |
      (between(BeginLat, 29.8, 30) & between(BeginLon, -85.2, -85)) |
      (between(BeginLat, 28.9, 29.1) & between(BeginLon, -86.7, -86.5)) |
      (between(BeginLat, 26, 26.3) & between(BeginLon, -81.7, -81.4)) |
      (between(BeginLat, 26.5, 27) & between(BeginLon, -82.2, -81.9)) |
      (between(BeginLat, 25.5, 26) & between(BeginLon, -81, -80)) |
      Depth == 0
  )

EffortCoord <- EffortCoord %>%
  filter(!outlier)

plot(EffortCoord$BeginLon, EffortCoord$BeginLat)
### Removed 12 implausible points

## Save species list for MRIP code
species.list <- as.data.frame(unique(EffortCoord$Species))
write.csv(species.list, "data/processed/species.list.csv")
  # 288 species caught (there is one blank in the list)


## Summarize data to have each row be a fishing site where a fish was hooked 
  # Currently each row is an individual fish hooked during a trip
  # First select variables of interest 
  # Create a column that records fishing presence (1)
EffortCoord <- EffortCoord %>% 
  select("SeriesID", "YR","MON", "Species",
         "Depth", "BeginLat", "BeginLon", "Project2") %>% 
  group_by(SeriesID, BeginLat, BeginLon) %>% 
  summarize(
            Year = unique(YR, na.rm = TRUE),
            Depth = median(Depth),
            Month = unique(MON, na.rm =TRUE),
            Region = unique(Project2, na.rm = TRUE)
            ) %>% 
  mutate(Present = 1) 
  
nrow(EffortCoord) #13302, number of fishing locations

### Plot number of sites recorded per month 
group <- EffortCoord %>% 
  group_by(Month) %>% 
  summarise(Count = sum(Present))
ggplot(group, aes(Month, Count)) +
  geom_bar(stat = "identity")
  # how many sites were recorded between months 6 and 9
filtered <- EffortCoord %>% 
  filter(between(Month, 6,9))

  nrow(filtered) #5532
  (5532/13302)*100 # 41.6% of trips


## Classifying Regions
  # Can either use the regions recorded by ASOS or manually define regions with boundaries
  # Plotting the differences 
# ASOS regions
ggplot(EffortCoord, aes(BeginLon, BeginLat, color = factor(Region))) + 
  geom_point()
  # There is a lot of overlap between regions and some incorrectly recorded regions

## Manually defining regions 
EffortCoord <- EffortCoord %>% 
  mutate(RegionBound = case_when(
    between(BeginLat, 29, 31) & BeginLon < -84.5 ~ "NW",
    BeginLat > 28.69 & BeginLon >= -84.5 ~ "BB",
    between(BeginLat, 27, 28.69) ~ "TB",
    between(BeginLat, 25.2, 27) ~ "SW",
    BeginLat <= 25.2 ~ "KY",
    TRUE ~ "Other"  # fallback for rows that don't meet the condition
  ))
# Plotting to check classification 
ggplot(EffortCoord, aes(BeginLon, BeginLat, color = factor(RegionBound))) + 
  geom_point() 
  # No overlap and all regions correctly specified

## Comparing the original region count (Project2) to manually classified region count (RegionBound)
RegionCount <- EffortCoord %>% 
  group_by(Region) %>% 
  summarize(OG_n = sum(Present)) 

ClassifiedCount <- EffortCoord %>% 
  group_by(RegionBound) %>%  
  summarize(RegionBoundCount = sum(Present)) %>% 
  rename(Region = RegionBound)
CompareUnboundBound <- left_join(RegionCount, ClassifiedCount, by = "Region") 
# Ignore AL (Alabama, only 2 points)

## Calculate the percent difference between original region definition and defined regions
CompareUnboundBound <- CompareUnboundBound %>% 
  mutate(OGPercent = (OG_n/ 13572)*100,
         BoundPercent = (RegionBoundCount/13572)*100,
         ChangePercent = OGPercent-BoundPercent)
CompareUnboundBound
  # use bound region instead of original region, all less than 0.5% change 

## Rename columns for clarity
EffortCoord$Region <- NULL
EffortCoord <- EffortCoord %>% 
  rename(Region = RegionBound,
         Longitude = BeginLon,
         Latitude = BeginLat) 
nrow(EffortCoord)
  # 13302

## Write csv that will be used in ASOS model
write.csv(EffortCoord, "data/processed/ASOSPresence.csv")

