## Code used to clean NOAA shark BLL data for Depredating shark models ##

# Set working directory 
setwd("~/Research Materials/DepredationPotential_Overlay/DepredationCode")

## load necessary libraries
library(tidyverse)
library(readxl)
library(lubridate)

## Read in BLL data and rename columns 
# rename species to common name and select only species of interest
BLL <- read_xlsx("data/raw/2023_SEFSC_BLL_raw.xlsx") %>% rename ( 
  spinner = "CARCHARHINUS_BREVIPINNA", bull = "CARCHARHINUS_LEUCAS", dusky = "CARCHARHINUS_OBSCURUS",
  sandbar = "CARCHARHINUS_PLUMBEUS", tiger ="GALEOCERDO_CUVIER", nurse = "GINGLYMOSTOMA_CIRRATUM",
  lemon = "NEGAPRION_BREVIROSTRIS", scallopedham = "SPHYRNA_LEWINI", greatham = "SPHYRNA_MOKARRAN", 
  reef ="CARCHARHINUS_PEREZII", sandtiger = "CARCHARIAS_TAURUS", blacktip = "CARCHARHINUS_LIMBATUS", 
  Longitude = "FLYLAST_OUTLON", Latitude = "FLYLAST_OUTLAT", Depth = `FLYLAST_OUTDEPTH m`, Temp = `TEMP C`) %>% 
  select(!c(reef, sandtiger))


## Read in Lengths data and rename 
Lengths <- read_csv("data/raw/Lengths_BLL _SEFSC_raw.csv") %>% 
  rename(FL = `FORK mm`, TL = `NATURAL_TOTAL mm`, TL2 = `TOTAL mm`) %>%
  mutate(
    TAXON = recode(TAXON, `CARCHARHINUS PLUMBEUS` = 'sandbar', 
                   `CARCHARHINUS LEUCAS` = 'bull', `SPHYRNA LEWINI` = 'scallopedham',
                   `CARCHARHINUS LIMBATUS` = 'blacktip', `GINGLYMOSTOMA CIRRATUM` = 'nurse', 
                   `CARCHARHINUS OBSCURUS` = 'dusky', `GALEOCERDO CUVIER`= 'tiger',  
                   `SPHYRNA MOKARRAN` = 'greatham', `NEGAPRION BREVIROSTRIS` = 'lemon',
                   `CARCHARHINUS BREVIPINNA` = 'spinner', `CARCHARHINUS PEREZII` = 'reef',
                   `CARCHARIAS TAURUS` = 'sandtiger'
    )) %>% 
  filter(!TAXON == 'reef', !TAXON == 'sandtiger')
  # reef and sandtiger omitted due to low number

## Check for any data entry errors
Lengths %>% filter(FL > TL) %>% nrow()
  # 0

## Goal is to use Fork length (FL) rather than total length (TL) - but some individuals are missing TL measurements
  # FL is used because it helps to standardize the jaw size between species 
  # Run a linear regression model for each species to predict FL from TL 
  # Find number individuals per species to see if it is large enough to run linear regression 
LengthsCount <- Lengths
LengthsCount$TAXON <- factor(LengthsCount$TAXON)
LengthsCount <- LengthsCount %>%
  mutate(Count = 1) %>%
  group_by(TAXON) %>%
  mutate(SpeciesCount = sum(Count))
LengthsCount <- LengthsCount %>%
  distinct(TAXON,SpeciesCount)
LengthsCount
  #sample size high enough for all

## How many individuals are missing FL measurements
Lengths %>% filter(!TAXON == 'nurse') %>% filter(is.na(FL)) %>% nrow() 
nrow(Lengths)
(1438/8590)*100
  # 16.7%

SummaryMinMax <- Lengths %>% 
  filter(!is.na(FL)) %>% 
  filter(!is.na(TL)) %>% 
  group_by(TAXON) %>% 
  summarize(
    max = max(FL, na.rm = TRUE),
    min = min(FL, na.rm =TRUE),
    mean = mean(FL, na.rm =TRUE)
  )
SummaryMinMax


## calculate FL for all species using a linear regression model 
SandbarBLL <- Lengths %>% filter(TAXON == "sandbar")
SandbarModel <- lm(FL~TL, data=SandbarBLL)
summary(SandbarModel)
  # FL = TL*0.830 +10.949

BullBLL <- Lengths %>% filter(TAXON == "bull")
BullModel <- lm(FL~TL, data=BullBLL)
summary(BullModel)
  # FL = TL*0.843 +11.732

ScallopedHamBLL <- Lengths %>% filter(TAXON == "scallopedham")
ScallopedHamModel <- lm(FL~TL, data=ScallopedHamBLL)
summary(ScallopedHamModel)
  # FL = TL*0.806 -33.897

BlacktipBLL <- Lengths %>% filter(TAXON == "blacktip")
BlacktipModel <- lm(FL~TL, data=BlacktipBLL)
summary(BlacktipModel)
  # FL = TL*0.837 -2.675

DuskyBLL <- Lengths %>% filter(TAXON == "dusky")
DuskyModel <- lm(FL~TL, data=DuskyBLL)
summary(DuskyModel)
  # FL =TL*0.803 +58.541

TigerBLL <- Lengths %>% filter(TAXON == "tiger")
TigerModel <- lm(FL~TL, data=TigerBLL)
summary(TigerModel)
  # FL =TL*0.857 - 80.590

GreatHamBLL <- Lengths %>% filter(TAXON == "greatham")
GreatHamModel <- lm(FL~TL, data=GreatHamBLL)
summary(GreatHamModel)
  # FL =TL*0.798 - 32.895

LemonBLL <- Lengths %>% filter(TAXON == "lemon")
LemonModel <- lm(FL~TL, data=LemonBLL)
summary(LemonModel)
  # FL = TL*0.791 + 128.327

SpinnerBLL <- Lengths %>% filter(TAXON == "spinner")
SpinnerModel <- lm(FL~TL, data=SpinnerBLL)
summary(SpinnerModel)
  # FL =TL*0.845 - 12.502


## Use estimate of intercept and TL estimated by linear regression models to calculate FL for each species 
  # Calculates a FL for all individuals, even if a FL measurement was taken. 
  # TL was used for Nurse Sharks due to inconsistencies in measuring FL 

Lengths <- Lengths %>% mutate(NewFL = case_when(
  TAXON == "tiger" ~ -80.590 + 0.857*TL,
  TAXON == "bull" ~ 11.732 + 0.843*TL,
  TAXON == "sandbar" ~ 10.949 + 0.830*TL,
  TAXON == "dusky" ~ 58.541 + 0.803*TL,
  TAXON == "scallopedham" ~ -33.897 + 0.806*TL,
  TAXON == "greatham" ~ -32.895 + 0.798*TL,
  TAXON == "spinner" ~  -12.502 + 0.845*TL,
  TAXON == "lemon" ~ 128.327 + 0.791*TL,
  TAXON == "nurse" ~ TL,
  TAXON == "blacktip" ~ -2.678 + 0.835*TL
))


##  UseFL uses measured FL from survey and uses calculated FL when FL was originally NA
Lengths <- Lengths %>% 
  mutate(UseFL = ifelse(is.na(FL),NewFL,FL),
  )

## Check for blanks 
Lengths %>% filter(is.na(UseFL)) %>% nrow()
  # 346 rows missing TL

## Some sharks have no TL (natural total) but have TL2 (stretched total)
  # Determine difference (in millimeters) between TL (natural) and TL2 (stretched total)
DifTLand2 <- Lengths %>% 
  filter(!is.na(TL), !is.na(TL2)) %>% 
  group_by(TAXON) %>% 
  summarise(AverageTL = mean(TL),
            AverageTL2 = mean(TL2),
            Difference = AverageTL2 - AverageTL) %>% 
  select(TAXON, AverageTL, AverageTL2, Difference)

DifTLand2
  # Use TL2(stretched total), when natural total is unavailable because difference between the two is 5-10 cm
  # which is 3.33% - 6.67% of variance in sharks that are 150 cm 

## Code that calculates FL using TL2 if TL is not available 
Lengths <- Lengths %>% mutate(NewFL2 = case_when(
  TAXON == "tiger" ~ -80.59+0.86*TL2,
  TAXON == "bull" ~ 11.73+0.84*TL2,
  TAXON == "sandbar" ~ 10.95+0.83*TL2,
  TAXON == "dusky" ~ 58.54+0.80*TL2,
  TAXON == "scallopedham" ~ -33.90+0.81*TL2,
  TAXON == "greatham" ~ -32.89+0.80*TL2,
  TAXON == "spinner" ~  -12.50+0.84*TL2,
  TAXON == "lemon" ~ 128.33+0.79*TL2,
  TAXON == "nurse" ~ TL2,
  TAXON == "blacktip" ~ -2.86+0.84*TL2
))


## Code that will use original FL, and if original is not available it will use calculated FL from NewFL,
   # if that is not available, use TL2
Lengths <- Lengths %>% 
  mutate(UseFLFinal = ifelse(is.na(UseFL),NewFL2,UseFL),
  )
## Check if there are any NA's
Lengths %>% filter(is.na(UseFLFinal))%>% nrow() 
  # 1 
## Filter out of data 
Lengths <- Lengths %>% filter(!is.na(UseFLFinal))


## combine lengths and create present and  total columns based on size
  # Nurse sharks must exceed 200 cm (not 150 cm) since TL was used
PresentandTotal <- Lengths %>%
  mutate(
    Over150Present = case_when(
      TAXON == "nurse" ~ ifelse(UseFLFinal > 2000, 1, 0),
      TRUE ~ ifelse(UseFLFinal > 1500, 1, 0)
    )) %>%
  group_by(STATIONKEY) %>%
  summarise(
    TotalOver150 = sum(Over150Present)
  ) %>%
  mutate(
    Over150Present = ifelse(TotalOver150 > 0, 1, 0)
  )

## Total presence counts and total individuals > 150
sum(PresentandTotal$Over150Present) # 1425
sum(PresentandTotal$TotalOver150) # 2543

## Combine data sets to station information 
BLLJoin <- left_join(BLL, PresentandTotal, by="STATIONKEY")
BLLJoin <- BLLJoin %>% 
  mutate(Over150Present = replace_na(Over150Present,0), TotalOver150 = replace_na(TotalOver150,0))
summary(BLLJoin)

## Total stations, presences, and total used in analysis per year block
BLLGroup <- BLLJoin %>%   
  mutate(
    YearBlock = case_when(
    between(YEAR, 2000, 2002) ~ "2000-2002", 
    between(YEAR, 2003, 2005) ~ "2003-2005",  
    between(YEAR, 2006, 2008) ~ "2006-2008",
    between(YEAR, 2009, 2011) ~ "2009-2011",  
    between(YEAR, 2012, 2014) ~ "2012-2014",
    between(YEAR, 2015, 2017) ~ "2015-2017",
    between(YEAR, 2018, 2020) ~ "2018-2020",
    between(YEAR, 2021, 2023) ~ "2021-2023"),
    Stations = 1) %>% 
  group_by(YEAR) %>% 
  summarise(
    TotalStations = sum(Stations),
    TotalSharks = sum(TotalOver150),
    TotalPresence = sum(Over150Present)
  )
  # NA Yearblock is from years 1995-1999 that were not used in analaysis
view(BLLGroup)
#write BLL as CSV for presence absence model 
write.csv(BLLJoin, "data/processed/BLLPresence.csv")



#### Summarize species composition of catch ####
### Summarize counts by species and year
species_counts <- Lengths %>%
  group_by(YEAR, TAXON) %>%
  summarise(n = n(), .groups = "drop")
species_counts

### now lets only look at the larger individuals 
species_counts_Big <- Lengths %>%
  group_by(YEAR, TAXON) %>%
  filter(UseFLFinal > 1500) %>% 
  summarise(n = n_distinct(STATIONKEY), total_catch = n(), .groups = "drop")

# Step 2: compute effort (number of stations) per year
effort_per_year_Big <- Lengths %>%
  group_by(YEAR) %>%
  summarise(effort = n_distinct(STATIONKEY), .groups = "drop")

# Step 3: join and standardize
species_cpue_Big <- species_counts_Big %>%
  left_join(effort_per_year, by = "YEAR") %>%
  mutate(cpue = total_catch / effort)


blue_palette <- c(
  "#08306B",  # navy
  "#08519C",
  "#2171B5",
  "#4292C6",
  "#6BAED6",
  "#9ECAE1",
  "#C6DBEF",
  "#DEEBF7",
  "#3182BD",
  "#084594"
)

taxon_labels <- c(
  bull         = "Bull Shark",
  greatham     = "Great Hammerhead",
  nurse        = "Nurse Shark",
  sandbar      = "Sandbar Shark",
  scallopedham = "Scalloped Hammerhead",
  tiger        = "Tiger Shark",
  blacktip     = "Blacktip Shark",
  spinner      = "Spinner Shark",
  dusky        = "Dusky Shark",
  lemon        = "Lemon Shark"
)


ggplot(species_cpue_Big,
       aes(x = factor(YEAR), y = cpue, fill = TAXON)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(
    values = blue_palette,
    labels = taxon_labels
  ) +
  labs(
    x = "Year",
    y = "Number of individuals >150 cm FL per Station",
    fill = "Species"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    legend.text = element_text(size = 11)
  )

ggplot(species_cpue_Big,
       aes(x = factor(YEAR), y = cpue, fill = TAXON)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = blue_palette, guide = "none") +
  labs(
    x = NULL,
    y = "Number of individuals >150 cm FL"
  ) +
  facet_wrap(
    ~ TAXON,
    scales = "free_y",
    ncol = 2,
    labeller = labeller(TAXON = taxon_labels)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank()
  )

