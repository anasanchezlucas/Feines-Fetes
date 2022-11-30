# library
library(ggplot2)
library(ggridges)
library(ggplot2)
library(purrr)
library(dplyr)

#Data Load
HeightWeight <- read.csv("~/Downloads/SOCR-HeightWeight.csv")
temp <- read.csv("~/Downloads/archive/GlobalLandTemperaturesByCountry.csv")
spain <-temp[temp$Country == 'Spain',]
species <- read.csv("~/Downloads/species.csv")

# Scatter
ggplot(HeightWeight, aes(x=Height.Inches., y=Weight.Pounds.)) +
  geom_point() + xlab("Height (inches)") + ylab("Weight (pounds)")

# Temperature
spain$dt <- as.character(spain$dt)
year <- (map(strsplit(spain$dt, split = "-"), 1))
year_c <- t(as.data.frame(year))
colnames(year_c) <- "Year"
spain_fin <- cbind(spain, year_c)

spain_fin <- spain_fin %>% filter(spain_fin$Year %in% (c(1750, 1800, 1850, 1900, 1950, 2000, 2013) ))

ggplot(spain_fin, aes(x = AverageTemperature, y = Year, fill = Year)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

#Sunburn

library(plotly)

fig <- fig %>% add_trace(
  type='sunburst',
  ids=species$specy,
  labels=species$genus,
  parents=species$infraorder,
  domain=list(column=1),
  maxdepth=2,
  insidetextorientation='radial'
)
fig

