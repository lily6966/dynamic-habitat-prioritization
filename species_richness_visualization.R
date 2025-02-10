setwd("/Users/lily/Library/CloudStorage/Box-Box/Thesis/SDM/code/")
parent.dir <- "/Users/lily/Library/CloudStorage/Box-Box/Thesis/SDM/code/"
# -----------------------
source.dir <- paste(parent.dir,"source/",sep="")
source(paste(source.dir,"vMR.library.R",sep=""))
source(paste(source.dir,"vMR_baseModels.R",sep=""))
#------------
library(fields)
require(rpart)
library(mgcv)
library(sf)
library(raster)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)
library(verification)
library(ebirdst)
library(fields)
library(gridExtra)
library(tidyverse)
library(sf)
library(raster)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)
library(verification)
library(ebirdst)
library(fields)
library(gridExtra)
library(tidyverse)
library(randomForest)
library(dggridR)
library(mlr)
library(dplyr)
library(ISOweek)
library(rasterVis)
library(exactextractr)

# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection
#seperate species for different rults of counting presence
#common
species <- c("ameavo", "bknsti", "dunlin", "lobdow", "margod","mouplo", "lobcur", "leasan",  "killde", "greyel","whimbr")
#rare
species1 <- c("amgplo",  "lesyel", "pecsan", "uplsan", "uplsan")

scenarios<- c("cdl", "dust", "bbau")
patterns_seasons<-c("^week([6-9]|1[0-9]|2[0-2])\\.tif$", "^week(2[3-9]|3[2-6])\\.tif$",
                    "^week(3[7-9]|4[0-4])\\.tif$",  "^week(4[5-9]|5[0-2]|[1-5])\\.tif$")
stack_list_scenarios<-NULL #create empty list to store raster list for each scenario each season
stack_list_scenarios1<-NULL #create empty list to store raster list for each scenario each season




#create list for common species

for(b in 1:3){
  # Initialize an empty list for each scenario
  raster_stack_list_scenario <- list()
  for(k in 1:4) {
    
    read_mean_raster <- function(folder_path) {
      file_list <- list.files(folder_path, pattern =  patterns_seasons[k], full.names = TRUE) #pattern reference "^week(3[8-9]|4[0-1])\\.tif$"
      raster_stack <- stack(file_list)
      mean_raster<- calc(raster_stack, fun = mean)
      return(mean_raster)
    }
    
    # Define a vector to store folder paths
    folder_paths <- character(length(species))
    
    # Populate folder_paths with folder names
    for (i in 1:length(species)) {
      folder_paths[i] <- paste0(species[i], "_abundance_10year_", scenarios[b])
    }
    
    # Create an empty list to store raster bricks
    raster_stack_list <- list()
    
    
    # Iterate over each folder path
    for (folder_path in folder_paths) {
      raster_stack_mean <- read_mean_raster(folder_path)
      raster_stack_list[[folder_path]] <- raster_stack_mean
    }
    raster_stack_list_scenario[[paste0("season", k)]]<- raster_stack_list
    
  }
  stack_list_scenarios[[paste0("scenario", b)]] <- raster_stack_list_scenario
} 

# Define a function to count the number of layers with values greater than 0.3
count_over_threshold <- function(x) {
  return(sum(x > 0.30, na.rm = TRUE))
}
#creat list of raster for rare species
for(b in 1:3){
  # Initialize an empty list for each scenario
  raster_stack_list_scenario <- list()
  for(k in 1:4) {
    
    read_mean_raster <- function(folder_path) {
      file_list <- list.files(folder_path, pattern =  patterns_seasons[k], full.names = TRUE) #pattern reference "^week(3[8-9]|4[0-1])\\.tif$"
      raster_stack <- stack(file_list)
      mean_raster<- calc(raster_stack, fun = mean)
      return(mean_raster)
    }
    
    # Define a vector to store folder paths
    folder_paths <- character(length(species1))
    
    # Populate folder_paths with folder names
    for (i in 1:length(species1)) {
      folder_paths[i] <- paste0(species1[i], "_abundance_10year_", scenarios[b])
    }
    
    # Create an empty list to store raster bricks
    raster_stack_list <- list()
    
    
    # Iterate over each folder path
    for (folder_path in folder_paths) {
      raster_stack_mean <- read_mean_raster(folder_path)
      raster_stack_list[[folder_path]] <- raster_stack_mean
    }
    raster_stack_list_scenario[[paste0("season", k)]]<- raster_stack_list
    
  }
  stack_list_scenarios1[[paste0("scenario", b)]] <- raster_stack_list_scenario
} 




count_over_threshold1 <- function(x) {
  return(sum(x > 0, na.rm = TRUE))
}

map_proj <- map_proj <- "ESRI:102003"
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
cv <- read_sf("data/gis-data.gpkg", "cv") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()



# Get the bounding box of the trend_pu geometry
bbox <- st_bbox(occurance_pu)

# Center the trend_pu on the map
center_x <- (bbox$xmin + bbox$xmax) / 2
center_y <- (bbox$ymin + bbox$ymax) / 2

# Transform the bounding box to the same CRS as ne_land
bbox_buffered<-st_buffer(st_as_sfc(bbox), dist=50000)
ne_land_transformed <- st_transform(ne_land, crs = st_crs(bbox_buffered))

# Crop ne_land using the transformed bounding box
ne_land_cropped <- st_crop(ne_land_transformed, bbox_buffered)
#define color pallete for each level of richness
library(colorspace)

# Define base color
base_color <- "Viridis"  # Green color

# Define number of shades
n_shades <- 17

# Generate shades with different lightness levels
color_palette <- rev(sequential_hcl(n = n_shades, h = base_color, l = c(30, 90)))
show_col(color_palette)


for(b in 1:3){
  scenario_name <- paste0("scenario", b)
  
  
  for(k in 1:4) {
    season_name <- paste0("season", k)
    # Apply the function to the RasterStack
    count_raster <- calc(stack(stack_list_scenarios[[scenario_name]][[season_name]]), fun = count_over_threshold)
    count_raster1 <- calc(stack(stack_list_scenarios1[[scenario_name]][[season_name]]), fun = count_over_threshold1)
    count_raster<-count_raster+count_raster1
    occurance_pu <- NULL
    
    # get the buffered checklists for a given year and extract values and calculate mean
    regions <- pu_multi 
    occurance_pu <- exact_extract(count_raster, regions, progress = FALSE) %>% 
      map_dfr(~ tibble(occurance_mean = mean(.$value, na.rm = TRUE))) %>% 
      # join to lookup table to get id
      bind_cols(regions, .)
    
    occurance_pu$occurance_mean<-round(occurance_pu$occurance_mean) #remove non-integer occurance
    occurance_pu <- occurance_pu %>%
      mutate(occurance_mean = ifelse(is.nan(occurance_mean), 0, occurance_mean))
    occurance_pu$occurance_mean <- factor(occurance_pu$occurance_mean, levels = sort(unique(occurance_pu$occurance_mean)))
    
    occurance<- ggplot() +
      geom_sf(data = ne_land_cropped, fill = "#dddddd", color = "#888888", lwd = 0.5)+
      geom_sf(data = occurance_pu, aes(fill = occurance_mean), color = "#888888")  +  # Add base map
      #geom_sf(data = occurance_pu[occurance_pu$occurance_mean == 0, ], fill = "#fb8072",  color = "#888888") +  # Add white layer for zero value
      scale_fill_manual(values = (setNames(color_palette, levels(occurance_pu$occurance_mean))), na.value = "transparent") +
      labs(title = paste0("scenario",b, "_season", k), caption = NULL) +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      # Center the trend_pu on the map
      coord_sf(crs = st_crs(ne_land_cropped), xlim = c(center_x - 150000, center_x + 150000), ylim = c(center_y - 350000, center_y + 350000)) 
    # Define the directory path
    tif_dir <- file.path("occurance2")
    
    # Check if the directory exists, create it if not
    if (!dir.exists(tif_dir)) {
      dir.create(tif_dir, recursive = TRUE)
    }
    
    ggsave(paste0("occurance2/",  "scenario",b, "season", k, ".pdf"), plot = occurance, width = 8, height = 6, units = "in")
    
    
    
  }
}

