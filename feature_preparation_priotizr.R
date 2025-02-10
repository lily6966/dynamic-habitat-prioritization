library(prioritizrdata)
library(prioritizr)
library(sf)
library(terra)
library(vegan)
library(cluster)
library(gurobi)
library(tidyverse)
library(sf)
library(raster)
library(gridExtra)
library(tidyverse)
library(raster)
library(tiff)
library(viridis)
library(dggridR)
library(ggplot2)
library(dplyr)
library(sp)
library(RColorBrewer)
library(terra)

library(withr)



wd <- ("/Users/lily/Library/CloudStorage/Box-Box/Thesis/SDM/code/")
setwd(wd)

# # Define a function to read raster files from one folder
# read_raster_files <- function(folder_path, crs = NULL) {
#   file_list <- list.files(folder_path, pattern = "^week[1-4]\\.tif$", full.names = TRUE)
#   
#   if (is.null(crs)) {
#     raster_stack <- stack(file_list)
#   } else {
#     raster_stack <- stack(file_list, crs = crs)
#   }
#   
#   return(raster_stack)
# }

#-----------------------------------Create birds abundance feature for the month you want---------------
# Define a function to read raster files from one folder
read_raster_files <- function(folder_path) {
  file_list <- list.files(folder_path, pattern = "^week3[4-7]\\.tif$", full.names = TRUE)
  raster_stack <- stack(file_list)
  return(raster_stack)
}

species <- c("ameavo", "bknsti", "dunlin", "lobdow", "mouplo", "lobcur", "leasan", "lesyel", "killde", "greyel")

# Define a vector to store folder paths
folder_paths <- character(length(species))

# Populate folder_paths with folder names
for (i in 1:length(species)) {
  folder_paths[i] <- paste0(species[i], "_abundance_10year_dust")
}

# Create an empty list to store raster bricks
raster_stack_list <- list()


# Iterate over each folder path
for (folder_path in folder_paths) {
  raster_stack <- read_raster_files(folder_path)
  raster_stack_list[[folder_path]] <- raster_stack
}

# Combine all raster bricks into a single SpatRaster
birds <- stack(raster_stack_list)
birds_spat<-rast(birds) #go to line 369  to add other features 


#read the shapefile of land use scenarios generated in python 
crop_coords_dust <- read_csv("data/pland-elev_nort_east_prediction-surface_greyel_cdl.csv")
#pland-elev_nort_east_prediction-surface_greyel_cdl, pland_pmp_bbau_lucas
#transform spatial dataframe to spatial features of a normal year
land_cover_pmp_dust <- crop_coords_dust %>% 
  filter(year == 2014) %>% #change to 2054 for climate change scenarios
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"),crs =  4326) 
land_cover_dust<-land_cover_pmp_dust[, "geometry"]
# Calculate the convex hull of land_cover_pmp_dust
convex_hull1 <- st_convex_hull(land_cover_pmp_dust)

# Buffer the convex hull
land_cover_pmp <- st_buffer(convex_hull1, dist = 10000)

land_cover_pmp<-st_transform(land_cover_pmp, crs=st_crs(pu_multi))


#create planning unit file

# #load cv boundary and crop elev 
# cv <-read_sf("data/gis-data.gpkg", "cv") %>% 
#   # project to the cdl projection
#   st_transform(crs = 4326)
# 
# # crop, buffer cv_cov by 10 km to provide a little wiggly room
# cv_extent <- cv %>% 
#   st_buffer(dist = 10000) 

#---------- Calculate the convex hull of land_cover_pmp_dust to create hexagon---------
convex_hull <- st_convex_hull(land_cover_dust)

# Buffer the convex hull
cv_extent <- st_buffer(convex_hull, dist = 10000)
# Combine all polygons into one and simply
combined_polygon <- st_union(cv_extent)
plot(combined_polygon)

valid <- st_is_valid(combined_polygon)

# If not valid, attempt to fix issues
if (!valid) {
  combined_polygon <- st_make_valid(combined_polygon)
}
set.seed(100)
#creating hexagon
pts <- st_sample(combined_polygon, 3000) %>% 
  st_sf(as.data.frame(st_coordinates(.)), geometry = .) %>% #converts the sampled points into an sf object with a geometry column representing the spatial points. 
  rename(lat = Y, lon = X)
dggs <- dgconstruct(spacing=10)
pts$cell <- dgGEO_to_SEQNUM(dggs, pts$lon, pts$lat)$seqnum #his assigns a cell identifier to each point in the pts dataset. 
#It uses the dgGEO_to_SEQNUM() function to convert the latitude and longitude coordinates of the points to cell identifiers (seqnum) based on the DGGS grid.
hexagon <- dgcellstogrid(dggs, unique(pts$cell)) %>%  #unique(pts$cell): This extracts the unique cell identifiers from the pts dataset.
  st_as_sf() #This converts the hexagonal grid polygons to an sf object.
#The resulting pts and hexagon objects will contain the sampled points and the hexagonal grid polygons

ggplot() +
  geom_sf(data = hexagon, colour = "red", fill = NA ) +
  theme_bw()
st_write(hexagon, "planning_uni.shp", rewrite=TRUE)

cost <- st_read("/Users/lily/Library/CloudStorage/Box-Box/calvin/cost_water_dust.shp")
cost <- st_read("c:/Users/liliy/Box/calvin/cost_water_dust.shp")
# Perform the join without Z and M coordinates
cost <- st_transform(cost, crs = st_crs(pu_multi))
cost_pu <- st_join(pu_multi, cost, join = st_nearest_feature)
pu_multi$cost_dust<-cost_pu$cost_of_wa
cost1 <- st_read("c:/Users/liliy/Box/calvin/cost_water_his.shp")
cost1 <- st_transform(cost1, crs = st_crs(cost_pu))
cost_pu1 <- st_join(cost_pu, cost1, join = st_nearest_feature)
pu_multi$cost_his<-cost_pu1$cost_of_wa.y
cost2 <- st_read("c:/Users/liliy/Box/calvin/cost_water_bbau.shp")
cost2 <- st_transform(cost2, crs = st_crs(cost_pu))
cost_pu2 <- st_join(cost_pu1, cost2, join = st_nearest_feature)
pu_multi$cost_bbau<-cost_pu2$cost_of_wa

cost <- st_transform(cost, crs = st_crs(hexagon))
cost_v<-vect(cost[,-which(names(cost) == "prmname")])#transfrom to spatial vector to map
CA <- read_sf("data/gis-data.gpkg", "ne_land") %>% 
  # project to the cdl projection
  st_transform(crs =  st_crs(hexagon))

p

num_colors <- 10
color_palette <- colorRampPalette(brewer.pal(9, "Blues"))(num_colors)

# Create breaks for color bins
breaks <- seq(min(cost$node_dual), max(cost$node_dual), length.out = num_colors + 1)
# Plot the sf object with color-coding based on the "node_dual" column

# Plot the sf object with color-coding based on the "node_dual" column
plot(vect(hexagon), lwd=1, border="dark red")
plot(CA, col="light gray", lwd=4, border="dark gray")
plot(cost_v,  col = color_palette, legend=TRUE, pch=16, cex = 1, breaks = breaks, add=TRUE, plg=list(x=-117.27, y=41.54))
lines(hexagon)

# Remove Z and M coordinates from hexagon and cost objects
hexagon_no_zm <- st_zm(hexagon, drop = TRUE)
cost_no_zm <- st_zm(cost, drop = TRUE)

# Perform the join without Z and M coordinates
cost_pu <- st_join(hexagon_no_zm, cost_no_zm, join = st_nearest_feature)


#cost_pu <- st_join(hexagon, cost, join = st_nearest_feature)

# # Extract coordinates from the sp object, same me as line82
# coords <- st_coordinates(cost_pu)
# # Create a data frame with longitude and latitude columns
# cost_df <- data.frame(cost=cost_pu$node_dual,
#   longitude = coords[, 1],
#   latitude = coords[, 2]
# )



d <- data.frame(geom(vect(cost_pu))[,c("x", "y")], as.data.frame(vect(cost_pu)))
# Set the number of nearest neighbors to consider (e.g., 5)

# Write the updated sf object to a new shapefile
st_write(cost_pu, "cost_pu_filled.shp")

cvpia<-st_read("data/CVPIA_Refuges/CVPIA_Refuges_[ds2636].shp") 
cvpia <- st_transform(cvpia, crs = st_crs(hexagon))
sagbi<-st_read("data/SAGBI/sagbi_shape/sagbi_shape.shp")
sagbi <- st_transform(sagbi, crs = crs(birds))


# Transform sagbi to the CRS of the first layer in the raster stack
sagbi <- st_transform(sagbi, crs = crs(birds[[1]]))

# Choose a specific layer from the raster stack (e.g., the first layer)
birds_layer <- birds[[1]]

# Convert sagbi to a raster based on the sagbi field and the spatial extent and resolution of birds_layer
sagbi_raster <- rasterize(sagbi, birds_layer, field = "sagbi")

# Define the file path where you want to save the raster
output_path <- "data/sagbi.tif"

# Save the raster to the specified file path
writeRaster(sagbi_raster, filename = output_path, format = "GTiff", overwrite = TRUE)


birds1<-stack(birds, sagbi_raster)
birds_spat1<-rast(birds1)

# Plot hexagon geometries with cvpia
ggplot() +
  geom_sf(data = hexagon, fill = "lightblue") +
  geom_sf(data = cvpia, color = "red", fill = NA) +
  theme_minimal()
#locked_out <- sample(c(TRUE, FALSE), nrow(cost_pu), replace = TRUE)
pu <- hexagon%>% 
  transmute(id=row_number()) %>% 
  mutate (cost = cost_pu$node_dual, locked_in = st_intersects(geometry, cvpia)
         
          )
#Convert the list to logical (TRUE/FALSE) indicating intersection
pu$locked_in <- sapply(pu$locked_in, function(x) length(x) > 0)
#cast polygon to multipolygon
pu_multi <- st_cast(pu, "MULTIPOLYGON")







##get the hexagon resolution
st_area(hexagon)

crops<-c("rice", "corn", 'pasture', 'alfalfa', 'grain', "perennial", "field_row")

multipliers<-c(rice=469, corn=30, alfalfa=468, pasture=298, perennial=2796, grain=141, field_row=1765)





pu_multi$cost_bbau

# Define a function to calculate the land cost
calculate_land_cost <- function(cell, land_cover_data, multipliers) {
  # Check if points in land_cover_data are within the multipolygon cell
  is_within <- sapply(st_intersects(land_cover_data, cell), function(x) length(x) > 0)
    
    # Extract the rice value from land_cover_data for the corresponding hexagon cell
    crop_values <- na.omit(land_cover_pmp[is_within, crops])  # Assuming columns 6 to 14 contain the crop values
    crop_values<-st_drop_geometry(crop_values)
    # Extract the values from land_cover_pmp_dust for the corresponding hexagon cell
    # Multiply each column by corresponding multiplier
    scaled_values <- crop_values * multipliers[names(crop_values)]
    
    # Calculate the mean of each row
    land_values <- colMeans(scaled_values)
    # Calculate the sum of values multiplied by their corresponding multipliers
    land_cost <- sum(land_values)*23715
    return(land_cost)
  
  
}  
  
# Initialize an empty vector to store land costs
pu_multi$land_cdl <- rep(NA, nrow(pu_multi))

# Iterate over each cell in pu_multi to calculate land cost
for (i in 1:nrow(pu_multi)) {
  cell <- pu_multi[i, "geometry"]
  pu_multi$land_cdl[i] <- calculate_land_cost(cell, land_cover_pmp,  multipliers)
}

#pu_multi$cost_water<-pu_multi$cost*7783 #convert price $1000/kAF to $ by multiply 10cm depth and 23722 acres per cell

# Create a directory
dir.create("data/pu_multi", showWarnings = FALSE, recursive = TRUE)



# Write the spatial dataframe to a GeoPackage file
st_write(pu_multi, "data/pu_multi/pu_multi_cropareas_cdl1.shp")
pu_multi<-st_read("data/pu_multi/pu_multi_processed.shp")
pu_multi<-st_transform(pu_multi, crs=crs(birds_spat))
#--------------calculate farmlands------------------
# Define a function to calculate the carbon sequestion potential
calculate_cropareas <- function(cell, land_cover_data) {
  # Check if points in land_cover_data are within the multipolygon cell
  is_within <- sapply(st_intersects(land_cover_data, cell), function(x) length(x) > 0)
  
  # Extract the rice value from land_cover_data for the corresponding hexagon cell
  crop_areas <- na.omit(land_cover_pmp[is_within, crops])  # Assuming columns 6 to 14 contain the crop values
  crop_areas<-st_drop_geometry(crop_areas)
  # Extract the values from land_cover_pmp_dust for the corresponding hexagon cell
  # Multiply each column by corresponding multiplier
  #cropacres <- crop_areas*23722 # 0.58t/acre/year convert square meters to acre and divieded by 12 months
  
  # Calculate the mean of each row
  cropacres_values <- sum(colMeans(crop_areas))
  
  return(cropacres_values)
  
  
}  

# Initialize an empty vector to store land costs
pu_multi$cropars_cdl <- rep(NA, nrow(pu_multi))

# Iterate over each cell in pu_multi to calculate land cost
for (i in 1:nrow(pu_multi)) {
  cell <- pu_multi[i, "geometry"]
  pu_multi$cropars_cdl[i] <- calculate_cropareas(cell, land_cover_pmp)
}


summary(pu_multi$cropareas)

# Assuming pu_multi is your simple feature collection and it has a column named 'cropareas'

# Create the histogram
hist_data <- hist(pu_multi$cropareas,
                  breaks = 100,  # You can adjust the number of breaks as needed
                  plot = FALSE)  # Create the histogram object without plotting

# Plot the histogram with custom limits and labels
plot(hist_data,
     col = "lightblue",  # Set the color of the bars
     xlim = c(0, 1),     # Set the x-axis limit
     ylim = c(0, 40),   # Set the y-axis limit
     xlab = "Cropland coverage (%)",  # X-axis label
     ylab = "Frequency",                         # Y-axis label
     main = "Feature representation by existing protected areas")  # Main title

#calculate water use based on flooding half of the farmlands in each unit

pu_multi$water_costk<-pu_multi$cropars*23715*0.3281/2*pu_multi$cost_dust/1000#($)convert price $1000/kAF to $ by multiply 10cm depth and 23722 acres per cell
pu_multi$flodd_cdl <-pu_multi$cropars_cdl*23715/2/1000 #in thousand acres
pu_multi$carbon_cdl<-pu_multi$cropars_cdl*23715*0.4*0.58/2
pu_multi$carbSeq<-pu_multi$carbon
#---------------create a raster feature of carbon seqestration---------------------------------------------------
# Define a function to calculate the carbon sequestion potential
calculate_carbSeq <- function(cell, land_cover_data) {
  # Check if points in land_cover_data are within the multipolygon cell
  is_within <- sapply(st_intersects(land_cover_data, cell), function(x) length(x) > 0)
  
  # Extract the rice value from land_cover_data for the corresponding hexagon cell
  crop_areas <- na.omit(land_cover_pmp[is_within, crops])  # Assuming columns 6 to 14 contain the crop values
  crop_areas<-st_drop_geometry(crop_areas)
  # Extract the values from land_cover_pmp_dust for the corresponding hexagon cell
  # Multiply each column by corresponding multiplier
  carbSeq <- colMeans(crop_areas)*0.58*23722*0.4 # 0.58t/hactare/year convert square acres to hactare
  
  # Calculate the mean of each row
  carbSeq_values <- sum(carbSeq)
  
  return(carbSeq_values)
  
  
}  

# Initialize an empty vector to store land costs
pu_multi$carbSeq <- rep(NA, nrow(pu_multi))

# Iterate over each cell in pu_multi to calculate land cost
for (i in 1:nrow(pu_multi)) {
  cell <- pu_multi[i, "geometry"]
  pu_multi$carbSeq[i] <- calculate_carbSeq(cell, land_cover_pmp)
}
birds_layer<-birds[[1]]
# Convert sagbi to a raster based on the sagbi field and the spatial extent and resolution of birds_layer
carbSeq_raster <- rasterize(pu_multi, birds_layer, field = "carbSeq")
# Define the file path where you want to save the raster
output_path <- "data/carbSeq_dust1.tif"

# Save the raster to the specified file path
writeRaster(carbSeq_raster, filename = output_path, format = "GTiff", overwrite = TRUE)

# nomalize carbSeq_raster 

# Define the new minimum and maximum values for normalization
new_min <- 0  # New minimum value
new_max <- 1  # New maximum value

# Normalize the raster values
normalized_carbSeq <- (carbSeq_raster - minValue(carbSeq_raster)) / (maxValue(carbSeq_raster) - minValue(carbSeq_raster)) * (new_max - new_min) + new_min

# Print the summary of the normalized raster
print(normalized_carbSeq)

# Normalize the raster values
normalized_sagbi <- (sagbi_raster - minValue(sagbi_raster)) / (maxValue(sagbi_raster) - minValue(sagbi_raster)) * (new_max - new_min) + new_min

birds3<-stack(birds, normalized_carbSeq, normalized_sagbi)
birds_spat3<-rast(birds3) # go to line 392 to modify cost to line to prioritize







birds2<-stack(birds1, carbSeq_raster)
birds_spat2<-rast(birds2)

pu_multi<-st_read("data/pu_multi/pu_multi_cost.shp")
pu_multi$locked_in<-as.logical(pu_multi$locked_in)


pu_multi<-st_transform(pu_multi, crs=crs(birds_spat))

#detemine cost column based on month
pu_multi$costk<-pu_multi$land_cost/1000+pu_multi$cost/1000

# build problem
p1 <-
  problem(pu_multi, birds_spat3, cost_column = "costk") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = 0.05) %>%
  add_relative_targets(0.3) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions() 

print(presolve_check(p1))

# create vector to store plot names
n <- c()

# create empty list to store solutions
s <- c()

if (require("gurobi")) {
  p2 <- p1 %>% add_gurobi_solver(verbose = FALSE)
  n <- c(n, "gurobi")
  s <- c(s, solve(p2))
}

# solve problem
s1 <- solve(p2)
# plot map of prioritization
plot(
  st_as_sf(s1[, "solution_1"]), main = "Jan(0.005)",
  pal = c("grey90", "darkgreen")
)

# create column for making a map of the prioritization
s1$map_1 <- case_when(
  s1$locked_in > 0.5 ~ "CVPIA",
  s1$solution_1 > 0.5 ~ "Flooded dynamic Habitat",
  TRUE ~ "other"
)

# plot map of prioritization
plot(
  s1[, "map_1"], pal = c("skyblue", "salmon", "grey90"),
  main = NULL, key.pos = 1
)



