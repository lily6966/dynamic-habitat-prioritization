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
library(topsis)

setwd("your_wd")

#------read the carbon sequestration feature for the scenario you want construct prioritization-----
carbpath<-"data/carbSeq_cdl1.tif" #carbSeq_cdl1.tif, carbSeq_dust1.tif, or carbSeq_bbau1.tif
carbSeq_raster<-raster(carbpath)
# Define the new minimum and maximum values for normalization
new_min <- 0  # New minimum value
new_max <- 1  # New maximum value
# Normalize the raster values
normalized_carbSeq <- (carbSeq_raster - minValue(carbSeq_raster)) / (maxValue(carbSeq_raster) - minValue(carbSeq_raster)) * (new_max - new_min) + new_min

# Define the file path where you saved the raster
sagbi_path <- "data/sagbi.tif"
sagbi_raster <- raster(sagbi_path)
# Normalize the raster values
normalized_sagbi <- (sagbi_raster - minValue(sagbi_raster)) / (maxValue(sagbi_raster) - minValue(sagbi_raster)) * (new_max - new_min) + new_min
# Write the spatial dataframe of planning unit from a GeoPackage file
pu_multi<-st_read("data/pu_multi/pu_multi_cropareas_cdl.shp")
pu_multi<-st_transform(pu_multi, crs=crs(birds_spat))
pu_multi$locked_in1<-as.logical(pu_multi$locked_n)
# Define a function to read raster files from one folder
patterns<-c("^week([1-5])\\.tif$", "^week([6-9])\\.tif$","^week(1[0-4])\\.tif$","^week(1[5-8])\\.tif$",
            "^week(1[9]|2[0-2])\\.tif$", "^week(2[3-6])\\.tif$","^week(2[7-9]|3[0-1])\\.tif$", "^week(3[2-6])\\.tif$",
            "^week(3[7-9]|4[0])\\.tif$", "^week(4[1-4])\\.tif$", "^week(4[5-8])\\.tif$" ,"^week(4[8-9]|[1-2])\\.tif$")

d<-c(50,40,50,40,40,40,50,50,40,40,40,40)
penalties <- c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,0.001, 0.001)

tar<-c(0.4, 0.4, 0.4, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4)
pu_multi$water_costk<-pu_multi$cropars *23715*0.3281/2*pu_multi$cost_hs/1000
#detemine cost column based on month

cost<-list(pu_multi$water_costk/4, pu_multi$water_costk/2, pu_multi$water_costk/2)
index_cost<-c(1,2,2,3,3,3,3,3,3,3,2,1)

result_plot_group<-NULL

# generate boundary length data for the planning units
pu_bd <- boundary_matrix(pu_multi)

# manually re-scale the boundary length values
pu_bd <- rescale_matrix(pu_bd)
pu_multi$zeros<-0

cost_cov<-NULL
tc<-NULL
tc1<-NULL
scenarios<-c("cdl", "dust", "bbau")
#prioritize over each month
#mmm

for(b in 1:1){
  conservation_cost<-NULL
  for(k in 1:12) {
    
    # Create the list of targets
    #t <- c(rep(0.30, d[k]))
 
    #detemine cost column based on month
    pu_multi$costk<-unlist(cost[index_cost[k]])
    
    
    read_raster_files <- function(folder_path) {
      file_list <- list.files(folder_path, pattern =  patterns[k], full.names = TRUE) #pattern reference "^week(3[8-9]|4[0-1])\\.tif$"
      raster_stack <- stack(file_list)
      return(raster_stack)
    }
    
    species <- c("ameavo","amgplo", "bknsti", "dunlin", "lobdow", "margod","mouplo", "lobcur", "leasan", "lesyel", "killde", "greyel","pecsan", "uplsan","shbdow", "whimbr")
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
      raster_stack <- read_raster_files(folder_path)
      raster_stack_list[[folder_path]] <- raster_stack
    }
    # Combine all raster bricks into a single SpatRaster
    birds <- stack(raster_stack_list)
    birds3<-stack(birds, normalized_carbSeq, normalized_sagbi)
    birds_spat6<-rast(birds3) # go to line 392 to modify cost to line to prioritize
    
    pu_multi<-st_transform(pu_multi, crs=crs(birds_spat6))
    pu_multi$pa <- round(pu_multi$locked_in)

    t<-tar[k]
    # build problem
    p3 <-
      problem(pu_multi, birds_spat6,cost_column = "costk") %>%
      add_min_set_objective() %>%
      add_boundary_penalties(penalty = penalties[k]) %>%
      add_relative_targets(t)  %>%
      add_locked_in_constraints("locked_in1") %>%
      add_binary_decisions() %>%
      add_default_solver(gap = 0)



    # create vector to store plot names
    n <- c()

    # create empty list to store solutions
    s <- c()

    if (require("gurobi")) {
      p4 <- p3 %>% add_gurobi_solver(verbose = TRUE)
      n <- c(n, "gurobi")
      s <- c(s, solve(p4))
    }

    # solve problem
    s2 <- solve(p4)


    # create column for making a map of the prioritization
    s2$map_1 <- case_when(
      s2$locked_in1 > 0.5 ~ "Protected areas",
      s2$solution_1 > 0.5 ~ "Flooded dynamic Habitat",
      TRUE ~ "other"
    )


    # Start PDF device
    pdf(paste0("month",k,"_varyingT_", scenarios[b], ".pdf"))

    # plot map of prioritization
    plot(
      s2[, "map_1"], pal = c("salmon", "grey90","skyblue"),
      main = NULL, key.pos = 1
    )
    s2_costk<-eval_cost_summary(p4, s2[, "solution_1"])$cost
    flooded_areas_p2 <- s2[, "solution_1"] %>%
      st_drop_geometry() %>%
      bind_cols(pu_multi) %>%
      filter(solution_1 == 1) %>%
      summarise(sum_flodd_r = sum(fldd_cd)) %>%
      pull(sum_flodd_r)
    # Add text with the cost value to the plot
    # Get plot extent for appropriate text placement
    plot_extent <- par("usr")  # Get the current plot dimensions (x1, x2, y1, y2)
    
    text(x = plot_extent[1] + (plot_extent[2] - plot_extent[1]) * 1,  # 10% from the left
         y = plot_extent[3] + (plot_extent[4] - plot_extent[3]) * 0.8,  # 90% from the bottom
         labels = paste0("Cost =", s2_costk,"k","\nFlooded areas=", flooded_areas_p2), pos = 2)
    dev.off()
    # calculate cost
    conservation_cost[k] <- eval_cost_summary(p3, s2[, "solution_1"])$cost

    # calculate irreplaceability
    irrep_s1 <- eval_replacement_importance(p3, s2["solution_1"])

    # manually coerce values for planning units not selected in prioritization
    # to NA, so that they are shown in white
    irrep_s1$plot_total <- irrep_s1$rc
    irrep_s1$plot_total[s2$solution_1 < 0.5] <- NA_real_


    # Open a PDF graphics device
    pdf(paste0("month", k, "_importance_", scenarios[b], ".pdf"))

    # Plot the map of overall importance scores
    plot(st_as_sf(irrep_s1[, "plot_total"]), main = "Overall importance")

    # Close the PDF device
    dev.off()
    #calculate feature representation
    tc_pa<-  eval_target_coverage_summary(p1, pu_multi[,'pa'])$relative_held*100
    tc1<-c(tc1, tc_pa)
    tc_s1 <- eval_target_coverage_summary(p1, s1[,'solution_1'])$relative_held*100
    ## visualize representation  (values show percent coverage)
    tc<-c(tc, tc_s1)
#---------------------when don't need penalty evaluation--------------------------
  }
  # Open a PDF graphics device
pdf(paste0("hist_month_rep_s1.pdf"))
hist(
  tc,
  col = "salmon",
  main = "Feature representation by prioritization",
  xlim = c(0, 100),
  xlab = "Percent coverage of features (%)"
)
  # Close the PDF device
dev.off()

pdf(paste0("hist_month_rep_pa.pdf"))
  hist(
    tc1,
    col = "skyblue",
    main = "Feature representation by prioritization",
    xlim = c(0, 100),
    xlab = "Percent coverage of features (%)"
  )
  # Close the PDF device
dev.off()
  # Save to CSV
  #write.csv(conservation_cost, paste0(scenarios[b], "_cost_varyingT.csv"), row.names = FALSE)
}
#-----------------------------------------------------------------------------------------------------
    
    #
    #----------blended approach for detemining penalty, if not skip------------------

    # # define a range of different penalty values
    # #prelim_penalty<-c(0.00001, 0.0001, 0.001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.3)
    # prelim_penalty<-c(0.001, 0.01, 0.1, 0.3, 1, 10, 20, 50)
    #
    # # define a problem without boundary penalties
    # p0 <-
    #   problem(pu_multi, birds_spat6, cost_column = "costk") %>%
    #   add_min_set_objective() %>%
    #   add_relative_targets(t) %>%
    #   add_locked_in_constraints("locked_in") %>%
    #   add_binary_decisions()
    #
    # # generate preliminary prioritizations based on each penalty
    # ## note that we specify a relaxed gap and time limit for the solver
    # prelim_blended_results <- lapply(prelim_penalty, function(x) {
    #   s <-
    #     p0 %>%
    #     add_boundary_penalties(penalty = x, data=pu_bd) %>%
    #     add_default_solver(gap = 0.2, time_limit = 10 * 60) %>%
    #     solve()
    #   s <- data.frame(s = s$solution_1)
    #   names(s) <- with_options(list(scipen = 30), paste0("penalty_", x))
    #   s
    # })
    #
    # # format results as a single spatial object
    # prelim_blended_results <- cbind(
    #   pu_multi, do.call(bind_cols, prelim_blended_results)
    # )
    #
    #
    #
    # # calculate metrics for blended method prioritizations
    # ## note that we use p0 and not p1 so that cost calculations are based
    # ## on the cost values and not zeros
    # blended_metrics <- lapply(
    #   grep("penalty_", names(prelim_blended_results)), function(x) {
    #     x <- prelim_blended_results[, x]
    #     data.frame(
    #       total_cost = eval_cost_summary(p0, x)$cost,
    #       total_boundary_length = eval_boundary_summary(p0, x)$boundary
    #     )
    #   }
    # )
    # blended_metrics <- do.call(bind_rows, blended_metrics)
    # blended_metrics$penalty <- prelim_penalty
    # blended_metrics <- as_tibble(blended_metrics)
    #
    # # create data for plotting
    # result_data1 <-
    #   blended_metrics %>%
    #   ## rename threshold column to value column
    #   rename(value = "penalty") %>%
    #   ## add column with column names that contain candidate prioritizations
    #   mutate(name = grep(
    #     "penalty_", names(prelim_blended_results), value = TRUE, fixed = TRUE
    #   )) %>%
    #   ## add column with labels for plotting
    #   mutate(label = paste("penalty =", value)) %>%
    #   ## add column to keep track prioritizations selected by different methods
    #   mutate(method = "none")
    #
    #
    # # specify prioritization selected by visual method
    # result_data1$method[3] <- "visual"
    #
    # # calculate TOPSIS scores
    # topsis_results <- topsis(
    #   decision =
    #     blended_metrics %>%
    #     dplyr::select(total_cost, total_boundary_length) %>%
    #     as.matrix(),
    #   weights = c(1, 1),
    #   impacts = c("-", "-")
    # )
    #
    #
    # # add column indicating prioritization selected by TOPSIS method
    # result_data1$method[which.max(topsis_results$score)] <- "TOPSIS"
    #
    # # generate ideal prioritization based on cost criteria
    #
    # p0 <-
    #   problem(pu_multi, birds_spat6, cost_column = "costk") %>%
    #   add_min_set_objective() %>%
    #   add_relative_targets(t) %>%
    #   add_locked_in_constraints("locked_in") %>%
    #   add_binary_decisions() %>%
    #   add_default_solver(gap = 0)
    #
    # # solve problem
    # s0 <- solve(p0)
    #
    # # generate ideal prioritization based on spatial fragmentation criteria
    # ## note that any non-zero penalty value would work here,
    # ## so we just use a penalty  of 1
    #
    # p00 <-
    #   problem(pu_multi, birds_spat6, cost_column = "zeros") %>%
    #   add_min_set_objective() %>%
    #   add_boundary_penalties(penalty = penalties[k]) %>%
    #   add_relative_targets(t) %>%
    #   add_locked_in_constraints("locked_in") %>%
    #   add_binary_decisions() %>%
    #   add_default_solver(gap = 0)
    # s00<-solve(p00)
    #
    # # generate problem formulation with costs and boundary penalties for
    # # calculating performance metrics
    # p_metrics <-
    #   problem(pu_multi, birds_spat6, cost_column = "costk") %>%
    #   add_min_set_objective() %>%
    #   add_boundary_penalties(penalty = penalties[k]) %>%
    #   add_relative_targets(t) %>%
    #   add_locked_in_constraints("locked_in") %>%
    #   add_binary_decisions()
    #
    # # calculate performance metrics for ideal cost prioritization
    # s0_metrics <- tibble(
    #   total_cost = eval_cost_summary(p_metrics, s0[, "solution_1"])$cost,
    #   total_boundary_length =
    #     eval_boundary_summary(p_metrics, s0[, "solution_1"])$boundary
    # )
    #
    # # calculate performance metrics for ideal boundary length prioritization
    # s00_metrics <- tibble(
    #   total_cost = eval_cost_summary(p_metrics, s00[, "solution_1"])$cost,
    #   total_boundary_length =
    #     eval_boundary_summary(p_metrics, s00[, "solution_1"])$boundary
    # )
    #
    # # calculate penalty value based on Cohon et al. 1979
    # cohon_penalty <- abs(
    #   (s0_metrics$total_cost - s00_metrics$total_cost) /
    #     (s0_metrics$total_boundary_length - s00_metrics$total_boundary_length)
    # )
    #
    # # round to 5 decimal places to avoid numerical issues during optimization
    # cohon_penalty <- round(cohon_penalty, 5)
    #
    #
    # # generate prioritization using penalty value calculated using Cohon et al. 1979
    # p6 <-
    #   problem(pu_multi, birds_spat6, cost_column = "costk") %>%
    #   add_min_set_objective() %>%
    #   add_boundary_penalties(penalty = cohon_penalty) %>%
    #   add_relative_targets(t) %>%
    #   add_locked_in_constraints("locked_in") %>%
    #   add_binary_decisions()
    #
    # # solve problem
    # s6 <- solve(p6)
    #
    # # add new row with data for prioritization generated following Cohon et al. 1979
    # result_data1 <- bind_rows(
    #   result_data1,
    #   tibble(
    #     total_cost = eval_cost_summary(p6, s6[, "solution_1"])$cost,
    #     total_boundary_length =
    #       eval_boundary_summary(p6, s6[, "solution_1"])$boundary,
    #     value = cohon_penalty,
    #     name = paste0("penalty_", cohon_penalty),
    #     label = paste0("Penalty = ",  cohon_penalty),
    #     method = "Cohon"
    #   )
    # )
    #
    #
    #
    # # create plot to visualize trade-offs and show selected candidate prioritization
    # result_plot1 <-
    #   ggplot(
    #     data =
    #       result_data1 %>%
    #       mutate(vjust = if_else(method == "Cohon", -1, 0.5)),
    #     aes(x = total_boundary_length, y = total_cost, label = label)
    #   ) +
    #   geom_line() +
    #   geom_point(aes(color = method), size = 3) +
    #   geom_text(aes(vjust = vjust, color = method), hjust = -0.1) +
    #   scale_color_manual(
    #     name = "Method",
    #     values = c(
    #       "visual" = "#984ea3",
    #       "none" = "#000000",
    #       "TOPSIS" = "#e41a1c",
    #       "Cohon" = "#377eb8"
    #     )
    #   ) +
    #   xlab("Total boundary length of prioritization") +
    #   ylab("Total cost of prioritization") +
    #   scale_x_continuous(expand = expansion(mult = c(0.05, 0.4)))
    #
    # # Open a PDF graphics device
    # pdf(paste0("month", k, "_penalty_newcost_", scenarios[b], ".pdf"))
    #
    # # Save the ggplot object
    # print(result_plot1)
    #
    # # Close the PDF device
    # dev.off()



#   }
# 
# # Save to CSV
# write.csv(conservation_cost, paste0(scenarios[b], "_newcost_varyingT.csv"), row.names = FALSE)
# }


#mmm
#----------each month each zone----------
#---------------create a raster pu of month 5-10---------------------------------------------------

pu_multi$water_costk<-pu_multi$cropars_cdl *23715*0.3281/2*pu_multi$cost_his/1000
#detemine cost column based on month


# calculate costs
pu_multi$costkMay <- pu_multi$water_costk/2
pu_multi$costkDec <- pu_multi$water_costk/4
pu_multi$costkFeb <- pu_multi$water_costk/2

# Convert  to a raster 
birds_layer <- birds[[1]]
pu_multi<-st_transform(pu_multi, crs=crs(birds_layer))
costkMay <- rasterize(pu_multi, birds_layer, field = "costkMay")
costkDec <- rasterize(pu_multi, birds_layer, field = "costkDec")
costkFeb <- rasterize(pu_multi, birds_layer, field = "costkFeb")

# create planning unit data with costs
pu_multi_zones <- stack(
  costkDec, costkFeb, costkFeb, costkFeb, 
  costkMay, costkMay, costkMay, costkMay, costkMay, costkMay,
  costkFeb, costkDec
)
names(pu_multi_zones) <- c("Jan", "Feb", "Mar", "Apr",
                           "May", "Jun", "Jul", "Aug", "Sep", "Oct",
                           "Nov", "Dec")
pu_multi_zones <- rast(pu_multi_zones)



# Define a function to read raster files from one folder
patterns<-c("^week([1-5])\\.tif$", "^week([6-9]|1[0])\\.tif$","^week(1[0-4])\\.tif$","^week(1[5-9])\\.tif$",
            "^week(1[9]|2[0-3])\\.tif$", "^week(2[3-7])\\.tif$","^week(2[7-9]|3[0-1])\\.tif$", "^week(3[2-6])\\.tif$",
            "^week(3[7-9]|4[0-1])\\.tif$", "^week(4[1-5])\\.tif$", "^week(4[5-9])\\.tif$" ,"^week([1-3]|4[8-9])\\.tif$")
#"^week(1|4[9]|5[0-2])\\.tif$" replace for dec for bbau and dust
features<-list()
#mmm
for(k in 1:12) {
 
  read_raster_files <- function(folder_path) {
    file_list <- list.files(folder_path, pattern =  patterns[k], full.names = TRUE) #pattern reference "^week(3[8-9]|4[0-1])\\.tif$"
    raster_stack <- stack(file_list)
    return(raster_stack)
  }
  
  species <- c("ameavo","amgplo", "bknsti", "dunlin", "lobdow", "margod","mouplo", "lobcur", "leasan", "lesyel", "killde", "greyel","pecsan", "uplsan","shbdow", "whimbr")
  
  # Define a vector to store folder paths
  folder_paths <- character(length(species))
  
  # Populate folder_paths with folder names
  for (i in 1:length(species)) {
    folder_paths[i] <- paste0(species[i], "_abundance_10year_cdl")
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
  #add otheer features
  birds3<-stack(birds, normalized_carbSeq, normalized_sagbi)
  birds_spat3<-rast(birds3) # go to line 392 to modify cost to line to prioritize
  names(birds_spat3) <-1:82
  features[[k]]<-birds_spat3
}


#targets with each zone contribute equally
t2 <- matrix(NA, ncol = 12, nrow = 82)
t2[, 5:10] <- 0.025
t2[, 1:4] <- 0.025
t2[, 11:12] <- 0.025


# build multi-zone problem
p2 <-
  problem(pu_multi_zones, zones(
    features[[1]],features[[2]], features[[3]],features[[4]],
    features[[5]],features[[6]], features[[7]],features[[8]],
    features[[9]],features[[10]], features[[11]],features[[12]]
  )) %>%
  add_min_set_objective() %>%
  add_relative_targets(t2) %>%
  add_binary_decisions()

# create vector to store plot names
n <- c()

# create empty list to store solutions
s <- c()

if (require("gurobi")) {
  p2 <- p2 %>% add_gurobi_solver(time_limit = 9000, verbose = TRUE)
  n <- c(n, "gurobi")
  s <- c(s, solve(p2))
}


# solve problem
s2 <- solve(p2)






# create targets that all the zone don't need to contribute the save
t4 <- tibble(
  feature = names(features[[1]]),
  zone = list(names(pu_multi_zones))[rep(1, 82)],
  target = rep(0.025, 82),
  type = rep("relative", 82)
)

# print targets
print(t4)

# create problem
p4 <-
  problem(
    pu_multi_zones,
    zones(
      features[[1]],features[[2]], features[[3]],features[[4]],
      features[[5]],features[[6]], features[[7]],features[[8]],
      features[[9]],features[[10]], features[[11]],features[[12]],
      zone_names = names(pu_multi_zones),
      feature_names = names(features[[1]])
    )
  ) %>%
  add_min_set_objective() %>%
  add_manual_targets(t4) %>%
  add_binary_decisions() 

# create vector to store plot names
n <- c()

# create empty list to store solutions
s <- c()

if (require("gurobi")) {
  p4 <- p4 %>% add_gurobi_solver(time_limit = 9000, verbose = TRUE)
  n <- c(n, "gurobi")
  s <- c(s, solve(p4))
}

# solve problem
s4 <- solve(p4)


# plot solution
plot(category_layer(s4), main = "solution", axes = FALSE)


# calculate feature representation
r4 <- eval_feature_representation_summary(p4, s4)
print(r4)


# #lockout planning units if needed
# # create multi-layer raster with locked in units
# 
# locked_inout_raster <- rasterize(pu_multi, birds_layer, field = "locked_in")
# 
# # Create an empty stack
# stacked_rasters <- stack()
# 
# # Loop through each raster file and stack it
# for (i in 1:12) {
#   stacked_rasters <- stack(stacked_rasters, locked_inout_raster)
# }
# 
# locked_raster<-rast(stacked_rasters)
# 
# p5 <- p4 %>% add_locked_out_constraints(locked_raster)
# p6<-  p2 %>% add_locked_out_constraints(locked_raster)
# # preview locked data
# # solve problem
# s5 <- solve(p5)
# 
# # plot solution
# plot(category_layer(s5), main = "solution", axes = FALSE)
# # solve problem
# s6 <- solve(p6)
# 
# # plot solution
# plot(category_layer(s6), main = "solution", axes = FALSE)

# Define a function to create the matrix
create_month_matrix <- function(n) {
  # Initialize an empty matrix with dimensions n x n
  mat <- matrix(0, nrow = n, ncol = n)
  
  # Loop through each row and column
  for (i in 1:n) {
    for (j in 1:n) {
      # Set diagonal elements to 1
      if (i == j) {
        mat[i, j] <- 1
      }
      # Set elements to 1 for pairs of consecutive months
      if (j == i + 1 || j == i - 1) {
        mat[i, j] <- 1
      }
    }
  }
  
  # Return the resulting matrix
  return(mat)
}

# Call the function to create a matrix with 12 rows and columns
month_matrix <- create_month_matrix(12)
#modify pair of Jan and Dec

month_matrix[1,12]<-1
month_matrix[12,1]<-1

colnames(month_matrix) <- names(pu_multi_zones)
rownames(month_matrix) <- colnames(month_matrix)
# Print the resulting matrix
# Convert the matrix to a data frame in long format

library(reshape2)
month_df <- melt(month_matrix)

# Plot the matrix
ggplot(month_df, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = factor(value)), color = "white") +
  geom_text(aes(label = value), color = "black") +
  scale_fill_manual(values = c("white", "blue"), name = "Value") +
  labs(x = "Month", y = "Month", title = "Month Matrix") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#don't do locked out
p7<-p4 %>% add_boundary_penalties(penalty = 0.001, zones = month_matrix)
p8<-p2 %>% add_boundary_penalties(penalty = 0.001, zones = month_matrix)
# preview locked data
# solve problem
s7 <- solve(p7)

# plot solution
plot(category_layer(s7), main = "solution", axes = FALSE)
# solve problem
s8 <- solve(p8)# set seed for reproducibility

# plot solution
plot(category_layer(s8), main = "solution", axes = FALSE)
# calculate feature representation in the solution
r8 <- eval_feature_representation_summary(p8, s8)
print(r8)
# plot solution
plot(category_layer(s8), main = "solution", axes = FALSE)


#---------plot selections in planning unit-------------------------



  s8_pu <- NULL
  # Define a function to calculate the most frequent value
  
  calculate_most_frequent <- function(x) {
    
    tbl <-table(x[!is.na(x) & x != 0])
    if (length(tbl) == 0) return(NA)
    most_frequent_value <- as.numeric(names(tbl[which.max(tbl)]))
    return(most_frequent_value)
  }
  
  library(exactextractr)
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- pu_multi
  regions<-st_transform(regions, crs=st_crs(category_layer(s7))) #change s7 to s8 for fixed targets
  pu_cell <- exact_extract(category_layer(s7), regions, progress = FALSE) %>% 
    map_dfr(~ {
      tibble(prio = calculate_most_frequent(.$value))
    }) %>% 
    # join to lookup table to get id
    bind_cols(regions, .)
  # bind to results
  
  s8_pu <- bind_rows(s8_pu, pu_cell)
  
  

  # Convert prio to a factor if needed
  s8_pu$prio <- factor(s8_pu$prio)
  
  
  scenario<-"even_targets"
  
  # Define the default color palette with twelve colors
  default_palette <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
                       "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
  
  abundance_change<-ggplot() +
    # Add base map
    #geom_sf(data = ne_land_cropped, fill = "#dddddd", color = "#888888", lwd = 0.5) +
    
    # Plot trend_pu with fill color based on trend_mean
    geom_sf(data = s8_pu, aes(fill = prio), color = "#888888") +
    
    # Add white layer for zero value in trend_mean
    geom_sf(data = s8_pu[s8_pu$prio == 0, ], fill = "transparent", color = "transpar") +
    # Gradient color scale for trend_mean
    # Create a gradient color scale with the default palette
    # scale_fill_gradientn(colors = default_palette, name = "Months", 
    #                      limits = c(0, 12), na.value = "transparent") +
    scale_fill_manual(values = (setNames(default_palette, levels(s8_pu$prio))), na.value = "transparent") +
    # Gradient color scale for trend_mean
    #scale_fill_gradient(low = "#D7191C", high = "#3C7EA6", name = "Change in Relative Abundance", na.value = "transparent") +
    
    # Set plot title and caption
    labs(caption = paste0("prioritization_", scenario)) +
    
    # Use minimal theme
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "white", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) 
  print(abundance_change)
  ggsave(paste0("prioritization/", scenario, ".pdf"), plot = abundance_change, width = 8, height = 6, units = "in")
  
  # calculate  feature representation statistics based on the prioritization
  tc_s8 <- eval_target_coverage_summary(p8, s8)
  summary(tc_s8$relative_held * 100)
  ## calculate number of features adequately represented by the prioritization
  sum(tc_s8$met)
  ## visualize representation  (values show percent coverage)
  hist(
    tc_s8$relative_held * 100,
    main = "Feature representation by prioritization",
    xlim = c(0, 5),
    xlab = "Percent coverage of features (%)"
  )
  # calculate irreplaceability
  irrep_s8 <- eval_replacement_importance(p8, s8)
  eval_cost_summary(p8, s8)

  
  #-------------explore the representation of current protected areas---------------------
cced<-st_read("data/CCED/CCED_California_Conservation_Easement_Database.shp")
cced <- st_transform(cced, crs = st_crs(pu_multi))
cpad<-st_read("data/cpad_2023a/CPAD_2023a_Holdings.shp")
capd <- st_transform(cpad, crs = crs(pu_multi))
# Define a function to calculate the carbon sequestion potential
add_cced <- function(cell) {
  # Check if points in land_cover_data are within the multipolygon cell
  is_within <- sapply(st_intersects(cced, cell), function(x) length(x) > 0)
  # Extract the protect areas valuefor the corresponding hexagon cell
  gis_areas <- cced[is_within, "gis_acres"]
  gis_areas<-st_drop_geometry(gis_areas)
  # Extract the values from land_cover_pmp_dust for the corresponding hexagon cell
  # Multiply each column by corresponding multiplier
  protect_area <- sum(gis_areas) # 0.58t/acre/year convert square meters to acre and divieded by 12 months

  return(protect_area)
  
  
}  



# Define a function to calculate the carbon sequestion potential
add_cced <- function(cell) {
  # Ensure cell and cpad have the same CRS
  if (st_crs(cell) != st_crs(cced)) {
    cced <- st_transform(cced, crs = st_crs(cell))
  }
  # Check if points in land_cover_data are within the multipolygon cell
  is_within <- sapply(st_intersects(cced, cell), function(x) length(x) > 0)
  # Extract the protect areas valuefor the corresponding hexagon cell
  gap1_areas <- cced[is_within, "GAP1_acres"]
  gap2_areas <- cced[is_within, "GAP2_acres"]
  gis_areas<-st_drop_geometry(gap1_areas)+ st_drop_geometry(gap2_areas)
  # Extract the values from land_cover_pmp_dust for the corresponding hexagon cell
  # Multiply each column by corresponding multiplier
  protect_area <- sum(gis_areas) # 0.58t/acre/year convert square meters to acre and divieded by 12 months
  # Cap values greater than 23715 cell area to 23715 
  protect_area <- pmin(protect_area, 23715)
  return(protect_area)
}  



add_cpad <- function(cell) {
  # Ensure cell and cpad have the same CRS
  if (st_crs(cell) != st_crs(cpad)) {
    cpad <- st_transform(cpad, crs = st_crs(cell))
  }
  # Check if points in land_cover_data are within the multipolygon cell
  is_within <- sapply(st_intersects(cpad, cell), function(x) length(x) > 0)
  # Extract the protect areas valuefor the corresponding hexagon cell
  gap1_areas <- capd[is_within, "GAP1_acres"]
  gap2_areas <- capd[is_within, "GAP2_acres"]
  gis_areas<-st_drop_geometry(gap1_areas)+ st_drop_geometry(gap2_areas)
  # Extract the values from land_cover_pmp_dust for the corresponding hexagon cell
  # Multiply each column by corresponding multiplier
  protect_area <- sum(gis_areas) # 0.58t/acre/year convert square meters to acre and divieded by 12 months
  # Cap values greater than 23715 cell area to 23715 
  protect_area <- pmin(protect_area, 23715)
  return(protect_area)
  
  
}  



# Initialize an empty vector to store land costs
pu_multi$cced <- rep(NA, nrow(pu_multi))

# Iterate over each cell in pu_multi to calculate land cost
for (i in 1:nrow(pu_multi)) {
  cell <- pu_multi[i, "geometry"]
  pu_multi$cced[i] <- add_cced(cell)
}




# Initialize an empty vector to store land costs
pu_multi$cpad <- rep(NA, nrow(pu_multi))

# Iterate over each cell in pu_multi to calculate land cost
for (i in 1:nrow(pu_multi)) {
  cell <- pu_multi[i, "geometry"]
  pu_multi$cpad[i] <- add_cpad(cell)
}


# Create a copy of the original capd values
cpad_transformed <- pu_multi$cpad

# Apply the transformation
cpad_transformed[cpad_transformed > 23715] <- 100
cpad_transformed[cpad_transformed <= 23715] <- (cpad_transformed[cpad_transformed <= 23715] / 23715) * 100

# If you want to replace the original capd values in pu_multi
pu_multi$cpad_percent <- cpad_transformed
pu_multi$cpad_percent<-round(pu_multi$cpad_percent)

# Create a copy of the original capd values
cced_transformed <- pu_multi$cced

# Apply the transformation
cced_transformed[cced_transformed > 23715] <- 100
cced_transformed[cced_transformed <= 23715] <- (cced_transformed[cced_transformed <= 23715] / 23715) * 100

# If you want to replace the original capd values in pu_multi
pu_multi$cced_percent <- cced_transformed
pu_multi$cced_percent<-round(pu_multi$cced_percent)



# Create a copy of the original capd values
pa_transformed <- pu_multi$cced+pu_multi$cpad

# Apply the transformation
pa_transformed[pa_transformed > 23715] <- 100
pa_transformed[pa_transformed <= 23715] <- (pa_transformed[cced_transformed <= 23715] / 23715) * 100

# If you want to replace the original capd values in pu_multi
pu_multi$pa_percent <- pa_transformed
pu_multi$pa_percent<-round(pu_multi$pa_percent)




# Plotting the data
ggplot() +
  # Plot trend_pu with fill color based on pa_percent
  geom_sf(data = pu_multi, aes(fill = cpad), color = "#888888") +
  
  # Add white layer for zero value in capd
  geom_sf(data = pu_multi[pu_multi$cpad <= 1, ], fill = "transparent", color = "transparent") +
  
  # Gradient color scale for capd with custom bins and colors
  scale_fill_gradientn(colors = color_palette, breaks = color_breaks, name = "CAPD Areas",
                       na.value = "transparent") +
  
  # Set plot title and caption
  labs(caption = paste0("CAPD")) +
  
  # Use minimal theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


library(sf)
library(dplyr)

# Calculate quantile-based breaks
quantiles <- quantile(pu_multi$cpad, probs = seq(0, 1, length.out = 9))  # 8 bins
# Ensure breaks are unique
quantiles <- unique(quantiles)

# Ensure that the breaks are unique by adding a small epsilon
epsilon <- 1e-6
for (i in 2:length(quantiles)) {
  if (quantiles[i] == quantiles[i - 1]) {
    quantiles[i] <- quantiles[i] + epsilon
  }
}

# Define a color palette with white for 0-60 and shades of blue for other values
color_palette <- c("white", colorRampPalette(c("#FFFFBF", "#D7191C"))(length(quantiles) - 2))
# Create the plot
ggplot() +
  # Plot trend_pu with fill color based on capd
  geom_sf(data = pu_multi, aes(fill = cut(cpad, breaks = quantiles, include.lowest = TRUE)), color = "#888888") +
  
  # Add white layer for zero value in capd
  geom_sf(data = pu_multi[pu_multi$cpad == 0, ], fill = "transparent", color = "transparent") +
  
  # Gradient color scale for capd with custom bins and colors
  scale_fill_manual(values = setNames(color_palette, levels(cut(pu_multi$cced, breaks = quantiles, include.lowest = TRUE))),
                    name = "CPAD Areas", na.value = "transparent") +
  
  # Set plot title and caption
  labs(caption = "CPAD") +
  
  # Use minimal theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

## Define a color palette with white for zero and a gradient for other values
#color_palette <- c("white", "#D7191C", "#FDAE61", "#FFFFBF", "#ABDDA4", "#2B83BA", "#2C7BB6")
breaks<-c(0, 300, 500, 800, 1200,5000,  10000, 23715)
# Define a color palette with white for 0-60 and shades of blue for other values
color_palette <- c("white", colorRampPalette(c("#FFFFBF", "salmon"))(length(breaks) - 2))
# Create the plot
ggplot() +
  # Plot trend_pu with fill color based on capd
  geom_sf(data = pu_multi, aes(fill = cut(cpad, breaks = breaks, include.lowest = TRUE)), color = "#888888") +
  
  # Add white layer for zero value in capd
  geom_sf(data = pu_multi[pu_multi$cpad == 0, ], fill = "transparent", color = "transparent") +
  
  # Gradient color scale for capd with custom bins and colors
  scale_fill_manual(values = setNames(color_palette, levels(cut(pu_multi$cpad, breaks = breaks, include.lowest = TRUE))),
                    name = "CPAD Areas", labels = c("0-300", "301-500", "501-800", "801-1200", "1201-5000", "5001-10000", "10001-23715"),
                    na.value = "transparent") +
  
  # Set plot title and caption
  labs(caption = "CPAD") +
  
  # Use minimal theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


## Define a color palette with white for zero and a gradient for other values
#color_palette <- c("white", "#D7191C", "#FDAE61", "#FFFFBF", "#ABDDA4", "#2B83BA", "#2C7BB6")
breaks<-c(0, 300, 500, 800, 1200,5000,  10000, 23715)
# Define a color palette with white for 0-60 and shades of blue for other values
color_palette <- c("white", colorRampPalette(c("#FFFFBF", "#D7191C"))(length(breaks) - 2))
# Create the plot
ggplot() +
  # Plot trend_pu with fill color based on capd
  geom_sf(data = pu_multi, aes(fill = cut(cced, breaks = breaks, include.lowest = TRUE)), color = "#888888") +
  
  # Add white layer for zero value in capd
  geom_sf(data = pu_multi[pu_multi$cced == 0, ], fill = "transparent", color = "transparent") +
  
  # Gradient color scale for capd with custom bins and colors
  scale_fill_manual(values = setNames(color_palette, levels(cut(pu_multi$cced, breaks = breaks, include.lowest = TRUE))),
                    name = "CCED Areas", labels = c("0-300", "301-500", "501-800", "801-1200", "1201-5000", "5001-10000", "10001-23715"),
                    na.value = "transparent") +
  
  # Set plot title and caption
  labs(caption = "CCED") +
  
  # Use minimal theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# Ensure breaks are unique
quantiles <- unique(quantiles)

# Calculate quantile-based breaks
quantiles <- quantile(pu_multi$cpad, probs = seq(0, 1, length.out = 9), na.rm = TRUE)  # 8 intervals

# Define a color palette with white for 0 and shades of blue for other values
color_palette <- c("white", colorRampPalette(c("#FFFFBF", "#D7191C"))(length(quantiles) - 2))



# Define breaks and custom color palette
breaks <- c(0, 1, 30, 100)
labels <- c("0", "30", "100")

color_palette <- c(
  "0" = "white",
  "1-30" = colorRampPalette(c("#FFFFBF", "#D7191C"))(30),
  "31-100" = colorRampPalette(c("lightblue",  "#2C7BB6"))(70)
)
# Plotting the data
ggplot() +
  # Plot trend_pu with fill color based on pa_percent
  geom_sf(data = pu_multi, aes(fill = pa_percent), color = "#888888") +
  
  # Add white layer for zero value in capd
  geom_sf(data = pu_multi[pu_multi$pa_percent <= 1, ], fill = "transparent", color = "transparent") +
  
  # Gradient color scale for capd with custom bins and colors
  scale_fill_gradientn(colors = color_palette, breaks = breaks, name = "Protected areas coverage (%)",
                       limit=c(0,100),na.value = "transparent") +
  
  # Set plot title and caption
  labs(caption = paste0("Total Protected areas coverage (%)")) +
  
  # Use minimal theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# Define breaks and custom color palette
breaks <- c(0, 1, 30, 100)
labels <- c("0", "30", "100")
# Plotting the data
ggplot() +
  # Plot trend_pu with fill color based on pa_percent
  geom_sf(data = pu_multi, aes(fill = cced_percent), color = "#888888") +
  
  # Add white layer for zero value in capd
  geom_sf(data = pu_multi[pu_multi$cced_percent <= 1, ], fill = "transparent", color = "transparent") +
  
  # Gradient color scale for capd with custom bins and colors
  scale_fill_gradientn(colors = color_palette, breaks = breaks, name = "CCED coverage (%)",
                       limit=c(0,100),na.value = "transparent") +
  
  # Set plot title and caption
  labs(caption = paste0("CCED coverage (%)")) +
  
  # Use minimal theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
tc_pa <- eval_target_coverage_summary(p2, pu_multi[, "pa_category"])

hist(tc_pa$relative_held * 100,
     main = "Feature representation by existing protected areas",
     col = "skyblue", 
     xlim = c(0, 100),
     ylim=c(0,400),
     xlab = "Percent coverage of features (%)")
# Define breaks and custom color palette
breaks <- c(0, 1, 30, 100)
labels <- c("0", "30", "100")
# Plotting the data
ggplot() +
  # Plot trend_pu with fill color based on pa_percent
  geom_sf(data = pu_multi, aes(fill = cpad_percent), color = "#888888") +
  
  # Add white layer for zero value in capd
  geom_sf(data = pu_multi[pu_multi$cpad_percent <= 1, ], fill = "transparent", color = "transparent") +
  
  # Gradient color scale for capd with custom bins and colors
  scale_fill_gradientn(colors = color_palette, breaks = breaks, name = "CPAD coverage (%)",
                       limit=c(0,100),na.value = "transparent") +
  
  # Set plot title and caption
  labs(caption = paste0("CPAD coverage (%)")) +
  
  # Use minimal theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#create a category with 30% protected as the threshold
pu_multi <- pu_multi %>%
  mutate(pa_category = ifelse(pa_percent < 30, 0, 1))
#create a category with 5% protected as the threshold
pu_multi <- pu_multi %>%
  mutate(pa_category5 = ifelse(pa_percent < 5, 0, 1))

# calculate feature representation statistics based on existing protected areas
tc_pa <- eval_target_coverage_summary(p2, pu_multi[, "pa_category5"])

hist(tc_pa$relative_held * 100,
     main = "Feature representation by existing protected areas",
     col = "skyblue", 
     xlim = c(0, 100),
     ylim=c(0,400),
     xlab = "Percent coverage of features (%)")



# calculate  feature representation statistics based on the prioritization
tc_s8 <- eval_target_coverage_summary(p8, s8)
## visualize representation  (values show percent coverage)
hist(
  tc_s8$relative_held * 100*12,
  col = "#bc80bd", 
  main = "Feature representation by prioritization",
  xlim = c(0, 100),
  xlab = "Percent coverage of features (%)"
)

# calculate  feature representation statistics based on the prioritization
tc_s7$zone <- eval_target_coverage_summary(p7, s7)
## visualize representation  (values show percent coverage)
hist(
  tc_s7$relative_held * 100*12,
  col="#8dd3c7",
  main = "Feature representation by prioritization",
  xlim = c(0, 100),
  xlab = "Percent coverage of features (%)"
)
#print the frame
hist_tc_s8 <- hist(tc_s7$relative_held * 100, plot = FALSE)
ylim_max <- max(hist_tc_s8$counts, na.rm = TRUE)
# Create a blank plot with the desired range of x-axis and y-axis
plot(NULL, xlim =c(0, 100), ylim = c(0, ylim_max*12),
     xlab = "Percent coverage of features (%)", ylab = "Frequency",
     main = "Feature representation comparison")

# Add a legend
legend("topright", legend = c("Existing CCED or CPAD Protected Areas", "Selected dynamic habitat"),
       fill = c("#FAE7AE", "#8dd3c7"))
