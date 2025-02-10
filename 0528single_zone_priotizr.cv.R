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
library(irr)

wd <- ("/Users/lily/Library/CloudStorage/Box-Box/Thesis/SDM/code/")
setwd(wd)
setwd("c:/Users/liliy/Box/Thesis/SDM/code/")

carbpath<-"data/carbSeq_cdl1.tif"
carbSeq_raster<-raster(carbpath)
# Define the new minimum and maximum values for normalization
new_min <- 0  # New minimum value
new_max <- 1  # New maximum value
# Normalize the raster values
normalized_carbSeq <- (carbSeq_raster - minValue(carbSeq_raster)) / (maxValue(carbSeq_raster) - minValue(carbSeq_raster)) * (new_max - new_min) + new_min
plot(normalized_carbSeq, axes = FALSE, ann = FALSE )
# Define the file path where you saved the raster
sagbi_path <- "data/sagbi.tif"
sagbi_raster <- raster(sagbi_path)
# Normalize the raster values
normalized_sagbi <- (sagbi_raster - minValue(sagbi_raster)) / (maxValue(sagbi_raster) - minValue(sagbi_raster)) * (new_max - new_min) + new_min

pu_multi<-st_read("data/pu_multi/pu_multi_cost.shp")

pu_multi<-st_read("data/pu_multi/pu_multi_cropareas_cdl1.shp")
pu_multi$locked_in<-as.logical(pu_multi$lockd_n)
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

#-----------------------------------Create birds abundance feature for full year circle---------------
# Define a function to read raster files from one folder
read_raster_files <- function(folder_path) {
  file_list <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)
  raster_stack <- stack(file_list)
  return(raster_stack)
}

species <- c("ameavo","amgplo", "bknsti", "dunlin", "lobdow", "margod","mouplo", "lobcur", "leasan", "lesyel", "killde", "greyel","pecsan", "uplsan","shbdow", "whimbr")
scenarios<- c("cdl", "dust", "bbau")

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
birds<- stack(raster_stack_list)
birds_spat6<-rast(birds) #go to line 369  to add other features 
#add otheer features
birds7<-stack(birds, normalized_carbSeq)
birds_spat7<-rast(birds7) # go to line 392 to modify cost to line to prioritize

birds8<-stack(birds, normalized_sagbi)
birds_spat8<-rast(birds8) # go to line 392 to modify cost to line to prioritize

birds9<-stack(birds, normalized_carbSeq, normalized_sagbi)
birds_spat9<-rast(birds9) # go to line 392 to modify cost to line to prioritize

birds_spat10<-rast(normalized_sagbi) # go to line 392 to modify cost to line to prioritize
birds_spat11<-rast(normalized_carbSeq)

pu_multi<-st_transform(pu_multi, crs=crs(birds_spat6))
pu_multi$water_costk<-pu_multi$crprs_c *23715*0.3281/2*pu_multi$cost_his/1000
#detemine cost column based on month
pu_multi$costk<-pu_multi$lnd_cdl/2/1000+pu_multi$water_costk*6

pu_multi$costm<-pu_multi$costk/1000

hist(pu_multi$cropars_bbau)
# build minimum set problem
p1 <-
  problem(pu_multi, birds_spat3[[1]], cost_column = "costk") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty =0.03) %>%#0.003
  add_relative_targets(0.25) %>%
  add_locked_in_constraints("locked_in") %>%
  #add_linear_constraints(threshold = 3026, sense = "<=", data = "costk")  %>%
  add_binary_decisions() 
# create vector to store plot names
n <- c()

# create empty list to store solutions
s <- c()

if (require("gurobi")) {
  p1 <- p1 %>% add_gurobi_solver(verbose = FALSE)
  n <- c(n, "gurobi")
  s <- c(s, solve(p1))
}

# solve problem
s1 <- solve(p1)
# Save as a shapefile
#st_write(s1[, "solution_1"], "solution_1_p0.05_t0.25.shp")
# create column for making a map of the prioritization
s1$map_1 <- case_when(
  s1$locked_in > 0.5 ~ "CVPIA",
  s1$solution_1 > 0.5 ~ "Flooded dynamic habitat",
  TRUE ~ "other"
)

# plot map of prioritization
plot(
  s1[, "map_1"], pal = c("skyblue", "salmon", "grey90"),
  main = NULL, key.pos = 1
)


s1_costk<-eval_cost_summary(p1, s1[, "solution_1"])$cost
flooded_areas_p1 <- s1[, "solution_1"] %>%
  st_drop_geometry() %>%
  bind_cols(pu_multi) %>%
  filter(solution_1 == 1) %>%
  summarise(sum_flodd_r = sum(flodd_r)) %>%
  pull(sum_flodd_r)
# Add text with the cost value to the plot
# Get plot extent for appropriate text placement
plot_extent <- par("usr")  # Get the current plot dimensions (x1, x2, y1, y2)

text(x = plot_extent[1] + (plot_extent[2] - plot_extent[1]) * 1,  # 10% from the left
     y = plot_extent[3] + (plot_extent[4] - plot_extent[3]) * 0.8,  # 90% from the bottom
     labels = paste0("Cost =", s1_costk,"\nFlooded areas=", flooded_areas_p1), pos = 2)

# calculate irreplaceability
irrep_s1 <- eval_replacement_importance(p1, s1["solution_1"])
print(irrep_s1)

# manually coerce values for planning units not selected in prioritization
# to NA, so that they are shown in white
irrep_s1$plot_total <- irrep_s1$rc
irrep_s1$plot_total[s1$solution_1 < 0.5] <- NA_real_

# plot map of overall importance scores
plot(st_as_sf(irrep_s1[, "plot_total"]), main = "Overall importance")
s1_costk<-eval_cost_summary(p1, s1[, "solution_1"])$cost

tc_s1 <- eval_target_coverage_summary(p1, s1[,'solution_1'])
## visualize representation  (values show percent coverage)

hist(
  tc_s1$relative_held * 100,
  breaks=(max(tc_s1$relative_held * 100)-min(tc_s1$relative_held * 100))/2,
  col = "salmon", 
  main = "Feature representation by prioritization",
  xlim = c(0, 100),
  ylim=c(0,250),
  xlab = "Percent coverage of features (%)"
)

# Create a blank plot with the desired range of x-axis and y-axis
plot(NULL, xlim =c(0, 100), ylim = c(0, 250),
     xlab = "Percent coverage of features (%)", ylab = "Frequency",
     main = "Feature representation comparison")

# Add a legend
legend("topright", legend = c("CVPIA refugia", "Selected dynamic habitat"),
       fill = c("skyblue", "salmon"))


w1.1 = rep(0, 786)
w1.1[785]<-10000
w1.1[786]<-10000


#build max features problem
p2 <-
  problem(pu_multi, birds_spat3[[1]], cost_column="costm") %>%
  add_max_features_objective(budget = 1475) %>%
  add_relative_targets(0.3) %>%
  add_boundary_penalties(penalty=0.001)%>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions() %>%
  add_gurobi_solver(verbose = TRUE)%>% add_feature_weights(w1.1)



# create vector to store plot names
n <- c()

# create empty list to store solutions
s <- c()

if (require("gurobi")) {
  p2 <- p2 %>% add_gurobi_solver(verbose = FALSE)
  n <- c(n, "gurobi")
  s <- c(s, solve(p2))
}

# solve problem
s2 <- solve(p2)
st_write(s1[, "solution_1"], "s2_1_p0.03_t0.3_1475.shp")
# create column for making a map of the prioritization
s2$map_1 <- case_when(
  s2$locked_in > 0.5 ~ "CVPIA",
  s2$solution_1 > 0.5 ~ "Flooded dynamic habitat",
  TRUE ~ "other"
)

# plot map of prioritization
plot(
  s2[, "map_1"], pal = c("skyblue", "salmon", "grey90"),
  main = NULL, key.pos = 1
)

s2_costk<-eval_cost_summary(p2, s2[, "solution_1"])$cost
flooded_areas_p2 <- s2[, "solution_1"] %>%
  st_drop_geometry() %>%
  bind_cols(pu_multi) %>%
  filter(solution_1 == 1) %>%
  summarise(sum_flodd_r = sum(flodd_r)) %>%
  pull(sum_flodd_r)
plot_extent <- par("usr")  # Get the current plot dimensions (x1, x2, y1, y2)

text(x = plot_extent[1] + (plot_extent[2] - plot_extent[1]) * 1,  # 10% from the left
     y = plot_extent[3] + (plot_extent[4] - plot_extent[3]) * 0.8,  # 90% from the bottom
     labels = paste0("Cost =", s2_costk,"\nFlooded areas=", flooded_areas_p2), pos = 2)
# calculate irreplaceability
irrep_s2 <- eval_replacement_importance(p2, s2["solution_1"])


# manually coerce values for planning units not selected in prioritization
# to NA, so that they are shown in white
irrep_s2$plot_total <- irrep_s2$rc
irrep_s2$plot_total[s2$solution_1 < 0.5] <- NA_real_

# plot map of overall importance scores
plot(st_as_sf(irrep_s2[, "plot_total"]), main = "Overall importance")
s2_costk<-eval_cost_summary(p2, s2[, "solution_1"])$cost

tc_s2 <- eval_target_coverage_summary(p2, s2[,'solution_1'])
## visualize representation  (values show percent coverage)
hist(
  tc_s2$relative_held * 100,
  breaks=(max(tc_s2$relative_held * 100)-min(tc_s2$relative_held * 100))/2,
  col = "salmon", 
  main = "Feature representation by prioritization",
  xlim = c(0, 100),
  ylim=c(0,400),
  xlab = "Percent coverage of features (%)"
)
tc_s3 <- eval_target_coverage_summary(p2, pu_multi[,'pa'])
## visualize representation  (values show percent coverage)
hist(
  tc_s3$relative_held * 100,
  breaks=(max(tc_s3$relative_held * 100)-min(tc_s3$relative_held * 100))/2,
  col = "skyblue", 
  main = "Feature representation by prioritization",
  xlim = c(0, 100),
  #ylim=c(0,250),
  xlab = "Percent coverage of features (%)"
)




#----------blended approach for detemining penalty------------------

# generate boundary length data for the planning units
pu_bd <- boundary_matrix(pu_multi)

# manually re-scale the boundary length values
pu_bd <- rescale_matrix(pu_bd)

# define a range of different penalty values

prelim_penalty<-c(0.001, 0.01, 0.03, 0.1, 5, 10, 15, 20)
#prelim_penalty<-c(0.1, 0.3, 1, 2, 3, 10, 20, 50, 100)
t=0.03
# define a problem without boundary penalties
p0 <-
  problem(pu_multi, birds_spat9, cost_column = "costk") %>%
  add_min_set_objective() %>%
  add_relative_targets(0.3) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions()

# generate preliminary prioritizations based on each penalty
## note that we specify a relaxed gap and time limit for the solver
prelim_blended_results <- lapply(prelim_penalty, function(x) {
  s <-
    p0 %>%
    add_boundary_penalties(penalty = x, data=pu_bd) %>%
    add_default_solver(gap = 0.1, time_limit = 10 * 60) %>%
    solve()
  s <- data.frame(s = s$solution_1)
  names(s) <- with_options(list(scipen = 30), paste0("penalty_", x))
  s
})

# format results as a single spatial object
prelim_blended_results <- cbind(
  pu_multi, do.call(bind_cols, prelim_blended_results)
)



# calculate metrics for blended method prioritizations
## note that we use p0 and not p1 so that cost calculations are based
## on the cost values and not zeros
blended_metrics <- lapply(
  grep("penalty_", names(prelim_blended_results)), function(x) {
    x <- prelim_blended_results[, x]
    data.frame(
      total_cost = eval_cost_summary(p0, x)$cost,
      total_boundary_length = eval_boundary_summary(p0, x)$boundary
    )
  }
)
blended_metrics <- do.call(bind_rows, blended_metrics)
blended_metrics$penalty <- prelim_penalty
blended_metrics <- as_tibble(blended_metrics)

# create data for plotting
result_data1 <-
  blended_metrics %>%
  ## rename threshold column to value column
  rename(value = "penalty") %>%
  ## add column with column names that contain candidate prioritizations
  mutate(name = grep(
    "penalty_", names(prelim_blended_results), value = TRUE, fixed = TRUE
  )) %>%
  ## add column with labels for plotting
  mutate(label = paste("penalty =", value)) %>%
  ## add column to keep track prioritizations selected by different methods
  mutate(method = "none")


# specify prioritization selected by visual method
result_data1$method[3] <- "visual"

# calculate TOPSIS scores
topsis_results <- topsis(
  decision =
    blended_metrics %>%
    dplyr::select(total_cost, total_boundary_length) %>%
    as.matrix(),
  weights = c(1, 1),
  impacts = c("-", "-")
)


# add column indicating prioritization selected by TOPSIS method
result_data1$method[which.max(topsis_results$score)] <- "TOPSIS"

# generate ideal prioritization based on cost criteria

p0 <-
  problem(pu_multi, birds_spat9, cost_column = "costk") %>%
  add_min_set_objective() %>%
  add_relative_targets(0.3) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0)

# solve problem
s0 <- solve(p0)

pu_multi$zeros<-0

# generate ideal prioritization based on spatial fragmentation criteria
## note that any non-zero penalty value would work here,
## so we just use a penalty  of 1
p00 <-
  problem(pu_multi, birds_spat9, cost_column = "zeros") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = t) %>%
  add_relative_targets(0.3) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0)
s00<-solve(p00)

# generate problem formulation with costs and boundary penalties for
# calculating performance metrics
p_metrics <-
  problem(pu_multi, birds_spat9, cost_column = "costk") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = t) %>%
  add_relative_targets(0.3) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions()

# calculate performance metrics for ideal cost prioritization
s0_metrics <- tibble(
  total_cost = eval_cost_summary(p_metrics, s0[, "solution_1"])$cost,
  total_boundary_length =
    eval_boundary_summary(p_metrics, s0[, "solution_1"])$boundary
)

# calculate performance metrics for ideal boundary length prioritization
s00_metrics <- tibble(
  total_cost = eval_cost_summary(p_metrics, s00[, "solution_1"])$cost,
  total_boundary_length =
    eval_boundary_summary(p_metrics, s00[, "solution_1"])$boundary
)

# calculate penalty value based on Cohon et al. 1979
cohon_penalty <- abs(
  (s0_metrics$total_cost - s00_metrics$total_cost) /
    (s0_metrics$total_boundary_length - s00_metrics$total_boundary_length)
)

# round to 5 decimal places to avoid numerical issues during optimization
cohon_penalty <- round(cohon_penalty, 0)

pu_multi$costm<-pu_multi$costk/1000
# generate prioritization using penalty value calculated using Cohon et al. 1979


p6 <-
  problem(pu_multi, birds_spat9, cost_column = "costk") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = cohon_penalty, data=pu_bd) %>%

  add_relative_targets(0.3) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions() %>%add_default_solver(gap = 0.2, time_limit = 10 * 60)

# solve problem
s6 <- solve(p6)

# add new row with data for prioritization generated following Cohon et al. 1979
result_data1 <- bind_rows(
  result_data1,
  tibble(
    total_cost = eval_cost_summary(p6, s6[, "solution_1"])$cost,
    total_boundary_length =
      eval_boundary_summary(p6, s6[, "solution_1"])$boundary,
    value = cohon_penalty,
    name = paste0("penalty_", cohon_penalty),
    label = paste0("Penalty = ",  cohon_penalty),
    method = "Cohon"
  )
)



# create plot to visualize trade-offs and show selected candidate prioritization
result_plot1 <-
  ggplot(
    data =
      result_data1 %>%
      mutate(vjust = if_else(method == "Cohon", -1, 0.5)),
    aes(x = total_boundary_length, y = total_cost, label = label)
  ) +
  geom_line() +
  geom_point(aes(color = method), size = 3) +
  geom_text(aes(vjust = vjust, color = method), hjust = -0.1) +
  scale_color_manual(
    name = "Method",
    values = c(
      "visual" = "#984ea3",
      "none" = "#000000",
      "TOPSIS" = "#e41a1c",
      "Cohon" = "#377eb8"
    )
  ) +
  xlab("Total boundary length of prioritization") +
  ylab("Total cost of prioritization") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.4)))

# Open a PDF graphics device
pdf("penalty_all_month_newcost.pdf")

# Save the ggplot object
print(result_plot1)

# Close the PDF device
dev.off()




#mmm

#-----------------determine the weights--------------------------------
s1_costk<-eval_cost_summary(p1, s1[, "solution_1"])$cost/1000 #convert to million dollars
pu_multi$costm<-pu_multi$costk/1000

weights<-c(0,2000, 6000, 10000, 12000, 20000, 30000, 40000, 50000)
weights1<-c(0,2000, 6000, 10000, 12000, 20000, 30000,40000, 50000)
w1.1 = rep(0, 834)
bird_met<-NULL
water_met<-NULL
bird_rep<-NULL
water_rep<-NULL
carbon_rep<-NULL
carbon_met<-NULL
for(w in 1:9){
  w1.1[833]<-weights[w]
  w1.1[834]<-weights1[w]
  #w1.1[1:730]<-weights[w]
  # create minimal problem that aims to maximize the number of features
  # adequately conserved given a total budget of s1_cost
  p1.1 <-problem(pu_multi, birds_spat9, cost_column="costm") %>%
    add_max_features_objective(budget = s1_costk) %>%
    add_relative_targets(0.3) %>%
    add_boundary_penalties(penalty=0.03)%>%
    add_locked_in_constraints("locked_in") %>%
    add_binary_decisions() %>%
    add_default_solver()%>% add_feature_weights(w1.1)
  
  # create vector to store plot names
  n <- c()
  
  # create empty list to store solutions
  s <- c()
  
  if (require("gurobi")) {
    p1.1 <- p1.1 %>% add_gurobi_solver(verbose = TRUE)
    n <- c(n, "gurobi")
    s <- c(s, solve(p1.1))
  }
  
  s1.1<-solve(p1.1)
  
  repr<-eval_feature_representation_summary(p1.1, s1.1[, "solution_1"])
  
  bird_rep[w]<-mean(repr[1:832, ]$relative_held*100)
  water_rep[w]<-mean(repr[834, ]$relative_held*100)
  carbon_rep[w]<-mean(repr[833, ]$relative_held*100)
  
  cover<-eval_target_coverage_summary(p1.1, s1.1[, "solution_1"])
  
  bird_met[w]<-sum(cover[1:832, ]$met)
  water_met[w]<-sum(cover[834, ]$met)
  carbon_met[w]<-sum(cover[833, ]$met)
}

weight_graph<-data.frame(weights=weights, Bird_abundance=bird_rep, 
                         Carbon_sequestration= carbon_rep, Groundwater_recharge=water_rep  )
# Reshape data for ggplot
weight_graph_long <- tidyr::gather(weight_graph, key = "Feature", value = "Value", -weights)

# Plot
ggplot(weight_graph_long, aes(x = weights, y = Value, color = Feature)) +
  geom_line() +
  labs(x = "Weights", y = "Feature representations", title = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "grey90", color = "black"),  # Grey background
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Black border
  )

s1_costk<-eval_cost_summary(p1, s1[, "solution_1"])$cost/1000 #convert to million dollars
pu_multi$costm<-pu_multi$costk/1000

weights<-c(0,2000, 6000, 10000, 12000, 20000, 30000, 40000, 50000)
weights1<-c(0,2000, 6000, 10000, 12000, 20000, 30000,40000, 50000)
w1.1 = rep(0, 786)
bird_met<-NULL
water_met<-NULL
bird_rep<-NULL
water_rep<-NULL
carbon_rep<-NULL
carbon_met<-NULL
for(w in 1:9){
  w1.1[785]<-weights[w]
  w1.1[786]<-weights1[w]
  #w1.1[1:730]<-weights[w]
  # create minimal problem that aims to maximize the number of features
  # adequately conserved given a total budget of s1_cost
  p1.1 <-problem(pu_multi, birds_spat9, cost_column="costm") %>%
    add_max_features_objective(budget = s1_costk) %>%
    add_relative_targets(0.3) %>%
    add_boundary_penalties(penalty=0.03)%>%
    add_locked_in_constraints("locked_in") %>%
    add_binary_decisions() %>%
    add_default_solver()%>% add_feature_weights(w1.1)
  
  # create vector to store plot names
  n <- c()
  
  # create empty list to store solutions
  s <- c()
  
  if (require("gurobi")) {
    p1.1 <- p1.1 %>% add_gurobi_solver(verbose = TRUE)
    n <- c(n, "gurobi")
    s <- c(s, solve(p1.1))
  }
  
  s1.1<-solve(p1.1)
  
  repr<-eval_feature_representation_summary(p1.1, s1.1[, "solution_1"])
  
  bird_rep[w]<-mean(repr[1:784, ]$relative_held*100)
  water_rep[w]<-mean(repr[786, ]$relative_held*100)
  carbon_rep[w]<-mean(repr[785, ]$relative_held*100)
  
  cover<-eval_target_coverage_summary(p1.1, s1.1[, "solution_1"])
  
  bird_met[w]<-sum(cover[1:784, ]$met)
  water_met[w]<-sum(cover[786, ]$met)
  carbon_met[w]<-sum(cover[785, ]$met)
}

weight_graph<-data.frame(weights=weights, Bird_abundance=bird_rep, 
                         Carbon_sequestration= carbon_rep, Groundwater_recharge=water_rep  )
# Reshape data for ggplot
weight_graph_long <- tidyr::gather(weight_graph, key = "Feature", value = "Value", -weights)

# Plot
ggplot(weight_graph_long, aes(x = weights, y = Value, color = Feature)) +
  geom_line() +
  labs(x = "Weights", y = "Feature representations", title = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "grey90", color = "black"),  # Grey background
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Black border
  )


#-----hierarchical methods to evaluate trade-off between cost, representation and targets-----------------------


# print cost
s1_cost<-eval_cost_summary(p1, s1[, "solution_1"])$cost
sequence <- s1_cost * seq(2e-1, 8e-1, length.out = 4)
# Calculate the sequence of different cost budgets
sequence <- s1_cost * seq(5e-1, 10e-1, length.out = 5)
sequence1 <- s1_cost * seq(1e-1, 3e-1, length.out = 3)

# Combine the values into a list
threshold <- c(s1_cost-rev(sequence), s1_cost, s1_cost + sequence)

threshold<-threshold/1000
# add a column with zeros, here we use zero cost values so that 
#the prioritization will focus exclusively on spatial fragmentation.
pu_multi$zeros <- 0

# define a problem with zero cost values and boundary penalties
## note that because all the costs are all zero, it doesn't actually
## matter what penalty value is used (as long as the value is > 0)
## and so we just use a value of 1
w1.1 = rep(0, 834)
w1.1[833]<-8000
w1.1[834]<-8000

w1.1 = rep(0, 786)
w1.1[785]<-10000
w1.1[786]<-10000

pu_multi$costm<-pu_multi$costk/1000
results_portfolio <- lapply(threshold, function(x) {
  ## generate solution by adding a constraint based on the threshold and
  ## using the "real" cost values (i.e., not zeros)
  p3 <-
    problem(pu_multi, birds_spat9, cost_column="costm") %>%
    add_max_features_objective(budget = x) %>%
    add_relative_targets(0.3) %>%
    add_boundary_penalties(penalty=0.001)%>%
    add_locked_in_constraints("locked_in") %>%
    add_binary_decisions() %>%
    add_gurobi_solver(verbose = TRUE)%>% add_feature_weights(w1.1)
  s<-solve(p3)
  ## return data frame with solution
  s <- data.frame(s = s$solution_1)
  names(s) <- paste0("threshold_", x)
  m<-data.frame(
    feature_representation = mean(eval_feature_representation_summary(p3, s)$relative_held[1:728]),
    conservation_cost = x
  )
  list(solution = s, summary = m)
})
library(dplyr)
#get the solution part of each element of  hirarchical_results
hierarchical_results <- cbind(pu_multi, do.call(bind_cols, lapply(results_portfolio, function(result) result$solution)))
# Assuming hierarchical_results is already loaded as an sf object
numeric_data <- hierarchical_results %>%
  st_drop_geometry() %>%
  dplyr::select(14:22) %>%
  mutate(across(everything(), as.numeric))

numeric_data <- hierarchical_results %>%
  st_drop_geometry() %>%
  dplyr::select(27:50) %>%
  mutate(across(everything(), as.numeric))

# Create a list to store the sums
target_columns <- grep("^threshold_", names(numeric_data), value = TRUE)
flooded_areas <- numeric(length(target_columns))
carbseq <- numeric(length(target_columns))
names(flooded_areas1) <- target_columns
names(carbseq1) <- target_columns
for (target_col in target_columns) {
  flooded_areas[target_col] <- numeric_data %>%
    filter(!!sym(target_col) == 1) %>%
    summarise(sum_flodd_r = sum(flodd_cdl)) %>%
    pull(sum_flodd_r)
  carbseq[target_col] <- numeric_data %>%
    filter(!!sym(target_col) == 1) %>%
    summarise(sum_flodd_r = sum(carbon_cdl)) %>%
    pull(sum_flodd_r)
}
hierarchical_metrics <- do.call(bind_rows, lapply(results_portfolio, function(result) as.data.frame(result$summary)))



# plot maps of prioritizations
plot(
  x =
    hierarchical_results %>%
    dplyr::select(starts_with("threshold_")) %>%
    mutate_if(is.numeric, function(x) {
      case_when(
        hierarchical_results$locked_in > 0.5 ~ "locked in",
        x > 0.5 ~ "priority",
        TRUE ~ "other"
      )
    }),
  pal = c("skyblue", "grey90", "salmon")
)

# calculate metrics for prioritizations


hierarchical_metrics <- as_tibble(hierarchical_metrics)

# preview metrics
print(hierarchical_metrics)
# create data for plotting
result_data <-
  hierarchical_metrics %>%
  ## rename threshold column to value column
  rename(value = "conservation_cost") %>%
  ## add column with column names that contain candidate prioritizations
  mutate(name = grep(
    "threshold_", names(hierarchical_results), value = TRUE, fixed = TRUE
  )) %>%
  ## add column with labels for plotting
  mutate(label = paste("Budget =", value)) %>%
  ## add column to keep track prioritizations selected by different methods
  mutate(method = "none")

# print table
print(result_data)
# specify prioritization selected by visual method
result_data$method[5] <- "cost_effective"

# create plot to visualize trade-offs and show selected candidate prioritization
result_plot <-
  ggplot(
    data = result_data,
    aes(x = value, y = feature_representation, label = label)
  ) +
  geom_line() +
  geom_point(size = 3) +
  geom_text(hjust = -0.15) +
  scale_color_manual(
    
    values = c("cost_effective" = "blue", "none" ="black")
  ) +
  xlab("Cost Budget") +
  ylab("Feature representation") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.4))) +
  theme(legend.title = element_blank())
# render plot
print(result_plot)

targets=seq(0.05, 0.5, length=10)
results_portfolio1 <- lapply(targets, function(x) {
  ## generate solution by adding a constraint based on the threshold and
  ## using the "real" cost values (i.e., not zeros)
  p3 <-
    problem(pu_multi, birds_spat9, cost_column="costm") %>%
    add_min_set_objective() %>%
    add_relative_targets(x) %>%
    add_boundary_penalties(penalty=0.001)%>%
    add_locked_in_constraints("locked_in") %>%
    add_binary_decisions() %>%
    add_gurobi_solver(verbose = TRUE)
  s<-solve(p3)
  ## return data frame with solution
  s <- data.frame(s = s$solution_1)
  names(s) <- paste0("target_", x)
  c<-data.frame(
    conservation_cost = eval_cost_summary(p3, s)$cost,
    conservation_target = x
  )
  list(solution = s, summary = c)
})


#get the solution part of each element of  hirarchical_results
hierarchical_results1 <- cbind(pu_multi, do.call(bind_cols, lapply(results_portfolio1, function(result) result$solution)))

hierarchical_metrics1 <- do.call(bind_rows, lapply(results_portfolio1, function(result) as.data.frame(result$summary)))

# Assuming hierarchical_results is already loaded as an sf object
numeric_data1 <- hierarchical_results1 %>%
  st_drop_geometry() %>%
  dplyr::select(27:49) %>%
  mutate(across(everything(), as.numeric))

# Create a list to store the sums
target_columns <- grep("^target_", names(numeric_data1), value = TRUE)
flooded_areas1 <- numeric(length(target_columns))
names(flooded_areas) <- target_columns
carbseq1 <- numeric(length(target_columns))
names(carbseq ) <- target_columns
for (target_col in target_columns) {
  flooded_areas1[target_col] <- numeric_data1 %>%
    filter(!!sym(target_col) == 1) %>%
    summarise(sum_flodd_r = sum(flodd_r)) %>%
    pull(sum_flodd_r)
  carbseq1[target_col] <- numeric_data1 %>%
    filter(!!sym(target_col) == 1) %>%
    summarise(sum_flodd_r = sum(carbon)) %>%
    pull(sum_flodd_r)
}


# plot maps of prioritizations
plot(
  x =
    hierarchical_results1 %>%
    dplyr::select(starts_with("target_")) %>%
    mutate_if(is.numeric, function(x) {
      case_when(
        hierarchical_results$locked_in > 0.5 ~ "locked in",
        x > 0.5 ~ "priority",
        TRUE ~ "other"
      )
    }),
  pal = c("skyblue", "grey90", "salmon")
)

# calculate metrics for prioritizations


hierarchical_metrics1 <- as_tibble(hierarchical_metrics1)


# create data for plotting
result_data1 <-
  hierarchical_metrics1 %>%
  ## rename threshold column to value column
  rename(value = "conservation_target") %>%
  ## add column with column names that contain candidate prioritizations
  mutate(name = grep(
    "target_", names(hierarchical_results1), value = TRUE, fixed = TRUE
  )) %>%
  ## add column with labels for plotting
  mutate(label = paste("target=", value)) %>%
  ## add column to keep track prioritizations selected by different methods
  mutate(method = "none")

# print table
print(result_data1)
# specify prioritization selected by visual method
result_data$method[5] <- "cost_effective"

# create plot to visualize trade-offs and show selected candidate prioritization
result_plot1 <-
  ggplot(
    data = result_data1,
    aes(x = conservation_cost, y = value, label = label)
  ) +
  geom_line() +
  geom_point(size = 3) +
  geom_text(hjust = -0.15) +
  scale_color_manual(
    values = c("cost_effective" = "blue", "none" ="black")
  ) +
  xlab("Cost") +
  ylab("Target") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.4))) +
  theme(legend.title = element_blank())
# render plot
print(result_plot1)








#--------climate scenarios cluster analysis------------------------------------------------------
birds_spat<-vector(mode = "list", length = 3)
birds_spat1<-vector(mode = "list", length = 3)
birds_spat2<-vector(mode = "list", length = 3)
birds_spat3<-vector(mode = "list", length = 3)
birds_spat4<-vector(mode = "list", length = 3)
birds_spat5<-vector(mode = "list", length = 3)
carbon<-vector(mode = "list", length = 3)


costk_scenarios<-c(pu_multi$lnd_cdl/2/1000+pu_multi$crprs_c*23715*0.3281/2*pu_multi$cost_his/1000*6,
                   pu_multi$lnd_cst/2/1000+pu_multi$cropars*23715*0.3281/2*pu_multi$cost_dust/1000*6,
                   pu_multi$land_bb/2/1000+pu_multi$crprs_b*23715*0.3281/2*pu_multi$cost_bbau/1000*6)
path_carbon<-c("data/carbSeq_cdl1.tif", "data/carbSeq_dust1.tif","data/carbSeq_bbau1.tif")

prt_costk=NULL



for (i in 1:3){
  carbpath<-path_carbon[i]
  carbSeq_raster<-raster(carbpath)
  # Define the new minimum and maximum values for normalization
  new_min <- 0  # New minimum value
  new_max <- 1  # New maximum value
  # Normalize the raster values
  carbon[[i]] <- (carbSeq_raster - minValue(carbSeq_raster)) / (maxValue(carbSeq_raster) - minValue(carbSeq_raster)) * (new_max - new_min) + new_min
  
  
}

# Define a function to read raster files from one folder
read_raster_files <- function(folder_path) {
  file_list <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)
  raster_stack <- stack(file_list)
  return(raster_stack)
}

scenarios<-c("cdl", "dust", "bbau")

# Populate folder_paths with folder names
for (j in 1:length(scenarios)) {
  # Define a vector to store folder paths
  folder_paths <- character(length(species))
  # Create an empty list to store raster bricks
  raster_stack_list <- list()
  for (i in 1:length(species)) {
    folder_paths[i] <- paste0(species[i], "_abundance_10year_", scenarios[j])
  }
  
  
  # Iterate over each folder path
  for (folder_path in folder_paths) {
    raster_stack <- read_raster_files(folder_path)
    raster_stack_list[[folder_path]] <- raster_stack
  }
  
  # Combine all raster bricks into a single SpatRaster
  birds <- stack(raster_stack_list)
  birds_spat[[j]]<-rast(birds) #go to line 369  to add other features 
  #add otheer features
  birds1<-stack(birds, carbon[[j]])
  birds_spat1[[j]]<-rast(birds1) # go to line 392 to modify cost to line to prioritize
  
  birds2<-stack(birds, normalized_sagbi)
  birds_spat2[[j]]<-rast(birds2) # go to line 392 to modify cost to line to prioritize
  
  birds3<-stack(birds, carbon[[j]], normalized_sagbi)
  birds_spat3[[j]]<-rast(birds3) # go to line 392 to modify cost to line to prioritize
  
  birds_spat4[[j]]<-rast(normalized_sagbi) # go to line 392 to modify cost to line to prioritize
  birds_spat5[[j]]<-rast(carbon[[j]])
  
}

pu_multi$locked_in<-as.logical(pu_multi$lockd_n)

prt_results<-vector(mode = "list", length = 3)
prt_costk<-NULL
# create new problem with a portfolio added to it
for (k in 1:length(scenarios)){
  #detemine cost column based on month
  pu_multi$costk<-costk_scenarios[k]
  pt_scenario <-
    problem(pu_multi, birds_spat3[[k]], cost_column = "costk") %>%
    add_min_set_objective() %>%
    add_boundary_penalties(penalty = 0.03) %>%
    add_relative_targets(0.3) %>%
    add_locked_in_constraints("locked_in") %>%
    add_binary_decisions()  %>%
    add_gap_portfolio(number_solutions = 30, pool_gap = 0.2)
  
  # generate prioritizations
  prt <- solve(pt_scenario)
  
  # extract solutions
  prt_result <- sf::st_drop_geometry(prt)
  prt_results[[k]] <- prt_result[, startsWith(names(prt_result), "solution_")]
  prt_costk[k]<-eval_cost_summary(pt_scenario, prt[, "solution_1"])$cost
}

# Rename columns for the first tibble
combined_tibble <- prt_results[[1]] %>%
  rename_with(~ paste0("his_", seq_along(.)), starts_with("solution"))

# Rename columns for the second tibble
combined_tibble <- bind_cols(
  combined_tibble,
  prt_results[[2]] %>%
    rename_with(~ paste0("hot_", seq_along(.)), starts_with("solution"))
)

# Rename columns for the third tibble
combined_tibble <- bind_cols(
  combined_tibble,
  prt_results[[3]] %>%
    rename_with(~ paste0("warm_", seq_along(.)), starts_with("solution"))
)


# Extract the solutions from prt_results as a matrix
solution_matrix <- as.matrix(combined_tibble)
# Initialize an empty matrix to store pairwise kappa coefficients
pairwise_kappa_matrix <- matrix(NA, ncol = ncol(solution_matrix), nrow = ncol(solution_matrix))
pairwise_kappa_matrix <- matrix(NA, ncol(solution_matrix), ncol(solution_matrix),
                                dimnames = list(colnames(solution_matrix), colnames(solution_matrix)))

# Loop through each pair of columns and calculate kappa coefficient
for (i in 1:(ncol(solution_matrix) - 1)) {
  for (j in (i + 1):ncol(solution_matrix)) {
    # Calculate kappa coefficient for the pair of columns
    kappa <- kappam.fleiss(cbind(solution_matrix[, i], solution_matrix[, j]))
    
    # Store the kappa coefficient in the pairwise kappa matrix
    pairwise_kappa_matrix[i, j] <- kappa$value
    pairwise_kappa_matrix[j, i] <- kappa$value  # Since it's symmetric, store in both positions
  }
}

# Calculate Fleiss' Kappa coefficient matrix
kappa_matrix <- kappam.fleiss(solution_matrix)


# Define color scale and corresponding values
color_scale <- colorRampPalette(c("red", "white", "blue"))(100)
scale_values <- seq(0, 1, length.out = 5)


# Plot the matrix as a heatmap
heatmap(pairwise_kappa_matrix,
        col = colorRampPalette(c("red", "white", "blue"))(100),  # Define red-white-blue color palette
        scale = "none",                                         # Do not scale values
        symm = TRUE,                                            # Symmetric color scale
        main = "Pairwise Kappa Coefficients",                   # Main title
        xlab = "Solutions", ylab = "Solutions")                 # Axis labels

# Add legend

legend("topright", legend = c("Low", "Medium", "High"), fill = colorRampPalette(c("red", "white", "blue"))(3))


# calculate pair-wise distances between different prioritizations for analysis
prt_dists <- vegan::vegdist(t(combined_tibble), method = "jaccard", binary = TRUE)

# run cluster analysis
prt_clust <- hclust(as.dist(prt_dists), method = "average")



library(factoextra)
library(plotly)
# Plot clusters in a "cloud" figure
fviz_dend(prt_clust, k = 6, rect = TRUE, cex = 0.6, horiz=FALSE)

# run k-medoids analysis
prt_med <- pam(prt_dists, k = 6)
# Extract cluster assignments
cluster_assignments <- prt_med$clustering

# Plot solutions in a 2D space
plot(prt_dists, col = cluster_assignments, main = "2D Plot of Clusters")

# extract names of prioritizations that are most central for each group.
prt_med_names <- prt_med$medoids
print(prt_med_names)

# create a copy of prt and set values for locked in planning units to -1
# so we can easily visualize differences between prioritizations
prt2 <- combined_tibble[, prt_med_names]
prt2[which(pu_multi$locked_in > 0.5), prt_med_names] <- -1
prt2<-bind_cols(prt2,pu_multi$geometry)
# plot a map showing main different prioritizations
# dark grey: locked in planning units
# grey: planning units not selected
# green: selected planning units
plot(st_as_sf(prt2), pal = c("skyblue", "grey90", "salmon"))


#-----------With varying targets------------------------
# generate boundary length data for the planning units
pu_bd <- boundary_matrix(pu_multi)

# manually re-scale the boundary length values
pu_bd <- rescale_matrix(pu_bd)

# Define a function to read raster files from one folder
patterns<-c("^week([0-5])\\.tif$", "^week([6-9])\\.tif$","^week(1[0-4])\\.tif$","^week(1[5-8])\\.tif$",
            "^week(1[9]|2[0-2])\\.tif$", "^week(2[3-6])\\.tif$","^week(2[7-9]|3[0-1])\\.tif$", "^week(3[2-6])\\.tif$",
             "^week(3[7-9]|4[0])\\.tif$", "^week(4[1-4])\\.tif$", "^week(4[5-8])\\.tif$" ,"^week(4[9]|5[0-2])\\.tif$")

d<-c(50,40,50,40,40,40,50,50,40,40,40,40)
tar<-c(0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.4)
#penalties <- c(0.00017,0.1, 0.1, 0.00035, 10, 10, 10, 20, 20, 10,0.0034, 0.0001)
#penalties <- c(0.00017,0.0003, 0.0003, 0.00035, 1, 1, 1, 1, 1, 1,0.00034, 0.0001)
cost<-list(pu_multi$cost/2000, pu_multi$cost/1000, pu_multi$cost/1000+pu_multi$land_cost/6000)
index_cost<-c(1,2,2,2,3,3,3,3,3,3,2,1)
conservation_cost<-NULL
result_plot_group<-NULL

#mmm
for(k in 1:12) {
  # Create the list of targets
  t <- c(rep(tar[k], d[k]))
  
  #detemine cost column based on month
  pu_multi$costk<-unlist(cost[index_cost[k]])
  
  
  read_raster_files <- function(folder_path) {
    file_list <- list.files(folder_path, pattern = patterns[k], full.names = TRUE) #pattern reference "^week(3[8-9]|4[0-1])\\.tif$"
    raster_stack <- stack(file_list)
    return(raster_stack)
  }
  
  
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
  #add otheer features
  birds3<-stack(birds, normalized_carbSeq, normalized_sagbi)
  birds_spat3<-rast(birds3) # go to line 392 to modify cost to line to prioritize
  
  
  # create problem
  
  # build problem
  p3 <-
    problem(pu_multi, birds_spat,cost_column = "costk") %>%
    add_min_set_objective() %>%
    add_boundary_penalties(penalty = penalties[k]) %>%
    add_relative_targets(t)  %>%
    add_locked_in_constraints("locked_in") %>%
    add_binary_decisions() %>%
    add_default_solver(gap = 0)
  
  
  
  # create vector to store plot names
  n <- c()
  
  # create empty list to store solutions
  s <- c()
  
  if (require("gurobi")) {
    p4 <- p3 %>% add_gurobi_solver(verbose = FALSE)
    n <- c(n, "gurobi")
    s <- c(s, solve(p2))
  }
  
  # solve problem
  s2 <- solve(p4)
  
  
  # create column for making a map of the prioritization
  s2$map_1 <- case_when(
    s2$locked_in > 0.5 ~ "CVPIA",
    s2$solution_1 > 0.5 ~ "Flooded dynamic Habitat",
    TRUE ~ "other"
  )
  
  # Start PDF device
  pdf(paste0("week",k,"prioritization_map1.pdf"))
  
  # plot map of prioritization
  plot(
    s2[, "map_1"], pal = c("skyblue", "salmon", "grey90"),
    main = NULL, key.pos = 1
  )
  dev.off()
  # calculate cost
  conservation_cost[k] <- eval_cost_summary(p3, s2[, "solution_1"])$cost
  
  # calculate irreplaceability
  irrep_s2 <- eval_ferrier_importance(p3, s2["solution_1"])
  
  # manually coerce values for planning units not selected in prioritization
  # to NA, so that they are shown in white
  irrep_s1$plot_total <- irrep_s1$total
  irrep_s1$plot_total[s1$solution_1 < 0.5] <- NA_real_
  
  
  # Open a PDF graphics device
  pdf(paste0("month", k, "_overall_importance_map1.pdf"))
  
  # Plot the map of overall importance scores
  plot(st_as_sf(irrep_s1[, "plot_total"]), main = "Overall importance")
  
  # Close the PDF device
  dev.off()
  
  # Convert to a data frame
  plot_total_df <- data.frame(plot_total = irrep_s1$plot_total)
  
  # Save to CSV
  write.csv(plot_total_df, paste0("month", k, "_overall_importance.csv"), row.names = FALSE)
  
}
#----------trade-off analysis----------




# print cost
print(s2_cost)

# calculate cost threshold values
threshold <- s2_cost + (s2_cost * seq(1e-5, 2, length.out = 12))
threshold <- ceiling(threshold)

# print cost thresholds
print(threshold)

# add a column with zeros, here we use zero cost values so that 
#the prioritization will focus exclusively on spatial fragmentation.
pu_multi$zeros <- 0

# define a problem with zero cost values and boundary penalties
## note that because all the costs are all zero, it doesn't actually
## matter what penalty value is used (as long as the value is > 0)
## and so we just use a value of 1
p2 <-
  problem(pu_multi, birds_spat3,cost_column = "zeros") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = 0.05) %>%
  add_relative_targets(t)  %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions() 
  


# generate prioritizations based on each cost threshold
## note that the prioritizations are solved to within 10% of optimality
## (the default gap) because the gap is not specified
hierarchical_results <- lapply(threshold, function(x) {
  ## generate solution by adding a constraint based on the threshold and
  ## using the "real" cost values (i.e., not zeros)
  s <-
    p2 %>%
    #add linear constraints to ensure that the total cost of the prioritization 
    #does not exceed a given cost threshold 
    add_linear_constraints(threshold = x, sense = "<=", data = "costk") %>%
    solve()
  ## return data frame with solution
  s <- data.frame(s = s$solution_1)
  names(s) <- paste0("threshold_", x)
  s
})

# format results as a single spatial object
hierarchical_results <- cbind(pu_multi, do.call(bind_cols, hierarchical_results))

# plot maps of prioritizations
plot(
  x =
    hierarchical_results %>%
    dplyr::select(starts_with("threshold_")) %>%
    mutate_if(is.numeric, function(x) {
      case_when(
        hierarchical_results$locked_in > 0.5 ~ "locked in",
        x > 0.5 ~ "priority",
        TRUE ~ "other"
      )
    }),
  pal = c("skyblue", "grey90", "salmon")
)

# calculate metrics for prioritizations
## note that we use p0 and not p1 so that cost calculations are based
## on the cost values and not zeros
hierarchical_metrics <- lapply(
  grep("threshold_", names(hierarchical_results)), function(x) {
    x <- hierarchical_results[, x]
    data.frame(
      total_cost = eval_cost_summary(p4, x)$cost,
      total_boundary_length = eval_boundary_summary(p4, x)$boundary
    )
  }
)
hierarchical_metrics <- do.call(bind_rows, hierarchical_metrics)
hierarchical_metrics$threshold <- threshold
hierarchical_metrics <- as_tibble(hierarchical_metrics)

# preview metrics
print(hierarchical_metrics)
# create data for plotting
result_data <-
  hierarchical_metrics %>%
  ## rename threshold column to value column
  rename(value = "threshold") %>%
  ## add column with column names that contain candidate prioritizations
  mutate(name = grep(
    "threshold_", names(hierarchical_results), value = TRUE, fixed = TRUE
  )) %>%
  ## add column with labels for plotting
  mutate(label = paste("Threshold =", value)) %>%
  ## add column to keep track prioritizations selected by different methods
  mutate(method = "none")

# print table
print(result_data)

# create plot to visualize trade-offs and show selected candidate prioritization
result_plot <-
  ggplot(
    data = result_data,
    aes(x = total_boundary_length, y = total_cost, label = label)
  ) +
  geom_line() +
  geom_point(size = 3) +
  geom_text(hjust = -0.15) +
  scale_color_manual(
    values = c("visual" = "blue", "not selected" ="black")
  ) +
  xlab("Total boundary length of prioritization") +
  ylab("Total cost of prioritization") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.4))) +
  theme(legend.title = element_blank())
# render plot
print(result_plot)
# specify prioritization selected by visual method
result_data$method[3] <- "visual"



library(topsis)
# calculate TOPSIS scores
topsis_results <- topsis(
  decision =
    hierarchical_metrics %>%
    dplyr::select(total_cost, total_boundary_length) %>%
    as.matrix(),
  weights = c(1, 1),
  impacts = c("-", "-")
)

# print results
print(topsis_results)
# add column indicating prioritization selected by TOPSIS method
result_data$method[which.max(topsis_results$score)] <- "TOPSIS"

#-----COHON method

# generate ideal prioritization based on cost criteria

p0 <-
  problem(pu_multi, birds_spat3, cost_column = "costk") %>%
  add_min_set_objective() %>%
  add_relative_targets(t) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0)

# solve problem
s0 <- solve(p0)

# generate ideal prioritization based on spatial fragmentation criteria
## note that any non-zero penalty value would work here,
## so we just use a penalty  of 1
p00 <-
  problem(pu_multi, birds_spat3, cost_column = "zeros") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = 0.05) %>%
  add_relative_targets(t) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0)
s00<-solve(p00)

# generate problem formulation with costs and boundary penalties for
# calculating performance metrics
p_metrics <-
  problem(pu_multi, birds_spat3, cost_column = "costk") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = 0.05) %>%
  add_relative_targets(t) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions()

# calculate performance metrics for ideal cost prioritization
s0_metrics <- tibble(
  total_cost = eval_cost_summary(p_metrics, s0[, "solution_1"])$cost,
  total_boundary_length =
    eval_boundary_summary(p_metrics, s0[, "solution_1"])$boundary
)

# calculate performance metrics for ideal boundary length prioritization
s00_metrics <- tibble(
  total_cost = eval_cost_summary(p_metrics, s00[, "solution_1"])$cost,
  total_boundary_length =
    eval_boundary_summary(p_metrics, s00[, "solution_1"])$boundary
)

# calculate penalty value based on Cohon et al. 1979
cohon_penalty <- abs(
  (s0_metrics$total_cost - s00_metrics$total_cost) /
    (s0_metrics$total_boundary_length - s00_metrics$total_boundary_length)
)

# round to 5 decimal places to avoid numerical issues during optimization
cohon_penalty <- round(cohon_penalty, 5)

# print penalty value
print(cohon_penalty)

# generate prioritization using penalty value calculated using Cohon et al. 1979
p6 <-
  problem(pu_multi, birds_spat3, cost_column = "costk") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = cohon_penalty) %>%
  add_relative_targets(t) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions()

# solve problem
s6 <- solve(p6)

# add new row with data for prioritization generated following Cohon et al. 1979
result_data <- bind_rows(
  result_data,
  tibble(
    total_cost = eval_cost_summary(p6, s6[, "solution_1"])$cost,
    total_boundary_length =
      eval_boundary_summary(p6, s6[, "solution_1"])$boundary,
    value = cohon_penalty,
    name = paste0("penalty_", cohon_penalty),
    label = paste0("Penalty = ",  cohon_penalty),
    method = "Cohon"
  )
)



# create plot to visualize trade-offs and show selected candidate prioritization
result_plot <-
  ggplot(
    data =
      result_data %>%
      mutate(vjust = if_else(method == "Cohon", -1, 0.5)),
    aes(x = total_boundary_length, y = total_cost, label = label)
  ) +
  geom_line() +
  geom_point(aes(color = method), size = 3) +
  geom_text(aes(vjust = vjust, color = method), hjust = -0.1) +
  scale_color_manual(
    name = "Method",
    values = c(
      "visual" = "#984ea3",
      "none" = "#000000",
      "TOPSIS" = "#e41a1c",
      "Cohon" = "#377eb8"
    )
  ) +
  xlab("Total boundary length of prioritization") +
  ylab("Total cost of prioritization") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.4)))

# render plot
print(result_plot)

# extract column names for creating the prioritizations
visual_name <- result_data$name[[which(result_data$method == "visual")]]
topsis_name <- result_data$name[[which(result_data$method == "TOPSIS")]]

# create object with selected prioritizations
solutions  <- bind_cols(
  pu_multi,
  hierarchical_results %>%
    st_drop_geometry() %>%
    dplyr::select(all_of(c(visual_name, topsis_name))) %>%
    setNames(c("Visual", "TOPSIS")),
  s6 %>%
    st_drop_geometry() %>%
    dplyr::select(solution_1) %>%
    rename(Cohon = "solution_1")
)

# plot maps of selected prioritizations
plot(
  x =
    solutions %>%
    dplyr::select(Visual, TOPSIS, Cohon) %>%
    mutate_if(is.numeric, function(x) {
      case_when(
        hierarchical_results$locked_in > 0.5 ~ "locked in",
        x > 0.5 ~ "priority",
        TRUE ~ "other"
      )
    }),
  pal = c("skyblue", "grey90", "salmon")
)

#-------------map different objectives-------------------








  p1 <-
    problem(pu_multi, birds_spat1[[1]], cost_column = "costk") %>%
    add_min_set_objective() %>%
    add_boundary_penalties(penalty = 0.03) %>%
    add_relative_targets(0.3) %>%
    add_locked_in_constraints("locked_in") %>%
    add_binary_decisions() 
  
  print(presolve_check(p1))
  
  # create vector to store plot names
  n <- c()
  
  # create empty list to store solutions
  s <- c()
  
  if (require("gurobi")) {
    p1 <- p1 %>% add_gurobi_solver(verbose = FALSE)
    n <- c(n, "gurobi")
    s <- c(s, solve(p1))
  }
  
  # solve problem
  s1 <- solve(p1)
  
  # create column for making a map of the prioritization
  s1$map_1 <- case_when(
    s1$locked_in > 0.5 ~ "CVPIA",
    s1$solution_1 > 0.5 ~ "Flooded dynamic habitat",
    TRUE ~ "other"
  )
  
  conservation_cost<-round(eval_cost_summary(p1, s1[, "solution_1"])$cost)
  # Open a PDF graphics device
  pdf(paste0("All_months_prio_cdl_birdcarbon.pdf"))
  
  plot(
    s1[, "map_1"], pal = c("skyblue", "salmon", "grey90"),
    main = NULL, key.pos = 1
  )
  plot_extent <- par("usr")  # Get the current plot dimensions (x1, x2, y1, y2)
  
  text(x = plot_extent[1] + (plot_extent[2] - plot_extent[1]) * 1,  # 10% from the left
       y = plot_extent[3] + (plot_extent[4] - plot_extent[3]) * 0.8,  # 90% from the bottom
       labels = paste0("Cost =", conservation_cost), pos = 2)
  
 
  
  # Close the PDF device
  dev.off()
  
  # calculate irreplaceability
  irrep_s1 <- eval_replacement_importance(p1, s1["solution_1"])
  
  
  # manually coerce values for planning units not selected in prioritization
  # to NA, so that they are shown in white
  irrep_s1$plot_total <- irrep_s1$rc
  irrep_s1$plot_total[s1$solution_1 < 0.5] <- NA_real_
  
  # plot map of overall importance scores
  
  
  # Open a PDF graphics device
  pdf(paste0("All_months_irrep_cdl_birdcarbon.pdf"))
  
  # Save the ggplot object
  plot(st_as_sf(irrep_s1[, "plot_total"]), main = "Overall importance")
  
  # Close the PDF device
  dev.off()
  


#--------------------------
  pu_multi$pa <- round(pu_multi$locked_in)
  
  # calculate feature representation statistics based on existing protected areas
  tc_pa <- eval_target_coverage_summary(p1, pu_multi[, "pa"])
  
  print(tc_pa)
  
  # calculate  feature representation statistics based on the prioritization
  tc_s1 <- eval_target_coverage_summary(p1, s1[, "solution_1"])
  print(tc_s1)
  
  # explore representation by existing protected areas
  ## calculate number of features adequately represented by existing protected
  ## areas
  sum(tc_pa$met)
  
  ## summarize representation (values show percent coverage)
  summary(tc_pa$relative_held * 100)
  
  ## visualize representation  (values show percent coverage)
  hist(tc_pa$relative_held * 100,
       main = "Feature representation by existing protected areas",
       col = "skyblue", 
       xlim = c(0, 100),
       ylim=c(0,400),
       xlab = "Percent coverage of features (%)")
  
  
  # explore representation by prioritization
  ## summarize representation (values show percent coverage)
  summary(tc_s1$relative_held * 100)
  
  ## calculate number of features adequately represented by the prioritization
  sum(tc_s1$met)
  
  ## visualize representation  (values show percent coverage)
  hist(
    tc_s1$relative_held * 100,
    col = "salmon", 
    main = "Feature representation by prioritization",
    xlim = c(0, 100),
    xlab = "Percent coverage of features (%)"
  )
  
  
  # Calculate the histograms
  hist_tc_pa <- hist(tc_pa$relative_held * 100, plot = FALSE)
  hist_tc_s1 <- hist(tc_s1$relative_held * 100, plot = FALSE)
  
  # Determine the y-axis limit as the maximum of the maximum counts
  ylim_max <- max(max(hist_tc_pa$counts, na.rm = TRUE), max(hist_tc_s1$counts, na.rm = TRUE))
  
  # Create a blank plot with the desired range of x-axis and y-axis
  plot(NULL, xlim =c(0, 100), ylim = c(0, ylim_max),
       xlab = "Percent coverage of features (%)", ylab = "Frequency",
       main = "Feature representation comparison")
  
  # Plot the histogram for existing protected areas with one color
  barplot(hist_tc_pa$counts, col = "skyblue", add = TRUE)
  
  # Determine the width of each bar
  bar_width <- diff(hist_tc_s1$breaks)[1]
  
  # Calculate the x-axis positions for the bars in the second histogram
  x_pos_s1 <- hist_tc_s1$mids + bar_width/2
  
  # Plot the histogram for prioritization with another color
  barplot(hist_tc_s1$counts, col = "salmon", add = TRUE)
  
  # Create a blank plot with the desired range of x-axis and y-axis
  plot(NULL, xlim =c(0, 100), ylim = c(0, ylim_max),
       xlab = "Percent coverage of features (%)", ylab = "Frequency",
       main = "Feature representation comparison")
  
  # Add a legend
  legend("topright", legend = c("Existing CCED or CPAD Protected Areas", "Prioritization"),
         fill = c("#FAE7AE", "#bc80bd"))
  
#map overlap with cced and cpad
results_year<- st_drop_geometry(s2[, "solution_1"])
results_tibble<-bind_cols(results_year, pu_multi$pa_category) 
names(results_tibble)<-c("Selected habitats", "Protected areas")


  # Extract the solutions from prt_results as a matrix
  results_matrix <- as.matrix(results_tibble)
  # Initialize an empty matrix to store pairwise kappa coefficients
  
  # Calculate Fleiss' Kappa coefficient matrix
  kappa_matrix <- kappam.fleiss(results_matrix)
  
  # Define the Kappa coefficient value
  kappa_value <- 0.12

pu_multi$selected_habitat<-results_year

#create a category with 30% protected as the threshold
pu_multi <- pu_multi %>%
  mutate(pa_category = ifelse(pa_percent < 30, 0, 1))
#create a category with 5% protected as the threshold
pu_multi <- pu_multi %>%
  mutate(pa_category5 = ifelse(pa_percent < 5, 0, 1))

# Create a new column 'combined'
pu_multi <- pu_multi %>%
  mutate(combined = case_when(
    `selected_habitat` == 1 & `p_ctgry` == 1 ~ "overlap",
    `selected_habitat` == 1 ~ "selected_habitat",
    `p_ctgry` == 1 ~ "pa_category",
    `selected_habitat` == 0 & `p_ctgry` == 0 ~ "not selected"
  ))


  # Define the color palette with distinct colors
  # Define the color palette with distinct colors
  color_palette <- c(
    "selected_habitat" = "salmon",
    "pa_category" = "skyblue",
    "overlap" = "#FAE7AE",
    "not selected" = "grey90"
  )
  
  # Assuming kappa_value is already defined
  kappa_value <- 0.12
  
  ggplot(data = pu_multi) +
    geom_sf(aes(fill = combined)) +
    scale_fill_manual(
      values = color_palette,
      name = "Category"
    ) +
    # labs(
    #   title = "Selected Habitats and Protected Areas",
    #   caption = "Blue: Selected Habitats, Green: Protected Areas, Purple: Overlap, Grey: None"
    # ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = -2200000, y = 2240000, label = paste("Kappa coefficient =", kappa_value), size = 3.5, color = "black", hjust = 0)

  color_palette <- c(
    
    "CVPIA" = "salmon",
    "Flooded dynamic habitat" = "salmon",
    "other" = "grey90"
  )
  ggplot(data = s2) +
    geom_sf(aes(fill = map_1)) +
    scale_fill_manual(
      values = color_palette,
      name = "Category"
    ) +
    # labs(
    #   title = "Selected Habitats and Protected Areas",
    #   caption = "Blue: Selected Habitats, Green: Protected Areas, Purple: Overlap, Grey: None"
    # ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = -2200000, y = 2240000, label = paste0("Cost =", s2_costk,"\nFlooded areas=", flooded_areas_p2), size = 3.5, color = "black", hjust = 0)
  
  
  
  results_year<- st_drop_geometry(s2[, "solution_1"])
  results_tibble5<-bind_cols(results_year, pu_multi$p_ctgr5) 
  names(results_tibble5)<-c("Selected habitats", "Protected areas")
  pu_multi$selected_habitat<-results_year
  
  # Extract the solutions from prt_results as a matrix
  results_matrix5 <- as.matrix(results_tibble5)
  # Initialize an empty matrix to store pairwise kappa coefficients
  
  # Calculate Fleiss' Kappa coefficient matrix
  kappa_matrix5 <- kappam.fleiss(results_matrix5)
  
  # Define the Kappa coefficient value
  kappa_value5 <- 0.162  

  # Create a new column 'combined'
  pu_multi <- pu_multi %>%
    mutate(combined5 = case_when(
      `selected_habitat` == 1 & `p_ctgr5` == 1 ~ "overlap",
      `selected_habitat` == 1 ~ "selected_habitat",
      `p_ctgr5` == 1 ~ "protected (5%)",
      `selected_habitat` == 0 & `p_ctgr5` == 0 ~ "not selected"
    ))
  
  
  # Define the color palette with distinct colors
  # Define the color palette with distinct colors
  color_palette <- c(
    "selected_habitat" = "salmon",
    "protected (5%)" = "skyblue",
    "overlap" = "#FAE7AE",
    "not selected" = "grey90"
  )
  
  
  
  ggplot(data = pu_multi) +
    geom_sf(aes(fill = combined5)) +
    scale_fill_manual(
      values = color_palette,
      name = "Category"
    ) +
    # labs(
    #   title = "Selected Habitats and Protected Areas",
    #   caption = "Blue: Selected Habitats, Green: Protected Areas, Purple: Overlap, Grey: None"
    # ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) 
  pu_multi$
  tc_p <- eval_target_coverage_summary(p2, pu_multi[,'p_ctgr5'])
  ## visualize representation  (values show percent coverage)
  hist(
    tc_p$relative_held * 100,
    breaks=(max(tc_p$relative_held * 100)-min(tc_p$relative_held * 100))/2,
    col = "skyblue", 
    main = "Feature representation by prioritization",
    xlim = c(0, 100),
    ylim=c(0,400),
    xlab = "Percent coverage of features (%)"
  )
  
  pu_multi <- pu_multi %>%
    mutate(overlap = case_when(
      `combined` ==   "overlap"~1,
      `combined` ==  "selected_habitat"~ 0,
      `combined` =="pa_category"~ 0,
      `combined` == "not selected"~ 0,
      
    ))
  pu_multi <- pu_multi %>%
    mutate(overlap5 = case_when(
      `combined5` ==   "overlap"~1,
      `combined5` ==  "selected_habitat"~ 0,
      `combined5` =="protected (5%)"~ 0,
      `combined5` == "not selected"~ 0,
      
    ))

  tc_o <- eval_target_coverage_summary(p2, pu_multi[, "overlap5"])
  

  ## visualize representation  (values show percent coverage)
  hist(
    tc_o$relative_held * 100,
    breaks=(max(tc_o$relative_held * 100)-min(tc_o$relative_held * 100))/2,
    col = "#FAE7AE", 
    main = "Feature representation by prioritization",
    xlim = c(0, 100),
    #ylim=c(0,250),
    xlab = "Percent coverage of features (%)"
  )

  # Create a new column 'combined'
  pu_multi <- pu_multi %>%
    mutate(locked_in1 = case_when(
      `lockd_n` == 1 | `overlap` == 1 ~ TRUE,
      #`lockd_n` == 1 & `overlap` == 0 ~ TRUE,
      #`lockd_n` == 0 & `overlap` == 1 ~ TRUE,
      `lockd_n` == 0 & `overlap` == 0 ~ FALSE
    ))
  
  st_write(pu_multi, "data/pu_multi/pu_multi_cropareas_overlap.shp")
  pu_multi$locked_in<-as.logical(pu_multi$lockd_n)
  pu_multi$locked_in1<-as.logical(pu_multi$locked_in1)
  
  