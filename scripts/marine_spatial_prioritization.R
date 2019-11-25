#########################################################
#Montserrat marine spatial prioritization - code to accompany paper Flower et al. 2019
#Purpose: Generate map of priority conservation areas within Montserrat's 100m shelf area
#Data from Waitt Institute Scientific Assessment of Montserrat conducted in October 2015
#Contact: Jason Flower, Sustainable Fisheries Group, University of California Santa Barbara (jflower@ucsb.edu)
#Date: 31 July 2019
#########################################################

#load required libraries
#Note: Use of the Gurobi library requires a license. More information about installing Gurobi https://prioritizr.net/articles/gurobi_installation.html
#Other optimizers can be used to solve the prioritizr conservation problems: https://prioritizr.net/articles/saltspring.html#solving-the-problem

# pacman will help us install any necessary packages
if (!require("pacman")) install.packages("pacman")
# pacman::p_load checks to see if these packages are installed, and installs them if not
pacman::p_load(readr, dplyr, raster, rgdal, sp, sf, prioritizr, gurobi, ipdw, tmap)

########################################################
#Set projections
########################################################

#projection for raw data (i.e. GPS points)
gps <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#set projection for Montserrat 1958 British West Indies grid. R does not set the projection correctly when it imports spatial data from this projection.
mni_proj <- crs("+proj=tmerc +lat_0=0 +lon_0=-62 +k=0.9995000000000001 +x_0=400000 +y_0=0 +ellps=clrk80 +towgs84=174,359,365,0,0,0,0 +units=m +no_defs")

########################################################
#Create cost layers
########################################################

#import planning grid = MNI 100m shelf, rasterized at 100m x 100m resolution
plangrid <- raster("data/plangrid")
crs(plangrid) <- mni_proj

#reclassify values to 0
plangrid[plangrid == 1] <- 0

#read in fish pots 100m buffered and dissolved shapefile
pots <- readOGR("data/fishpots/")
proj4string(pots) <- mni_proj

#create raster from the pots polygon, setting all fish pot locations to value of 1
potscost <- rasterize(pots, plangrid, field = 1, update =TRUE)

#import fishing value layer masked to 100m shelf and extent set to same as planning grid in ArcGIS
fishingvalue <- raster("data/fishing_value")
projection(fishingvalue) <- mni_proj

#scale all values 0-1 by dividing by max value using log-transformed +1 fishing values (+1 so that zero values are transformed to zero not NA)
fishing<- log(fishingvalue+1)/ cellStats(log(fishingvalue+1), max)

#add fish pots layer to fishing value layer to get a cost layer
cost <- potscost + fishing

#reclassify zero cost values to 0.01 as prioritzr won't work with zeros in cost file
cost <- reclassify(cost, c(-0.1, 0.1, 0.01))

##################################################
#Import habitat layers
#################################################

filelist <- list.files("data/habitat_rasters/", ".tif", full.names = TRUE)
habitatstack <- stack(filelist)

#For the purposes of the prioritization, we combined all algal reefs (hard bottom and mixed bottom) to one habitat. Sand and Hard botom and sand combined to another single habitat category

algal <- merge(habitatstack$Algal_Reef_.Hard_Bottom., habitatstack$Algal_Reef_.Mixed_Bottom.)
nonliving <- merge(habitatstack$Hard_Bottom_and_Sand, habitatstack$Sand)

habitatstack_simplified <- stack(algal, nonliving, habitatstack$Artificial_Reef, habitatstack$Colonized_Volcanic_Boulders, habitatstack$Coral_Reef, habitatstack$Sargassum_Forest, habitatstack$Seagrass)

#################################################
#Create species richness layer 
#################################################

#script for import and interpolation of site species richness data
source("scripts/interp_sprichness.R")

species <- raster("outputs/speciesrich.tif")

#############################################
#Spatial prioritization using prioritizr
#############################################

#create a stack with the habitat and species richness data
features_stack<- stack(habitatstack_simplified, species)

#create the conservation problem - protect 30% of each habitat, 50% of total species richness, while minimizing cost (overlap with fishing effort proxy)
prob <- problem(cost, features = features_stack) %>%
  add_min_set_objective() %>%
  add_relative_targets(c(rep(0.3, raster::nlayers(habitatstack_simplified)),0.5)) %>%
  add_binary_decisions() %>% 
  add_boundary_penalties(penalty = .001, edge_factor = 0.5)

#solve problem
sprob <- solve(prob)

plot(sprob,  main = c("30% of 7 habitats, 50% species richness (interpolated), penalty = 0.001, edge = 0.5"))

#################################################
#10 - 40% area protection targets comparison - results for Appendix Figure S3
#################################################

#read in the final zones map - this will be used here and in the locked-in no-take areas code chunk at the end of this script
zoning_plan <- st_read("data/2018_03_20_Option_5/")

#filter for only the no-take zones and transform to same projection as rest of spatial data
no_take_MPAs <- zoning_plan %>% 
  filter(ZONE_TYPE == "Sanctuary") %>% 
  st_transform(mni_proj)

#create the conservation problem - protect 10, 20, 30 or 40% of each habitat, 50% of total species richness, while minimizing cost (overlap with fishing effort proxy)
prob_10perc <- problem(cost, features = features_stack) %>%
  add_min_set_objective() %>%
  add_relative_targets(c(rep(0.1, raster::nlayers(habitatstack_simplified)),0.5)) %>%
  add_binary_decisions() %>% 
  add_boundary_penalties(penalty = .001, edge_factor = 0.5)

prob_20perc <- problem(cost, features = features_stack) %>%
  add_min_set_objective() %>%
  add_relative_targets(c(rep(0.2, raster::nlayers(habitatstack_simplified)),0.5)) %>%
  add_binary_decisions() %>% 
  add_boundary_penalties(penalty = .001, edge_factor = 0.5)

prob_30perc <- problem(cost, features = features_stack) %>%
  add_min_set_objective() %>%
  add_relative_targets(c(rep(0.3, raster::nlayers(habitatstack_simplified)),0.5)) %>%
  add_binary_decisions() %>% 
  add_boundary_penalties(penalty = .001, edge_factor = 0.5)

prob_40perc <- problem(cost, features = features_stack) %>%
  add_min_set_objective() %>%
  add_relative_targets(c(rep(0.4, raster::nlayers(habitatstack_simplified)),0.5)) %>%
  add_binary_decisions() %>% 
  add_boundary_penalties(penalty = .001, edge_factor = 0.5)

#solve problem
sprob_10perc <- solve(prob_10perc)
sprob_20perc <- solve(prob_20perc)
sprob_30perc <- solve(prob_30perc)
sprob_40perc <- solve(prob_40perc)

plotstack <- stack(sprob_10perc, sprob_20perc, sprob_30perc, sprob_40perc)

#plot outputs
tm_shape(plotstack)+
tm_raster(palette = c("#c6c5c5", "#409a00"), n=2, legend.show = FALSE) +
  tm_layout(title = c("(a)", "(b)", "(c)", "(d)")) +
  tm_shape(no_take_MPAs) +
  tm_borders("black")

tmap_save(filename = "outputs/varying_habitat_targets.png")

########################################################
#calculate percent of total no-take MPA area occupied by cells selected by prioritizr

#create single MPA layer for calculation
no_take_MPAs_dissolved <- no_take_MPAs %>% 
  group_by() %>% 
  summarise()

#percent for each habitat target scenario: 20, 30, 40, 50%
table_habitat_target_area_priortized <- extract(plotstack, no_take_MPAs_dissolved, fun=function(x, ...) sum(na.omit(x))/length(na.omit(x)))
colnames(table_habitat_target_area_priortized) <- c("20% target", "30%", "40%", "50%")
table_habitat_target_area_priortized

#percent for each habitat target scenario, broken down into each of the 4 MPA zones - see 
table_habitat_target_area_priortized_byzone <- extract(plotstack, no_take_MPAs, fun=function(x, ...) sum(na.omit(x))/length(na.omit(x)))
colnames(table_habitat_target_area_priortized_byzone) <- c("20% target", "30%", "40%", "50%")
rownames(table_habitat_target_area_priortized_byzone) <- pull(no_take_MPAs, NAME)
table_habitat_target_area_priortized_byzone

#################################################
#Varying species richness feature target - results for Appendix Figure S2
#################################################

#create the conservation problem - protect 10, 30, 50 or 70% of total (summed) species richness, 30% of each habitat, while minimizing cost (overlap with fishing effort proxy)
species_prob_10perc <- problem(cost, features = features_stack) %>%
  add_min_set_objective() %>%
  add_relative_targets(c(rep(0.3, raster::nlayers(habitatstack_simplified)), 0.1)) %>%
  add_binary_decisions() %>% 
  add_boundary_penalties(penalty = .001, edge_factor = 0.5)

species_prob_30perc <- problem(cost, features = features_stack) %>%
  add_min_set_objective() %>%
  add_relative_targets(c(rep(0.3, raster::nlayers(habitatstack_simplified)), 0.3)) %>%
  add_binary_decisions() %>% 
  add_boundary_penalties(penalty = .001, edge_factor = 0.5)

species_prob_50perc <- problem(cost, features = features_stack) %>%
  add_min_set_objective() %>%
  add_relative_targets(c(rep(0.3, raster::nlayers(habitatstack_simplified)), 0.5)) %>%
  add_binary_decisions() %>% 
  add_boundary_penalties(penalty = .001, edge_factor = 0.5)

species_prob_70perc <- problem(cost, features = features_stack) %>%
  add_min_set_objective() %>%
  add_relative_targets(c(rep(0.3, raster::nlayers(habitatstack_simplified)), 0.7)) %>%
  add_binary_decisions() %>% 
  add_boundary_penalties(penalty = .001, edge_factor = 0.5)

#solve problem
s_species_prob_10perc <- solve(species_prob_10perc)
s_species_prob_30perc <- solve(species_prob_30perc)
s_species_prob_50perc <- solve(species_prob_50perc)
s_species_prob_70perc <- solve(species_prob_70perc)

plotstack <- stack(s_species_prob_10perc, s_species_prob_30perc, s_species_prob_50perc, s_species_prob_70perc)

#plot outputs
tm_shape(plotstack)+
  tm_raster(palette = c("#c6c5c5", "#409a00"), n=2, legend.show = FALSE) +
  tm_layout(title = c("(a)", "(b)", "(c)", "(d)")) 

tmap_save(filename = "outputs/varying_sprich_targets.png")

#####################################################
#Lock-in no-take areas and re-run prioritization - results for Appendix Figure S5
#####################################################

#rasterize the no-take areas polygons from the final zoning plan
no_take_MPAs_raster <- rasterize(no_take_MPAs, plangrid, field = 1, update =TRUE)

#create the conservation problem - protect 30% of each habitat, 50% of total species richness, while minimizing cost (overlap with fishing effort proxy), with areas selected as no-take reserves in the final draft zoning plan locked-in as protected
#we only need to take the original conservation problem defined previously and add the locked-in constraints

prob_with_lock_in <- prob %>% add_locked_in_constraints(no_take_MPAs_raster)

#solve problem
sprob_lock_in <- solve(prob_with_lock_in)

#plot output with no-take area borders
tm_shape(sprob_lock_in)+
  tm_raster(palette = c("#c6c5c5", "#409a00"), n=2, legend.show = FALSE) +
tm_shape(no_take_MPAs) +
  tm_borders("black")

tmap_save(filename = "outputs/lock_in_no_takes.png")
