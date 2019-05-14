#########################################################
#Montserrat marine spatial prioritization - code to accompany paper Flower et al. 2019
#Purpose: Generate map of priority conservation areas within 100m shelf
#Data from Waitt Institute Scientific Assessment of Montserrat conducted in October 2015
#Contact: Jason Flower, Sustainable Fisheries Group, University of California Santa Barbara (jflower@ucsb.edu)
#Date: 13 May 2019
#########################################################

#load required libraries
#Note: Use of the Gurobi library requires a license. More information about installing Gurobi https://prioritizr.net/articles/gurobi_installation.html

library(readr)
library(dplyr)
library(raster)
library(rgdal)
library(sp)
library(prioritizr)
library(gurobi) 
library(tmap)

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
#10 - 40% area protection targets comparison
#################################################

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
  tm_layout(title = c("(a) 10%", "(b) 20%", "(c) 30%", "(d) 40%")) 

tmap_save(filename = "outputs/varying_habitat_targets.png")
#################################################
#Varying species richness feature target
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
  tm_layout(title = c("(a) 10%", "(b) 30%", "(c) 50%", "(d) 70%")) 

tmap_save(filename = "outputs/varying_sprich_targets.png")