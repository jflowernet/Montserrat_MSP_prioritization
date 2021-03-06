richness <- read_csv("data/fishcoralDiversity.csv")

coordinates(richness) <- ~ LonDD + LatDD

#set projection
proj4string(richness) <- gps

#transform to same projection as other data
richness <- spTransform(richness, mni_proj)

#using planning grid as cost raster, setting all NA values to high cost so they won't be interpolated through - same as creating a barrier to the interpolation

costras <- plangrid
costras[costras == 0] <- 1
costras[is.na(costras)] <- 10000

species <- ipdw(richness, costras, range = 500, "Total_Diversity")

species[species == 0] <- NA

#scale 0-1
species <- species/cellStats(species, max)

writeRaster(species, filename = "outputs/speciesrich.tif", format = "GTiff", datatype = "FLT4S")