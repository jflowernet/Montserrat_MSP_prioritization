################################
#Rasterize habitat layers
################################

habmap <- readOGR("data/benthic habitat/20161019_MNI_BenthicHabitatMap.shp")
habmap <- spTransform(habmap, mni_proj)

#map habitats if needed
spplot(habmap, "Habitat", col.regions = c("darkolivegreen4", "darkolivegreen3", "gray", "red4" , "coral", "goldenrod3", "yellow", "green4", "forestgreen"))

habitats <- unique(habmap$Habitat)

for(i in habitats) {
  print(paste0("Creating raster ",i))
  
  r <- rasterize(habmap[habmap$Habitat == i,], plangrid, field = 1)
  #output raster file
  writeRaster(r, filename = paste0("outputs/habitat_rasters/", i, ".tif"), format = "GTiff", datatype = "INT1U", overwrite = TRUE)
  
}