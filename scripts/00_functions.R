blood_slope <- function(input_object){
  
  # loop to get mean mass change per species
  species_list <- list()
  for(i in unique(input_object$species)){
    tmp <- input_object[input_object$species==i,]
    species <- unique(tmp$species)
    elev_range <- unique(tmp$elev_range)
    elev_guild  <- unique(tmp$elev_guild)
    avg_slope <- mean(tmp$slope)
    n <- nrow(tmp) # get sample size
    n.elev <- length(unique(tmp$elevation)) # get number of unique elevation localities
    mod.hb <- lm(hb ~ elevation, tmp) # simple regression for mass change
    mod.hb.sum <- summary(mod.hb)
    slope.hb <- mod.hb.sum$coefficients[2,1] # get hb slope
    r2.hb <- mod.hb.sum$r.squared # get hb r-squared
    se.hb <- mod.hb.sum$coefficients[2,2] # get hb SE
    mod.hct <- lm(hct ~ elevation, tmp) # simple regression for mass change
    mod.hct.sum <- summary(mod.hct)
    slope.hct <- mod.hct.sum$coefficients[2,1] # get hct slope
    r2.hct <- mod.hct.sum$r.squared # get hct r-squared
    se.hct <- mod.hct.sum$coefficients[2,2] # get hct SE
    species_list[[i]] <- cbind.data.frame(species,n,n.elev,slope.hb,r2.hb,se.hb,
                                          slope.hct,r2.hct,se.hct,
                                          elev_range,elev_guild,avg_slope)
  }
  
  # assemble dataframe
  data.blood.df <- do.call(rbind, species_list)
  rownames(data.blood.df) <- NULL
  colnames(data.blood.df) <- c("species","sample_size","unique_elevations","slope_hb","variance_hb","error_hb",
                               "slope_hct","variance_hct","error_hct","elev_range","elev_guild","avg_slope")
  return(data.blood.df)

}

convert_lat <- function(dataframe){
  lat_degrees <- dataframe$lat_degrees 
  lat_minutes <- dataframe$lat_minutes 
  lat <- lat_degrees + (lat_minutes/60)
  return(lat)
}

convert_long <- function(dataframe){
  long_degrees <- dataframe$long_degrees
  long_minutes <- dataframe$long_minutes
  long <- long_degrees + (long_minutes/60)
  return(long)
}



