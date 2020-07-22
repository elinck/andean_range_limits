# calculate slope and variance of blood parameters; species' range mean attributes
blood_slope <- function(input_object){
  
  # loop to get mean mass change per species
  species_list <- list()
  for(i in unique(input_object$species)){
    tmp <- input_object[input_object$species==i,]
    species <- unique(tmp$species)
    med_elev <- unique(tmp$elev_min) + unique(tmp$elev_range)/2
    elev_range <- unique(tmp$elev_range)
    mass <- mean(tmp$mass, na.rm = TRUE)
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
    mod.mchc <- lm(MCHC_calculated ~ elevation, tmp) # simple regression for mass change
    mod.mchc.sum <- summary(mod.mchc)
    slope.mchc <- mod.mchc.sum$coefficients[2,1] # get hct slope
    r2.mchc <- mod.mchc.sum$r.squared # get hct r-squared
    se.mchc <- mod.mchc.sum$coefficients[2,2] # get hct SE
    species_list[[i]] <- cbind.data.frame(species,n,n.elev,slope.hb,r2.hb,se.hb,
                                          slope.hct,r2.hct,se.hct,slope.mchc,r2.mchc,se.mchc,
                                          elev_range,med_elev,mass)
  }
  
  # assemble dataframe
  data.blood.df <- do.call(rbind, species_list)
  rownames(data.blood.df) <- NULL
  colnames(data.blood.df) <- c("species","sample_size","unique_elevations","slope_hb","r2_hb","error_hb",
                               "slope_hct","r2_hct","error_hct","slope_mchc","r2_mchc","error_mchc",
                               "elev_range","median_elevation","mass")
  return(data.blood.df)

}

# calculate variance of blood parameters + species' range mean attributes
blood_variance <- function(input_object){
  
  # loop to get mean mass change per species
  species_list <- list()
  for(i in unique(input_object$species)){
    tmp <- input_object[input_object$species==i,]
    species <- unique(tmp$species)
    med_elev <- unique(tmp$elev_min) + unique(tmp$elev_range)/2
    elev_range <- unique(tmp$elev_range)
    mass <- mean(tmp$mass, na.rm = TRUE)
    n <- nrow(tmp) # get sample size
    n.elev <- length(unique(tmp$elevation)) # get number of unique elevation localities
    tt <- table(tmp$binID)
    best_bin <- names(tt[which.max(tt)])
    tmp <- tmp[tmp$binID==best_bin,]
    range_position <- mean(tmp$range_position)
    edge_distance <- mean(tmp$edge_distance)
    var.hb <- var(tmp$hb)
    var.hct <- var(tmp$hct)
    var.mchc <- var(tmp$MCHC_calculated)
    if(nrow(tmp)>4){species_list[[i]] <- cbind.data.frame(species,n,n.elev,range_position,edge_distance,var.hb,var.hct,
                                          var.mchc,elev_range,med_elev,mass)}
  }
    
    # assemble dataframe
    data.blood.df <- do.call(rbind, species_list)
    rownames(data.blood.df) <- NULL
    colnames(data.blood.df) <- c("species","sample_size","unique_elevations","range_position","edge_distance",
                                 "variance_hb","variance_hct","variance_mchc","elev_range","median_elevation","mass")
    return(data.blood.df)
}

# convert latitudes to decimal degrees
convert_lat <- function(dataframe){
  lat_degrees <- dataframe$lat_degrees 
  lat_minutes <- dataframe$lat_minutes 
  lat <- lat_degrees + (lat_minutes/60)
  return(lat)
}

# convert longitudes to decimal degrees
convert_long <- function(dataframe){
  long_degrees <- dataframe$long_degrees
  long_minutes <- dataframe$long_minutes
  long <- long_degrees + (long_minutes/60)
  return(long)
}

# drop outliers based on cook's distance, whether individual records at range limits are far from others
outliers_cooks <- function(dataframe, param, cutoff1=4, cutoff2=3.5){
  newdf <- list()
  for(i in dataframe$species){
    tmp <- dataframe[dataframe$species==i,]
    min <- min(unique(tmp$elevation))
    max <- max(unique(tmp$elevation))
    range <- max-min
    range_quant <- range/4
    sample_size <- nrow(tmp)
    lm1 <- lm(tmp[, param] ~ tmp[, "elevation"])
    cooksd <- cooks.distance(lm1)
    influential_4 <- as.numeric(names(cooksd)[(cooksd > (cutoff1/sample_size))]) 
    influential_3 <- as.numeric(names(cooksd)[(cooksd > (cutoff2/sample_size))]) 
    tmp$include <- "yes"
    tmp$include[row.names(tmp) %in% influential_3 & !grepl("no",tmp$bursa) | 
                  row.names(tmp) %in% influential_4] <- "no"
    tmp$include[tmp$elevation==max & nrow(tmp[tmp$elevation>max-range_quant,])<2] <- "no"
    tmp$include[tmp$elevation==min & nrow(tmp[tmp$elevation<min+range_quant,])<2] <- "no"
    newdf[[i]] <- tmp[!tmp$include=="no",]
  }
  passing <- do.call(rbind, newdf)
  return(passing)
}

# drop species that don't hit minimum sample criteria
outliers_limits <- function(dataframe, min_sample, min_limit, min_range){
  elev_cutoff <- c()
  for(i in dataframe$species){
    tmp <- dataframe[dataframe$species==i,]
    num <- length(unique(tmp$elevation))
    min <- min(unique(tmp$elevation))
    max <- max(unique(tmp$elevation))
    sample_size <- nrow(tmp)
    range <- max-min
    range_quant <- range/4
    if(num >= min_sample  & range > min_range  &
       nrow(tmp[tmp$elevation<min+range_quant,])>=min_limit & 
       nrow(tmp[tmp$elevation>max-range_quant,])>=min_limit)
    {elev_cutoff[i] <- as.character(tmp$species[1])}
  }
  elev_cutoff <- as.vector(elev_cutoff)
  dataframe <- dataframe[dataframe$species %in% elev_cutoff,]
  return(dataframe)
}
