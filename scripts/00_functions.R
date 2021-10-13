# approximation for water vapor pressure
water_vapor_pressure <- function(temp_celsius, humidity){
  temp_kelvin <- temp_celsius + 273.15
  wvp <- exp(20.386-(5132/temp_kelvin))*(humidity/100)
  return(wvp)
}

# coefficient of variation function
cv <- function(data){
  mean <- mean(data)
  sd <- sd(data)
  coef <- sd/mean
  return(coef)
}
  
# calculate slope and variance of blood parameters; species' range mean attributes
blood_slope <- function(input_object){
  
  # loop to get mean mass change per species
  species_list <- list()
  for(i in unique(input_object$species)){
    tmp <- input_object[input_object$species==i,]
    species <- unique(tmp$species)
    med_elev <- unique(tmp$elev_min) + unique(tmp$elev_range)/2
    elev_range <- unique(tmp$elev_range)
    sampling_range  <- (max(tmp$elevation) - min(tmp$elevation))/elev_range 
    elev_min <- unique(tmp$elev_min)
    elev_max <- unique(tmp$elev_max)
    mass <- mean(tmp$mass, na.rm = TRUE)
    n <- nrow(tmp) # get sample size
    n.elev <- length(unique(tmp$elevation)) # get number of unique elevation localities
    mod.hb <- lm(hb ~ elevation, tmp) # simple regression for mass change
    mod.hb.sum <- summary(mod.hb)
    slope.hb <- mod.hb.sum$coefficients[2,1] # get hb slope
    r2.hb <- mod.hb.sum$r.squared # get hb r-squared
    se.hb <- mod.hb.sum$coefficients[2,2] # get hb SE
    mod.hct <- lm(hct_percent ~ elevation, tmp) # simple regression for mass change
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
                                          elev_range,sampling_range,elev_min,elev_max,med_elev,mass)
  }
  
  # assemble dataframe
  data.blood.df <- do.call(rbind, species_list)
  rownames(data.blood.df) <- NULL
  colnames(data.blood.df) <- c("species","sample_size","unique_elevations","slope_hb","r2_hb","error_hb",
                               "slope_hct","r2_hct","error_hct","slope_mchc","r2_mchc","error_mchc",
                               "elev_range","sampling_range","elev_min","elev_max","median_elevation","mass")
  return(data.blood.df)

}

# calculate variance of blood parameters + species' range mean attributes
blood_variance <- function(input_object, min_bin){ 

  # loop to get mean mass change per species
  species_list <- list()
  for(i in unique(input_object$species)){
    #print(i) 
    tmp <- input_object[input_object$species==i,]
    tt <- table(tmp$binID)
    best_bins <- tt[tt>min_bin] # use only the 100 m elevation bins passing the data cutoff
    tt <- tt[tt>min_bin] %>% as.data.frame() 
    if(dim(tt)[2]==1){ # need to add an extra column if there is only one bin
      tt$binID <- rownames(tt) 
      colnames(tt) <- c("Freq","binID")
    } else { # otherwise we just need to assign columns
      colnames(tt) <- c("binID","Freq")
    }
    tmp <- tmp[tmp$binID %in% names(best_bins),]
    tmp <- merge(tmp, tt, by.x = "binID", by.y = "binID")
    if(nrow(tt)>0){
      include <- "yes"
      range_position <- aggregate(range_position ~ binID, tmp, mean) %>% select(range_position) %>% as.vector() # realtive position of bin
      edge_distance <- aggregate(edge_distance ~ binID, tmp, mean) %>% select(edge_distance) %>% as.vector() # edge distance of bin
      var.hb <- aggregate(hb ~ binID, tmp, cv) %>% select(hb) %>% as.vector() # hb variance per bin
      var.hct <- aggregate(hct ~ binID, tmp, cv) %>% select(hct) %>% as.vector() # hct variance per bin
      var.mchc <- aggregate(MCHC_calculated ~ binID, tmp, cv) %>% select(MCHC_calculated) %>% as.vector() # mchc variance per bin
      species <- unique(tmp$species) %>% rep(., nrow(var.hb)) # species ID
      med_elev <- unique(tmp$elev_min) + unique(tmp$elev_range)/2  %>% rep(., nrow(var.hb)) # median species range elevation
      bin_elev <- aggregate(elevation ~ binID, tmp, mean) %>% select(elevation) %>% as.vector() # mean bin elevation
      elev_range <- unique(tmp$elev_range)  %>% rep(., nrow(var.hb))  # mean mass of species
      mass <- mean(tmp$mass, na.rm = TRUE)  %>% rep(., nrow(var.hb)) # mean mass of species
      n <- aggregate(Freq ~ binID, tmp, mean) %>% select(Freq) %>% as.vector() # sample size per bin
      n.elev <- length(unique(tmp$elevation))  %>% rep(., nrow(var.hb)) # number of unique elevation localities per species
    } else {
      include <- "no"
      }
    if(include=="yes")
      {species_list[[i]] <- cbind.data.frame(species,n,n.elev,range_position,
                                             edge_distance,var.hb,var.hct,
                                          var.mchc,elev_range,med_elev,bin_elev,mass)}
  }
    
    # assemble dataframe
    data.blood.df <- do.call(rbind, species_list)
    rownames(data.blood.df) <- NULL
    colnames(data.blood.df) <- c("species","sample_size","unique_elevations","range_position","edge_distance",
                                 "variance_hb","variance_hct","variance_mchc","elev_range","median_elevation","bin_elevation",
                                 "mass")
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


credibility_coder <- function(dataframe){
  dataframe$credible <- 0
  df <- list()
  for(i in unique(dataframe$.variable)){
    sub <- dataframe[dataframe$.variable==i,]
    sub <- sub %>% mutate(levels = case_when(
      .lower < 0 & .upper < 0 & .width==0.95 | .lower > 0 & .upper > 0 & .width==0.95  ~ 1,
      .lower < 0 & .upper < 0 & .width==0.80 | .lower > 0 & .upper > 0 & .width==0.80  ~ 2, 
      .lower < 0 & .upper > 0 | .lower > 0 & .upper < 0 ~ 0,
      .width==0.5 ~ 0
    ))
    if(1 %in% sub$levels){sub$credible <- 1} 
    else if(2 %in% sub$levels){sub$credible <- 2} 
    else {sub$credible <- 0} 
    df[[i]] <- sub
  }
  new_df <- do.call(rbind, df)
  return(new_df)
}

