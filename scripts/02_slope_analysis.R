library(sp)  # classes for spatial data
library(raster)  # grids, rasters
library(rasterVis)  # raster visualisation
library(maptools)
library(rgeos)
library(ggplot2)
library(viridis)
library(ggthemes)
library(cowplot)

peru_elev <- getData("alt", country = "PER", )
peru_terrain <- terrain(peru_elev, opt = c("slope", "aspect"), unit = "degrees")
peru_df  <- as.data.frame(peru_terrain, xy = TRUE)
peru_df$country <- "Peru"

ecuador_elev <- getData("alt", country = "ECU")
ecuador_terrain <- terrain(ecuador_elev, opt = c("slope", "aspect"), unit = "degrees")
ecuador_df  <- as.data.frame(ecuador_terrain, xy = TRUE)
ecuador_df <- ecuador_df[ecuador_df$x>-82,]
ecuador_df$country <- "Ecuador"

colombia_elev <- getData("alt", country = "COL")
colombia_terrain <- terrain(colombia_elev, opt = c("slope", "aspect"), unit = "degrees")
colombia_df  <- as.data.frame(colombia_terrain, xy = TRUE)
colombia_df$country <- "Colombia"

bolivia_elev <- getData("alt", country = "BOL")
bolivia_terrain <- terrain(bolivia_elev, opt = c("slope", "aspect"), unit = "degrees")
bolivia_df  <- as.data.frame(bolivia_terrain, xy = TRUE)
bolivia_df$country <- "Bolivia"

venezuela_elev <- getData("alt", country = "VEN")
venezuela_terrain <- terrain(venezuela_elev, opt = c("slope", "aspect"), unit = "degrees")
venezuela_df  <- as.data.frame(venezuela_terrain, xy = TRUE)
venezuela_df$country <- "Venezuela"


master_df <- rbind.data.frame(venezuela_df,colombia_df,ecuador_df, peru_df, bolivia_df)
master_df$country <- as.factor(master_df$country)
master_df$country <- factor(master_df$country, levels = c("Venezuela","Colombia","Ecuador","Peru","Bolivia"))

plot <- ggplot(master_df,aes(x=x, y=y, fill=slope)) +
  geom_raster() +
  scale_fill_viridis(na.value = "white", option="cividis") +
  coord_quickmap() +
  theme_bw() +
  facet_wrap(~country, scales="free") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill="Slope (°)") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

png("~/Desktop/pendiente_andina.png", width=8, height=6.5, units="in", res=500)
plot
dev.off()

ggplot(peru_df, aes(x=slope)) +
  theme(
    panel.grid = element_blank(),
  )  +
  theme_classic() +
  geom_histogram(binwidth = 1, fill="forestgreen",color="black",alpha=0.5) +
  xlab("Slope (°)") +
  ylab("Count")
  
