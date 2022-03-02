library("readxl")
library("janitor")
library("tidyverse")
#remotes::install_github("CIAT-DAPA/analogues")
library("analogues")
library("sf")
library("climatrends")
library("nasapower")

dat <- read_excel("data/EtNAM INFO.xls")

names(dat) <- make_clean_names(names(dat))

dat

p <- as.data.frame(dat[,c("long", "lat")])

tmean <- getData("worldclim", var = "tmean", res = 2.5, lon = p[1,1], lat = p[1,2])
prec <- getData("worldclim", var = "prec", res = 2.5, lon = p[1,1], lat = p[1,2])
eth <- getData("GADM", country = "ETH", level = 1)

tmean <- stack(crop(tmean, eth))
prec <- stack(crop(prec, eth))

plot(tmean[[1]])
plot(eth, add = TRUE)
points(p[,c("long","lat")], pch = "+")


r <- list()

for(i in seq_along(p[,1])){
  pars <- createParameters(x = p[i, 1], 
                           y = p[i, 2], 
                           vars = c("tmean","prec"),
                           weights = c(0.6,0.4),
                           ndivisions = c(12,12),
                           growing.season = c(6, 1),
                           rotation = "none",
                           threshold = 0.6,
                           env.data.ref = list(tmean, prec), 
                           env.data.targ = list(tmean, prec),
                           outfile="~/.",
                           fname = NA,
                           writefile=FALSE)
  
  sim <- calc_similarity(pars)
  
  r[[i]] <- sim
}

r2 <- stack(r)

r2

r2[r2[] < 0.3] <- 0
#r2[r2[] > 0.4] <- 1

r2

r2 <- calc(r2, fun = function(x){mean(x, na.rm = TRUE)})

r2

eth <- st_as_sf(eth)

r2 <- as.data.frame(r2, xy = TRUE)
r2 <- r2[!is.na(r2[, "layer"]), ]
r2 <- r2[!r2$layer < 0.2, ]

quantile(r2$layer)



# define colours for map of gradient of distribution
colpall <- colorRampPalette(c("#FFFF80", "#38E009",
                              "#1A93AB", "#0C1078"))

g <-
ggplot() + 
  #theme_void() +
  geom_sf(data = eth, color = "grey40", fill = "grey97") +
  geom_tile(r2, mapping = aes(x = x, y = y, fill = layer)) +
  geom_point(dat, mapping = aes(x = long, y = lat), 
             size = 3, col = "red", pch = 18) +
  scale_fill_gradientn(name = NULL, 
                       colours = colpall(10),
                       limits = c(0.2, 0.7),
                       breaks = seq(0.2, 0.7, 0.1),
                       labels = seq(0.2, 0.7, 0.1)) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c(.8, .7),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        panel.background = element_blank(),
        plot.margin = unit(rep(0, 4), "cm"))


g

ggsave(filename = "manuscript/similarity_map.png", 
       plot = g,
       dpi = 600,
       width = 20,
       height = 20,
       units = "cm")

ggsave(filename = "manuscript/similarity_map.pdf", 
       plot = g,
       dpi = 1000,
       width = 20,
       height = 20,
       units = "cm")



# now perform the characterization of trials during the experiment
temperature(dat[,c("long", "lat")], 
            day.one = dat$planting_date, 
            span = dat$harvesting_date - dat$planting_date)




