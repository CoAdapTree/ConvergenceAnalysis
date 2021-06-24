##pooja.singh09@gmail.com
##March2021
##importing climate data from worldclim: https://www.worldclim.org/data/bioclim.html

library(raster)
library(sp)
library(rgeos)
library(rgdal)
library(sf)
setwd("/Users/pooja/Desktop/Postdoc_Calgary_2019/Research_projects/CoAdapTree/Convergence/environmental_data/all_species_2021/worldclim")


worldclim <- getData("worldclim",var="bio",res=5)
#par(oma = c(0.1, 0.1, 0.1, 2.1))
#plot(worldclim[[c("bio1", "bio12")]])

names(worldclim) <-  c("AMT", "MDR", "ISO", "TSN", "MWMT", "MCMT", "TAR", "MWeQT", "MDQT", "MWQT", "MCQT", "AP", "PWM", "PDM", "PS", "PWeQ", "PDQ", "PWQ", "PCQ")

##plotworld
gain(worldclim)=c(rep(0.1,11),rep(1,7))

pdf("worldclim_all.pdf", h=8, w=8)
plot(worldclim)
dev.off()

##plot NA+Eurasia
r1 <- raster::crop(worldclim, extent(-180,150,20,100))
pdf("worldclim_all_cropped.pdf", h=8, w=8)
plot(r1)
dev.off()


### Psylvestris
a <- read.csv("/Users/pooja/Desktop/Postdoc_Calgary_2019/Research_projects/gentree/envs/Psylvestris_locations.csv", header=T)
psyl <- a[,c(4,3)] ##NB: input data
ex <- extract(worldclim,psyl)
psyl1 <- cbind(psyl, ex)
psyl1$ID1 <- "Psyl"
psyl1$ID2 <- a$ID2
write.table(psyl1, "Psylvestris_locations_woldclim.txt", quote=F, sep="\t")


##psyl_points=raster::extract(worldclim,psyl,df=T)  ## another method to get climates

gplot(worldclim[[1]])+
  geom_raster(aes(fill=value))+
  geom_point(
    data=as.data.frame(psyl),
    aes(x=x,y=y),col="red")+
  coord_equal()

##Pabies

a <- read.csv("/Users/pooja/Desktop/Postdoc_Calgary_2019/Research_projects/gentree/envs/Pabies_locations.csv", header=T)
pabi <- a[,c(4,3)]
ex <- extract(worldclim,pabi)
pabi1 <- cbind(pabi, ex)
pabi1$ID1 <- "Pabi"
pabi1$ID2 <- a$ID2
write.table(pabi1, "Pabies_locations_woldclim.txt", quote=F, sep="\t")



###Pobovata


a <- read.csv("/Users/pooja/Desktop/Postdoc_Calgary_2019/Research_projects/gentree/envs/Pobovata_locations.csv", header=T)
pobo <- a[,c(4,3)]
ex <- extract(worldclim,pobo)
pobo1 <- cbind(pobo, ex)
pobo1$ID1 <- "Pobo"
pobo1$ID2 <- a$ID2
write.table(pobo1, "Pobovata_locations_woldclim.txt", quote=F, sep="\t")




##### misc





#getData('worldclim', var="bio1", res = 10,lon=37.76667, lat=-2.816667)

## crop to a latitude/longitude box
r1 <- raster::crop(clim[[1]], extent(10,35,-35,-20))
## Crop using a Spatial polygon
r1 <- raster::crop(clim[[1]], bbox(za))







BIO1 = Annual Mean Temp. [AMT]

BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)) [MDR]

BIO3 = Isothermality (BIO2/BIO7) (×100) [ISO]

BIO4 = Temp. Seasonality (standard deviation ×100) [TSN]

BIO5 = Max Temp. of Warmest Month [MWMT]

BIO6 = Min Temp. of Coldest Month [MCMT]

BIO7 = Temp. Annual Range (BIO5-BIO6) [TAR]

BIO8 = Mean Temp. of Wettest Quarter [MWeQT]

BIO9 = Mean Temp. of Driest Quarter [MDQT]

BIO10 = Mean Temp. of Warmest Quarter [MWQT]

BIO11 = Mean Temp. of Coldest Quarter [MCQT]

BIO12 = Annual Precipitation [AP]

BIO13 = Precipitation of Wettest Month [PWM]

BIO14 = Precipitation of Driest Month [PDM]

BIO15 = Precipitation Seasonality (Coefficient of Variation) [PS]

BIO16 = Precipitation of Wettest Quarter [PWeQ]

BIO17 = Precipitation of Driest Quarter [PDQ]

BIO18 = Precipitation of Warmest Quarter [PWQ]

BIO19 = Precipitation of Coldest Quarter [PCQ]


