#### Code for simulating uncertainty in dems using geostatistical analysis of errors.
#### inputs needed: an original dem (uased as a grid template) and a set of measured errors as points
#### this version: 2015_03_31, based on wageningen model


rm(list = ls())  # clean memory

# install and Load libraries
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("gstat", "sp","rgdal", "RColorBrewer",
              "ggplot2", "automap","maptools")

ipak(packages)



setwd("/Volumes/MacData/odrive/GDrive_UC/geo/Gualaceo/AOI_R/Inputs")


# read error Data
herror <- readOGR(".", "gnd_pts_AOIR")

# a random sample subset
herror <- herror[sample(1:nrow(herror), 500,replace=FALSE),]

str(herror)
dim(herror)
names(herror)
summary(herror$res_st_gt)
proj4string(herror) <- CRS("+init=epsg:32717")
# histogram
hist(herror$res_st_gt, main="DEM sample error", nclass=25)

# make spatial if not
#coordinates(herror) = ~x+y
#proj4string(herror) <- CRS("+init=epsg:32717")



# read boundary study area
aoi = readShapePoly("AOI_R.shp")
str(aoi)

# bubble plot
bubble(herror, zcol="res_st_gt",
       maxsize=3, 
       sp.layout=list("sp.polygons", aoi))

#geostatistical analysis

# define gstat object and 
gs_herror = gstat(id = c("res_st_gt"), formula = res_st_gt~1, data = herror)

# compute experimental variogram
#vgs_herror = variogram(gs_herror,boundaries = c(2,4,5,10,15,20,40,60,80,100,150,200))
vgs_herror = variogram(gs_herror,cutoff = 100, width=6, map=F)
plot(vgs_herror, plot.numbers = TRUE)

###USING AUTOVARIOGRAM
#compute automatic variogram from automap library
autovg <- autofitVariogram(res_st_gt~1,herror)
autovg
plot(autovg)

autocv <- autoKrige.cv(res_st_gt~1, herror)
summary(autocv)

###IDW
krige.cv <- krige.cv(res_st_gt~1,herror)
summary(krige.cv)
cor(krige.cv$var1.pred,krige.cv$observed)

###exploring and defining variongrams
#List of variogram models
show.vgms()

# model initial variogram based on visual matching, try different models
vgm_herror_sph <- vgm(nugget=0.4, psill=0.6, range=40, model="Sph")
plot(vgs_herror,vgm_herror_sph, main="Variogram Model Sph")

vgm_herror_ste <- vgm(nugget=0.64, psill=0.4, range=26, kappa=10, model="Ste")
plot(vgs_herror,vgm_herror_ste, main="Variogram Model Ste")

vgm_herror_gau <- vgm(nugget=0, psill=0.46, range=50, model="Gau")

plot(vgs_herror,vgm_herror_gau,plot.numbers = TRUE, xlim=c(0,500), main="Variogram Model Exp")


# fit a gau variogram to experimental variogram
vgf_herror_gau <- fit.variogram(vgs_herror,vgm_herror_gau)
plot(vgs_herror,vgf_herror_gau, main="gau Variogram fitted")
vgm_herror_gau
vgf_herror_gau

# fit a sph variogram to experimental variogram
vgf_herror_sph <- fit.variogram(vgs_herror,vgm_herror_sph)
plot(vgs_herror,vgf_herror_sph, main="sph Variogram fitted")
vgm_herror_sph
vgf_herror_sph

# fit a ste variogram to experimental variogram
vgf_herror_ste <- fit.variogram(vgs_herror,vgm_herror_ste)
plot(vgs_herror,vgf_herror_ste, main="ste Variogram fitted")
vgm_herror_ste
vgf_herror_ste

#crossvalidation of different models. Zscores' mean and std should be near 0 and 1.
vgm_herror_gau.xv = krige.cv(res_gt_st~1, locations=herror, model=vgf_herror_gau,nmax=20, nfold=5)
hist(vgm_herror_gau.xv$var1.pred)
mean(vgm_herror_gau.xv$zscore, na.rm=T)
sd(vgm_herror_gau.xv$zscore, na.rm=T)

vgm_herror_ste.xv = krige.cv(res_gt_st~1, locations=herror, model=vgf_herror_ste,nmax=20, nfold=5)
hist(vgm_herror_ste.xv$var1.pred)
mean(vgm_herror_ste.xv$zscore, na.rm=T)
sd(vgm_herror_ste.xv$zscore, na.rm=T)

vgm_herror_sph.xv = krige.cv(res_gt_st~1, locations=herror, model=vgf_herror_sph,nmax=20, nfold=5)
hist(vgm_herror_sph.xv$var1.pred)
mean(vgm_herror_sph.xv$zscore, na.rm=T)
sd(vgm_herror_sph.xv$zscore, na.rm=T)

# compare results of crossvalidation.
#Spherical model is selected
xvresult <- rbind(c("model","mean","sd"))
xvresult <- rbind(xvresult, c("gau",mean(vgm_herror_gau.xv$zscore, na.rm=T),sd(vgm_herror_gau.xv$zscore, na.rm=T)))
xvresult <- rbind(xvresult, c("ste",mean(vgm_herror_ste.xv$zscore, na.rm=T),sd(vgm_herror_ste.xv$zscore, na.rm=T)))
xvresult <- rbind(xvresult, c("sph",mean(vgm_herror_sph.xv$zscore, na.rm=T),sd(vgm_herror_sph.xv$zscore, na.rm=T)))

hist(vgm_herror_ste.xv$zscore,nclass=15 )
hist(vgm_herror_ste.xv$residual,nclass=20)

bubble(vgm_herror_ste.xv, zcol="residual",
       maxsize=2,
       sp.layout=list("sp.polygons", aoi), 
       main = "cross-validation residuals")


#dem used as template
dem <- readGDAL("ST_DEM_AOIr.tif")
proj4string(dem) <- CRS("+init=epsg:32717")
spplot(dem)

# local kriging using STE model
dem_err.krig = krige(res_st_gt~1, herror, newdata=dem, vgf_herror_ste, nmax=20)
names(dem_err.krig)

spplot(dem_err.krig, zcol = "var1.pred", col.regions = bpy.colors())

dem_err.krig$var1.sd = sqrt(dem_err.krig$var1.var)
spplot(dem_err.krig, zcol = "var1.sd", col.regions = bpy.colors())
hist(dem_err.krig$var1.sd)

#DEM+error
newdem<-dem@data + dem_err.krig$var1.pred
coordinates(newdem) <- coordinates(dem)
gridded(newdem) = TRUE
spplot(newdem, col.regions = bpy.colors())
writeGDAL(drivername = "GTiff",newdem, "outputs/DEMdefault.tif")
summary(dem@data)
summary(newdem)

# simulation using ste model
nsims <- 1000
dem_err.sim = krige(res_st_gt~1, herror, newdata = dem, vgf_herror_ste, nsim = nsims, nmax = 20)
names(dem_err.sim)
#windows(width = 5, height = 6)
#spplot(dem_err.sim, zcol = c("sim1"), col.regions = bpy.colors())
#windows(width = 5, height = 6)
#spplot(demres.sim, col.regions = bpy.colors())

#summary (demres.sim[3])
#summary (dem_err.sim[1])



#save results
for(i in 1:nsims)
{
  outname <- paste("outputs/demsim_", i, ".tif", sep = "")
  outdem <- dem@data + dem_err.sim[i]@data
  # promote data.frame to SpatialPPixelsDataFrame
  coordinates(outdem) <- coordinates(dem)
  gridded(outdem) <- T
  writeGDAL(outdem, outname)
}

#load some examples

rm(list = ls())  # clean memory
