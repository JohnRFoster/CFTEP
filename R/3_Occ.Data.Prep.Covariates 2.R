#-------------------
#
# Add Site Level covariates
#
# By: Ryan Miller
#
#
# Last updated: 28 Oct 2020
#
#------------------

## Clean Workspace
rm(list = ls(all = TRUE))
gc()

options(scipen = 999)


#----Set Directories----
setwd("C:/Documents/Project Documents/CattleFeverTick/Data/")

code.path <- "C:/Documents/Project Documents/CattleFeverTick/Code/"
data.path <- "C:/Documents/Project Documents/CattleFeverTick/Model.data/"
raw.data.path <- "G:/CFTEP/Data/DatasetCleaning/"
sp.data.path <- "C:/Documents/Project Documents/CattleFeverTick/Data/SpatialData/"
write.path <- "C:/Documents/Project Documents/CattleFeverTick/Model.data/"


#----Load Packages for Session----
packages.vec <- c(
  "epitools",
  "rgdal",
  "sp",
  "plyr",
  "raster",
  "tidyr",
  "anytime",
  "lubridate",
  "operators",
  "geosphere",
  "deldir",
  "rgdal"
)

# Check if installed and install if not
if (
  FALSE %in% unique(is.element(packages.vec, rownames(installed.packages())))
) {
  install.packages(setdiff(packages.vec, rownames(installed.packages())))
}
# load
lapply(packages.vec, require, character.only = TRUE)

source(paste0(code.path, "Supporting.Functions.R"))

#----END Load libraries----

#---- Read Data ----
dat <- read.csv(
  paste0(data.path, "dat.occ.cov.2021-02-05.csv"),
  stringsAsFactors = FALSE
)


#---- Address issues with Lat Lon ----

#-- Exclude those with missing lat/lon data
dat <- dat[is.na(dat$pasture_latitude) == FALSE, ]
dat <- dat[is.na(dat$pasture_longitude) == FALSE, ]

#-- Remove Locations Outside Texas

#Texas Polygon
poly <- readOGR(
  dsn = file.path("C:/DATA/Cartography.Layers/"),
  layer = "Counties_dtl_geo"
)
poly <- poly[poly$STATE_NAME == "Texas", ]
colnames(poly@data) <- tolower(colnames(poly@data))

#Generate Points
pts <- dat

coordinates(pts) <- ~ pasture_longitude + pasture_latitude
crs(pts) <- "+proj=longlat +datum=WGS84 +no_defs"

#Intersect points and polygons
pts.dat <- raster::extract(x = poly, y = pts, method = "simple", df = TRUE)

dat <- cbind.data.frame(dat, pts.dat[, c("name", "state_name")])

#Remove observations with locations outside Texas
dat <- dat[is.na(dat$state_name) == FALSE, ]
nrow(dat)

#Remove locations not in appropriate county
dat <- dat[dat$county_name %in% dat$name, ]
nrow(dat)

#Drop processing columns
dat <- dat[, colnames(dat) %!in% c("name", "state_name")]

#---- END Address issues with lat / lon ----

#---- Derived Covariates ----

#--Add seasonal values based on Leon 2012 - doi.org/10.3389/fphys.2012.00195
dat[dat$month %in% c(1, 2, 3, 4, 10, 11, 12), "season"] <- 0
dat[dat$month %!in% c(1, 2, 3, 4, 10, 11, 12), "season"] <- 1

#-- END season values ----

#---- Generate Distance to Infested ----

#-- Median distance to infested pastures in each year
dat <- median.distance.to.infested(dat, type = "season")

#-- Convert distances to KM
dat$med.dist.pos <- dat$med.dist.pos / 1000
dat$sd.dist.pos <- dat$sd.dist.pos / 1000

#-- Make site level covariate
tmp <- aggregate(
  cbind(med.dist.pos, sd.dist.pos) ~ pasture_name,
  data = dat,
  FUN = mean
)
colnames(tmp) <- c("pasture_name", "site.med.dist.pos", "site.sd.dist.pos")

dat.dist <- merge(dat, tmp, by = "pasture_name", all.x = TRUE)
nrow(dat.dist)

dat <- dat.dist

#---- END Distance to Infested ----

#---- Generate Distance to Adjacent ----

#-- Distance to adjacent pastures in each year
tmp <- generate.adj.pasture.dist(dat)

tmp <- tmp[, c("pasture_name", "year", "mean.pas.dist", "sd.pas.dist")]

#-- Convert distances to KM
tmp$mean.pas.dist <- tmp$mean.pas.dist / 1000
tmp$sd.pas.dist <- tmp$sd.pas.dist / 1000

dat <- merge(dat, tmp, by = c("pasture_name", "year"), all.x = TRUE)

#-- Make site level covariate
tmp <- aggregate(
  cbind(mean.pas.dist, sd.pas.dist) ~ pasture_name,
  data = tmp,
  FUN = mean
)
colnames(tmp) <- c("pasture_name", "site.mean.pas.dist", "site.sd.pas.dist")

dat <- merge(dat, tmp, by = "pasture_name", all.x = TRUE)
nrow(dat)

#---- END Distance to Adjacent ----

#---- Generate infestation rate and count ----

lag.dat <- generate.lag.rate(dat)

lag.dat <- lag.dat[, c(
  "pasture_name",
  "year",
  "num.adj.pas",
  "lag.count",
  "lag.rate"
)]

dat <- merge(dat, lag.dat, by = c("pasture_name", "year"))

#Generate annual mean
tmp.agg <- aggregate(
  cbind(num.adj.pas, lag.count, lag.rate) ~ pasture_name,
  data = dat,
  FUN = mean
)
colnames(tmp.agg) <- c(
  "pasture_name",
  "site.num.adj.pas",
  "site.lag.count",
  "site.lag.rate"
)

dat <- merge(dat, tmp.agg, by = "pasture_name", all.x = TRUE)
nrow(dat)


#-- Generate infestation rate and count for previous year

#lag.dat<-generate.lag.1.rate(dat)

#lag.dat<-lag.dat[,c("pasture_name","year","lag.1.count","lag.1.rate")]

#dat<-merge(dat, lag.dat, by=c("pasture_name","year"))

#Generate annual mean
#tmp<-aggregate(cbind(num.adj.pas,lag.count,lag.rate)~pasture_name, data=dat, FUN=mean)
#colnames(tmp)<-c("pasture_name","site.num.adj.pas","site.lag.1.count","site.lag.1.rate")

#dat<-merge(dat,tmp, by="pasture_name", all.x=TRUE)
#nrow(dat)

#---- END Infestation Rate ----

#---- Add NASS Cattle Data ----

tmp.rds <- readRDS("C:/DATA/NASS/NASS.Census.Cattle.Data.2019-12-11.RDS")

# Limit to 2012 and 2017
tmp.rds <- tmp.rds[tmp.rds$year %in% c(2017, 2012), ]

# Limit to Texas
tmp.rds <- tmp.rds[tmp.rds$state_name == "TEXAS", ]

# Inventory
n.inv <- tmp.rds[tmp.rds$short_desc == "CATTLE, COWS - INVENTORY", ]
n.inv <- aggregate(Value.numeric ~ FIPS, data = n.inv, FUN = mean, rm.na = TRUE)
colnames(n.inv) <- c("FIPS", "inv")

# Number Farms
n.dat <- tmp.rds[
  tmp.rds$short_desc == "CATTLE, COWS - OPERATIONS WITH INVENTORY",
]
n.dat <- aggregate(Value.numeric ~ FIPS, data = n.dat, FUN = mean, rm.na = TRUE)
colnames(n.dat) <- c("FIPS", "frm")

# Merge Data
n.dat <- merge(n.dat, n.inv, by = "FIPS", all.x = TRUE)


# Load county boundries
sp.dat <- readOGR(
  dsn = file.path("C:/DATA/Cartography.Layers/"),
  layer = "Counties.Course.Lower48",
  stringsAsFactors = FALSE
)

sp.dat <- sp.dat[sp.dat$STATE_NAME == "Texas", ]
sp.dat <- sp.dat[, c("FIPS")]

# Merge spatial data and nass data
n.dat <- merge(sp.dat, n.dat, by = c("FIPS"))

#Generate Spatial Points for Intersect
pts <- dat

coordinates(pts) <- ~ pasture_longitude + pasture_latitude
crs(pts) <- "+proj=longlat +datum=WGS84 +no_defs"

#Intersect
pts.dat <- raster::extract(x = n.dat, y = pts, method = "simple", df = TRUE)

pts.dat <- cbind.data.frame(dat, pts.dat[, c("FIPS", "frm", "inv")])

#Add mean farm size
pts.dat$mean.farm.size <- pts.dat$inv / pts.dat$frm

#Return data
dat <- pts.dat

#---- END Add NASS data ----

#---- Assign suitability values to locations ----

#-- Generate raster stack
ras1 <- raster::raster(
  "C:/Documents/Project Documents/CattleFeverTick/Data/SpatialData/recode.giles.anulates.current.tif"
)
ras1[is.na(ras1) == TRUE] <- 0

ras2 <- raster::raster(
  "C:/Documents/Project Documents/CattleFeverTick/Data/SpatialData/recode.giles.microplus.current.tif"
)
ras2[is.na(ras2) == TRUE] <- 0

#Generate Points for Intersection
pts <- unique(dat[, c(
  "pasture_name",
  "pasture_longitude",
  "pasture_latitude",
  "site.pasture_qty_acres"
)])

#Area of pasture in acres
area <- pts$pasture_qty_acres * 4046.85642

#Use area to determine radius for extraction
rad <- sqrt(area / pi)


#-- Asign tick suitability
for (i in 1:nrow(pts)) {
  x1 <- raster::extract(
    x = ras1,
    y = pts[i, c("pasture_longitude", "pasture_latitude")],
    method = "simple",
    buffer = rad[i],
    na.rm = TRUE,
    df = TRUE,
    fun = mean
  )
  x2 <- raster::extract(
    x = ras2,
    y = pts[i, c("pasture_longitude", "pasture_latitude")],
    method = "simple",
    buffer = rad[i],
    na.rm = TRUE,
    df = TRUE,
    fun = mean
  )

  #Use max suitability among two tick species
  tmp <- cbind.data.frame(id = i, tick.suit = max(x1[, 2], x2[, 2]))

  if (i == 1) {
    out <- tmp
  }
  if (i > 1) {
    out <- rbind.data.frame(out, tmp)
  }
} #END Loop


#Merge results
tmp <- cbind.data.frame(
  pts[, c("pasture_name", "pasture_longitude", "pasture_latitude")],
  tick.suit = out$tick.suit
)

tmp <- merge(
  dat,
  tmp,
  by = c("pasture_name", "pasture_longitude", "pasture_latitude"),
  all.x = TRUE
)

dat <- tmp

#---- END Assign Tick Suitability ----

#---- Distance to Mexico-U.S. Border ----

#Read MX Border
sp.dat <- readOGR(dsn = file.path(sp.data.path), layer = "mx_border")

#Reduce to only unique pasture lat lon
pnt.dat <- unique(dat[, c(
  "pasture_name",
  "pasture_longitude",
  "pasture_latitude"
)])

#Calculate line distance (this takes awhile)
tmp <- dist2Line(
  p = as.matrix(pnt.dat[, c("pasture_longitude", "pasture_latitude")]),
  line = sp.dat
)
tmp <- as.data.frame(tmp)

#Convert to KM
tmp$distance <- tmp$distance / 1000

#Reshape and add pasture names
tmp$pasture_name <- pnt.dat$pasture_name
colnames(tmp) <- c(
  "distance",
  "pasture_longitude",
  "pasture_latitude",
  "id",
  "pasture_name"
)

tmp <- tmp[, c(
  "pasture_longitude",
  "pasture_latitude",
  "pasture_name",
  "distance"
)]
colnames(tmp)[ncol(tmp)] <- "mx.distance"

#Merge data
tmp.out <- tmp[, c("pasture_name", "mx.distance")]
tmp.out <- merge(dat, tmp.out, by = c("pasture_name"), all.x = TRUE)

dat <- tmp.out

#---- END MX Border Distance ----

#---- Add sample occasions ----

#--Order data
dat <- dat[order(dat$pasture_name, dat$jdate), ]

#--Make sampling occasions within each year/pasture
vec <- unique(dat$pasture_name)

year.vec <- unique(dat$year)

#--Loop Over Year and Pasture
for (y in 1:length(year.vec)) {
  pb <- txtProgressBar(min = 0, max = length(year.vec), style = 3)
  for (i in 1:length(vec)) {
    if (nrow(dat[dat$year == year.vec[y] & dat$pasture_name == vec[i], ]) > 0) {
      dat[
        dat$year == year.vec[y] & dat$pasture_name == vec[i],
        "sample"
      ] <- seq(
        1,
        nrow(dat[dat$year == year.vec[y] & dat$pasture_name == vec[i], ]),
        1
      )
    } #END Logical
  } #END Pasture
  setTxtProgressBar(pb, y)
} #END year

#--Add Sample Year
dat$sample.year <- paste0("y", dat$year, "s", dat$sample)

nrow(dat)
length(unique(dat$pasture_name))

#---- END Sample Occasions ----

#---- Write Data ----

write.csv(
  dat,
  paste0(write.path, "dat.occ.cov.sp.covariates.", Sys.Date(), ".csv"),
  row.names = FALSE
)

#----END
