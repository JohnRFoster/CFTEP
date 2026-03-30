#-------------------
#
# By: Ryan Miller
#
# Format data for Occupancy based modeling approaches.
#
# Last updated: 11 Sept 2020
#
#------------------

## Clean Workspace
rm(list = ls(all = TRUE))
gc()


#----Set Directories----
code.path <- "C:/Documents/Project Documents/CattleFeverTick/Code/"
data.path <- "C:/Documents/Project Documents/CattleFeverTick/Data/"
sp.data.path <- "C:/DATA/Cartography.Layers/"
write.path <- "C:/Documents/Project Documents/CattleFeverTick/Model.data/"


#----Load Packages for Session----
packages.vec <- c(
  "sp",
  "plyr",
  "raster",
  "tidyr",
  "anytime",
  "lubridate",
  "operators",
  "geosphere"
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

#---- Read Inspection Data ----
dat.tt <- read.csv(
  paste0(data.path, "CFT_InspectionDataInitial_8.9.20.csv"),
  stringsAsFactors = FALSE
)
dim(dat.tt)
names(dat.tt)

dat.tt$db.source <- "tt"


dat.scs <- read.csv(
  paste0(data.path, "scs.ontology.reshaped.2021-01-30.csv"),
  stringsAsFactors = FALSE
)
dim(dat.scs)
names(dat.scs)

dat.scs$db.source <- "scs"

#Keep only those records with data

col.names <- c(
  "Qty_Herds",
  "Qty_Inspected",
  "Qty_Infested",
  "Qty_Moved",
  "Qty_Added",
  "Qty_Adults_Added",
  "Qty_Calves_Added"
)

dat.scs <- completeFun(dat.scs, col.names)
dim(dat.scs)
names(dat.scs)

#--Merge SCS and Tick Tracker Data
dat <- rbind.data.frame(dat.tt, dat.scs)
dim(dat)
names(dat)

length(unique(dat$Pasture_Name))

#---- Fix Lat Lon Data ----

#--Assume missing Pasture lat/lon = Prem lat / lon
dat[is.na(dat$Pasture_Longitude) == TRUE, "Pasture_Longitude"] <- dat[
  is.na(dat$Pasture_Longitude) == TRUE,
  "Premises_Longitude"
]
dat[is.na(dat$Pasture_Latitude) == TRUE, "Pasture_Latitude"] <- dat[
  is.na(dat$Pasture_Latitude) == TRUE,
  "Premises_Latitude"
]


#--Address multple lat / lon for same pasture

tmp <- aggregate(
  Pasture_Longitude ~ Pasture_Name,
  data = dat,
  FUN = mean,
  na.rm = TRUE,
  na.action = na.pass
)
colnames(tmp)[ncol(tmp)] <- "value"
dat <- merge(dat, tmp, by = "Pasture_Name", all.x = TRUE)
dat$Pasture_Longitude <- dat$value
dat <- dat[-ncol(dat)]

tmp <- aggregate(
  Pasture_Latitude ~ Pasture_Name,
  data = dat,
  FUN = mean,
  na.rm = TRUE,
  na.action = na.pass
)
colnames(tmp)[ncol(tmp)] <- "value"
dat <- merge(dat, tmp, by = "Pasture_Name", all.x = TRUE)
dat$Pasture_Latitude <- dat$value
dat <- dat[-ncol(dat)]

dat <- dat[is.na(dat$Pasture_Longitude) == FALSE, ]
nrow(dat)

dat <- dat[is.na(dat$Pasture_Latitude) == FALSE, ]
nrow(dat)


#Fix Longitude value sign
dat[dat$Pasture_Longitude > 0, "Pasture_Longitude"] <- dat[
  dat$Pasture_Longitude > 0,
  "Pasture_Longitude"
] *
  -1

dat[abs(dat$Pasture_Longitude) < 1, "Pasture_Longitude"] <- dat[
  abs(dat$Pasture_Longitude) < 1,
  "Pasture_Longitude"
] *
  100


#---- END Lat Lon Data

#---- Fix Missing Areas ----

#--Assume missing Pasture area = Premises area
dat[
  dat$Pasture_Qty_Acres == 0 & is.na(dat$Premises_Qty_Acres) == FALSE,
  "Pasture_Qty_Acres"
] <- dat[
  dat$Pasture_Qty_Acres == 0 & is.na(dat$Premises_Qty_Acres) == FALSE,
  "Premises_Qty_Acres"
]

#---- End Missing Areas

#---- Assign Vacated Pasture ----

#--Set pasture vacated to yes based on inspection type
dat[dat$Inspection_Type %in% c("Vacated premises"), "Pasture_Vacated"] <- "Y"

#--Assume NA values in mean the pasture is not vacated
dat[is.na(dat$Pasture_Vacated) == TRUE, "Pasture_Vacated"] <- "N"

#---- END

#---- Group Some Categories ----

plyr::count(dat$Inspection_Type)

#--Remove empty and 'none' inspection types
val.list <- c("", "None", "Horse River Patrol")
dat <- dat[dat$Inspection_Type %!in% val.list, ]

#--Consolidate inspection types
val.list <- c(
  "14-day Prem-VPatrol",
  "14-day Prem-Range",
  "14-day Prem-Pen",
  "14-day Pass"
)
dat[dat$Inspection_Type %in% val.list, "Inspection_Type"] <- "14-day inspection"

val.list <- c(
  "Equip Inspection",
  "Horse Patrol Inspection",
  "Vehicle Patrol Inspection"
)
dat[dat$Inspection_Type %in% val.list, "Inspection_Type"] <- "14-day inspection"

val.list <- c("Wildlife-Range")
dat[dat$Inspection_Type %in% val.list, "Inspection_Type"] <- "14-day inspection"

val.list <- c(
  "Premises",
  "Pen",
  "One-time Movement",
  "Issue Quarantine",
  "Infested"
)
dat[dat$Inspection_Type %in% val.list, "Inspection_Type"] <- "Scratch"

val.list <- c("Wildlife-Scratch")
dat[dat$Inspection_Type %in% val.list, "Inspection_Type"] <- "Scratch"

#--Consolidate species
val.list <- c("Nilgai Bull", "Nilgai Cow", "Nilgai Unknown")
dat[dat$Species %in% val.list, "Species"] <- "Nilgai"

val.list <- c(
  "Whitetail Buck",
  "Whitetail Deer",
  "Whitetail Doe",
  "Whitetail Unknown",
  "Deer"
)
dat[dat$Species %in% val.list, "Species"] <- "Whitetail"

val.list <- c("Other Species", "Other_Wildlife")
dat[dat$Species %in% val.list, "Species"] <- "Other_Wildlife"

val.list <- c("Premises")
dat[dat$Species %in% val.list, "Species"] <- "Bovine"

dat[is.na(dat$Species) == TRUE, "Species"] <- "None"

#---- End Regrouping

#---- Drop Some Data Not Useful for Occupancy ----

plyr::count(dat$Inspection_Type)
plyr::count(dat$Species)
plyr::count(dat$Pasture_Qty_Acres)
plyr::count(dat$db.source)
plyr::count(dat$County_Name)

plyr::count(dat[dat$Qty_Inspected == 0, c("Species")])


nrow(dat)

#--Remove those with 0 inspected animals
dat <- dat[dat$Qty_Inspected != 0, ]
dat <- dat[is.na(dat$Qty_Inspected) == FALSE, ]
nrow(dat)

#--Remove those with 0 inspected animals
dat <- dat[dat$Species != "None", ]
dat <- dat[is.na(dat$Species) == FALSE, ]
nrow(dat)


#---- Modify Dates ----

#Convert to dates
date.vals.tt <- as.Date(
  dat[dat$db.source == "tt", "Inspection_Date"],
  "%m/%d/%y"
)
date.vals.scs <- as.Date(
  dat[dat$db.source == "scs", "Inspection_Date"],
  "%Y-%m-%d"
)

dat[dat$db.source == "tt", "tmp"] <- date.vals.tt
dat[dat$db.source == "scs", "tmp"] <- date.vals.scs

dat$Inspection_Date <- dat$tmp
dat <- dat[, -ncol(dat)]

#Add year
dat$year <- year(dat$Inspection_Date)

#Add month
dat$month <- month(dat$Inspection_Date)

#---- END modify dates

#---- Make Site Level Data ----

#--Pasture Size
site.size <- aggregate(
  cbind(Pasture_Qty_Acres, Pasture_Latitude, Pasture_Longitude) ~ Pasture_Name +
    County_Name,
  data = dat,
  FUN = mean,
  na.rm = TRUE,
  na.action = na.pass
)

#--Set those sites with no size to NA

site.size[is.nan(site.size$Pasture_Qty_Acres), "Pasture_Qty_Acres"] <- 0
site.size[site.size$Pasture_Qty_Acres == 0, "Pasture_Qty_Acres"] <- NA

nrow(site.size)
length(unique(site.size$Pasture_Name))


#--Herd size

#Herd size by year
site.herd.size <- aggregate(
  cbind(Qty_Inspected) ~ Pasture_Latitude +
    Pasture_Longitude +
    Pasture_Name +
    County_Name +
    year +
    Species,
  data = dat,
  FUN = mean,
  na.rm = TRUE,
  na.action = na.pass
)

#Annual herd size
site.herd.size <- aggregate(
  cbind(Qty_Inspected) ~ Pasture_Latitude +
    Pasture_Longitude +
    Pasture_Name +
    County_Name +
    Species,
  data = site.herd.size,
  FUN = mean,
  na.rm = TRUE,
  na.action = na.pass
)

nrow(site.herd.size)
length(unique(site.herd.size$Pasture_Name))

site.dat <- merge(
  site.herd.size,
  site.size,
  by = c(
    "Pasture_Latitude",
    "Pasture_Longitude",
    "Pasture_Name",
    "County_Name"
  ),
  all.x = TRUE,
  all.y = TRUE
)

site.dat <- spread(site.dat, Species, Qty_Inspected)
nrow(site.dat)
length(unique(site.dat$Pasture_Name))

site.dat$site.wildlife <- rowSums(site.dat[, c(
  "Nilgai",
  "Whitetail",
  "Other_Wildlife"
)])

site.dat <- site.dat[,
  colnames(site.dat) %!in%
    c(
      "Nilgai",
      "None",
      "Other Species",
      "Whitetail",
      "V1",
      "",
      "<NA>",
      "Deer",
      "Other_Wildlife"
    )
]


#Change colnames
colnames(site.dat)[which(colnames(site.dat) == "Bovine")] <- "site.inv.bovine"
colnames(site.dat)[which(colnames(site.dat) == "Equine")] <- "site.inv.equine"
colnames(site.dat)[which(
  colnames(site.dat) == "Pasture_Qty_Acres"
)] <- "site.pasture_qty_acres"

#Set NA inventories to zero
site.dat[is.na(site.dat$site.inv.bovine) == TRUE, "site.inv.bovine"] <- 0
site.dat[is.na(site.dat$site.inv.equine) == TRUE, "site.inv.equine"] <- 0

#--Adjust units
site.dat$site.inv.bovine <- round(site.dat$site.inv.bovine)
site.dat$site.inv.equine <- round(site.dat$site.inv.equine)


#---- Generate Neighborhood Mean for site.pasture_qty_acres with missing area

site.dat <- neighbor.mean(site.dat, col.name = "site.pasture_qty_acres")

site.dat[
  is.nan(site.dat$site.pasture_qty_acres) == TRUE,
  "site.pasture_qty_acres"
] <- NA

#Fill in those remaining
site.dat <- neighbor.mean(site.dat, col.name = "site.pasture_qty_acres")

#---- END

#--Make Site Densities
site.dat$site.dens.bovine <- round(
  site.dat$site.inv.bovine / site.dat$site.pasture_qty_acres,
  digits = 3
)
site.dat$site.dens.equine <- round(
  site.dat$site.inv.equine / site.dat$site.pasture_qty_acres,
  digits = 3
)

#--Rename columns

#--Merge site and obs data
dat <- merge(
  dat,
  site.dat,
  by = c("Pasture_Name", "Pasture_Latitude", "Pasture_Longitude", "County_Name")
)

nrow(dat)
length(unique(dat$Pasture_Name))

#---- END Site Level

#---- Derived variables ----

#--Add apparent prevalence of infestation
dat$apr.infest.prev <- dat$Qty_Infested / dat$Qty_Inspected

#--Add proportion new animals
dat$prop.new.animals <- dat$Qty_Added / dat$Qty_Inspected

#---- END Derived variables

#---- For Occupancy Models ----

#--Aggregate occasion needs to be a julian date
dat$jdate <- format(dat$Inspection_Date, "%j")

#---- Generate Observation Level Covariates

#--Aggregate data for same visit
obs.dat <- aggregate(
  cbind(Qty_Infested, Qty_Inspected, Qty_Added) ~ jdate +
    Pasture_Name +
    Inspection_Type +
    Species +
    month +
    year,
  data = dat,
  FUN = sum,
  na.rm = TRUE,
  na.action = na.pass
)
nrow(obs.dat)
length(unique(obs.dat$Pasture_Name))

obs.dat.mu <- aggregate(
  cbind(apr.infest.prev, prop.new.animals) ~ jdate +
    year +
    Pasture_Name +
    Inspection_Type +
    Species +
    month,
  data = dat,
  FUN = mean,
  na.rm = TRUE,
  na.action = na.pass
)
nrow(obs.dat.mu)
length(unique(obs.dat.mu$Pasture_Name))

obs.dat <- merge(
  obs.dat,
  obs.dat.mu,
  by = c("jdate", "year", "Pasture_Name", "Inspection_Type", "Species", "month")
)
nrow(obs.dat)
length(unique(obs.dat$Pasture_Name))

#--Convert infested to 0,1
obs.dat$y <- obs.dat$Qty_Infested
obs.dat[obs.dat$y > 0, "y"] <- 1

#site.dat$Qty_Added<-site.dat$Qty_Added+1

#---- Merge Observation and Site Level Data
out.dat <- merge(obs.dat, site.dat, by = c("Pasture_Name"))
nrow(out.dat)
length(unique(out.dat$Pasture_Name))

colnames(out.dat) <- tolower(colnames(out.dat))

# Formate for csvToUMF import function
#Requires data to be: Site ID, Date of Observation, Observations (y), Covariates....

#--Formate: Site level data, Covariates, Response....

col.names <- c(
  "pasture_name",
  "jdate",
  "y",
  "pasture_latitude",
  "pasture_longitude",
  "site.pasture_qty_acres",
  "county_name",
  "year",
  "month",
  "species",
  "qty_infested",
  "qty_inspected",
  "inspection_type",
  "qty_added",
  "apr.infest.prev",
  "prop.new.animals",
  "site.inv.bovine",
  "site.inv.equine",
  "site.dens.bovine",
  "site.dens.equine"
)

wrt.dat <- out.dat[, col.names]


summary(wrt.dat)

write.csv(
  wrt.dat,
  paste0(write.path, "dat.occ.cov.", Sys.Date(), ".csv"),
  row.names = FALSE
)

#----END
