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
data_path <- "Data"
# sp_data_path <- "C:/DATA/Cartography.Layers/" # not used
write_path <- "Data"


#----Load Packages for Session----

library(sp)
library(plyr)
library(raster)
library(tidyr)
library(anytime)
library(lubridate)
library(operators)
library(geosphere)

source("R/Supporting.Functions 2.R")

#----END Load libraries----

#---- Read Inspection Data ----
dat_tt <- read.csv(
  file.path(data_path, "CFT_InspectionDataInitial_8.9.20.csv"),
  stringsAsFactors = FALSE
)
dim(dat_tt)
names(dat_tt)
glimpse(dat_tt)

dat_tt$db.source <- "tt"

dat_scs <- read.csv(
  file.path(data_path, "scs.ontology.reshaped.csv"),
  stringsAsFactors = FALSE
)
dim(dat_scs)
names(dat_scs)
glimpse(dat_scs)

dat_scs$db.source <- "scs"

#Keep only those records with data

col_names <- c(
  "Qty_Herds",
  "Qty_Inspected",
  "Qty_Infested",
  "Qty_Moved",
  "Qty_Added",
  "Qty_Adults_Added",
  "Qty_Calves_Added"
)

dat_scs <- completeFun(dat_scs, col_names)
dim(dat_scs)
names(dat_scs)

#--Merge SCS and Tick Tracker Data
dat <- rbind.data.frame(dat_tt, dat_scs)
dim(dat)
names(dat)
glimpse(dat)

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
val_list <- c("", "None", "Horse River Patrol")
dat <- dat[dat$Inspection_Type %!in% val_list, ]

#--Consolidate inspection types
val_list <- c(
  "14-day Prem-VPatrol",
  "14-day Prem-Range",
  "14-day Prem-Pen",
  "14-day Pass"
)
dat[dat$Inspection_Type %in% val_list, "Inspection_Type"] <- "14-day inspection"

val_list <- c(
  "Equip Inspection",
  "Horse Patrol Inspection",
  "Vehicle Patrol Inspection"
)
dat[dat$Inspection_Type %in% val_list, "Inspection_Type"] <- "14-day inspection"

val_list <- c("Wildlife-Range")
dat[dat$Inspection_Type %in% val_list, "Inspection_Type"] <- "14-day inspection"

val_list <- c(
  "Premises",
  "Pen",
  "One-time Movement",
  "Issue Quarantine",
  "Infested"
)
dat[dat$Inspection_Type %in% val_list, "Inspection_Type"] <- "Scratch"

val_list <- c("Wildlife-Scratch")
dat[dat$Inspection_Type %in% val_list, "Inspection_Type"] <- "Scratch"

#--Consolidate species
val_list <- c("Nilgai Bull", "Nilgai Cow", "Nilgai Unknown")
dat[dat$Species %in% val_list, "Species"] <- "Nilgai"

val_list <- c(
  "Whitetail Buck",
  "Whitetail Deer",
  "Whitetail Doe",
  "Whitetail Unknown",
  "Deer"
)
dat[dat$Species %in% val_list, "Species"] <- "Whitetail"

val_list <- c("Other Species", "Other_Wildlife")
dat[dat$Species %in% val_list, "Species"] <- "Other_Wildlife"

val_list <- c("Premises")
dat[dat$Species %in% val_list, "Species"] <- "Bovine"

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
date_vals_tt <- as.Date(
  dat[dat$db.source == "tt", "Inspection_Date"],
  "%m/%d/%y"
)
date_vals_scs <- as.Date(
  dat[dat$db.source == "scs", "Inspection_Date"],
  "%Y-%m-%d"
)

dat[dat$db.source == "tt", "tmp"] <- date_vals_tt
dat[dat$db.source == "scs", "tmp"] <- date_vals_scs

dat$Inspection_Date <- dat$tmp
dat <- dat[, -ncol(dat)]

#Add year
dat$year <- year(dat$Inspection_Date)

#Add month
dat$month <- month(dat$Inspection_Date)

#---- END modify dates

#---- Make Site Level Data ----

#--Pasture Size
site_size <- aggregate(
  cbind(Pasture_Qty_Acres, Pasture_Latitude, Pasture_Longitude) ~ Pasture_Name +
    County_Name,
  data = dat,
  FUN = mean,
  na.rm = TRUE,
  na.action = na.pass
)

#--Set those sites with no size to NA

site_size[is.nan(site_size$Pasture_Qty_Acres), "Pasture_Qty_Acres"] <- 0
site_size[site_size$Pasture_Qty_Acres == 0, "Pasture_Qty_Acres"] <- NA

nrow(site_size)
length(unique(site_size$Pasture_Name))


#--Herd size

#Herd size by year
site_herd_size <- aggregate(
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
site_herd_size <- aggregate(
  cbind(Qty_Inspected) ~ Pasture_Latitude +
    Pasture_Longitude +
    Pasture_Name +
    County_Name +
    Species,
  data = site_herd_size,
  FUN = mean,
  na.rm = TRUE,
  na.action = na.pass
)

nrow(site_herd_size)
length(unique(site_herd_size$Pasture_Name))

site_dat <- merge(
  site_herd_size,
  site_size,
  by = c(
    "Pasture_Latitude",
    "Pasture_Longitude",
    "Pasture_Name",
    "County_Name"
  ),
  all.x = TRUE,
  all.y = TRUE
)

site_dat <- spread(site_dat, Species, Qty_Inspected)
nrow(site_dat)
length(unique(site_dat$Pasture_Name))

site_dat$site.wildlife <- rowSums(site_dat[, c(
  "Nilgai",
  "Whitetail",
  "Other_Wildlife"
)])

site_dat <- site_dat[,
  colnames(site_dat) %!in%
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
colnames(site_dat)[which(colnames(site_dat) == "Bovine")] <- "site.inv.bovine"
colnames(site_dat)[which(colnames(site_dat) == "Equine")] <- "site.inv.equine"
colnames(site_dat)[which(
  colnames(site_dat) == "Pasture_Qty_Acres"
)] <- "site.pasture_qty_acres"

#Set NA inventories to zero
site_dat[is.na(site_dat$site.inv.bovine) == TRUE, "site.inv.bovine"] <- 0
site_dat[is.na(site_dat$site.inv.equine) == TRUE, "site.inv.equine"] <- 0

#--Adjust units
site_dat$site.inv.bovine <- round(site_dat$site.inv.bovine)
site_dat$site.inv.equine <- round(site_dat$site.inv.equine)


#---- Generate Neighborhood Mean for site.pasture_qty_acres with missing area

site_dat <- neighbor.mean(site_dat, col.name = "site.pasture_qty_acres")

site_dat[
  is.nan(site_dat$site.pasture_qty_acres) == TRUE,
  "site.pasture_qty_acres"
] <- NA

#Fill in those remaining
ntmp <- length(which(is.na(site_dat$site.pasture_qty_acres)))
if (ntmp > 0) {
  site_dat <- neighbor.mean(site_dat, col.name = "site.pasture_qty_acres")
}

#---- END

#--Make Site Densities
site_dat$site.dens.bovine <- round(
  site_dat$site.inv.bovine / site_dat$site.pasture_qty_acres,
  digits = 3
)
site_dat$site.dens.equine <- round(
  site_dat$site.inv.equine / site_dat$site.pasture_qty_acres,
  digits = 3
)

#--Rename columns

#--Merge site and obs data
dat <- merge(
  dat,
  site_dat,
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
obs_dat <- aggregate(
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
nrow(obs_dat)
length(unique(obs_dat$Pasture_Name))

obs_dat_mu <- aggregate(
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
nrow(obs_dat_mu)
length(unique(obs_dat_mu$Pasture_Name))

obs_dat <- merge(
  obs_dat,
  obs_dat_mu,
  by = c("jdate", "year", "Pasture_Name", "Inspection_Type", "Species", "month")
)
nrow(obs_dat)
length(unique(obs_dat$Pasture_Name))

#--Convert infested to 0,1
obs_dat$y <- obs_dat$Qty_Infested
obs_dat[obs_dat$y > 0, "y"] <- 1

#site_dat$Qty_Added<-site_dat$Qty_Added+1

#---- Merge Observation and Site Level Data
out_dat <- merge(obs_dat, site_dat, by = c("Pasture_Name"))
nrow(out_dat)
length(unique(out_dat$Pasture_Name))

colnames(out_dat) <- tolower(colnames(out_dat))

# Formate for csvToUMF import function
#Requires data to be: Site ID, Date of Observation, Observations (y), Covariates....

#--Formate: Site level data, Covariates, Response....

col_names <- c(
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

wrt_dat <- out_dat[, col_names]

summary(wrt_dat)
glimpse(wrt_dat)

write.csv(
  wrt_dat,
  file.path(write_path, paste0("dat.occ.cov.", Sys.Date(), ".csv")),
  row.names = FALSE
)

#----END
