#------
#
# Extract Data from NASS By Different Commodities
#
# Generates RDS files for each commodity
#
# By: Ryan Miller
#
# Last updated: 16 April 2021
#
#-----

#---- Clean Up Work Space ----
rm(list = ls())
gc()
#---- END Clean Up ----

#---- Libraries ----
library(rnassqs)
library(operators)

source("C:/DATA/NASS/Code/NASS.Supporting.Functions.R")
#---- END Libraries ----

#---- Set Directories ----
setwd("C:/DATA/NASS/")
write.path <- "C:/Data/NASS/"
#---- END Directories ----

#---- Set System Parameters ----
#Request API Key here -> https://quickstats.nass.usda.gov/api
api_key <- "BC35C8D3-C647-39B5-A1C8-3D1164838CF4"
Sys.setenv(NASSQS_TOKEN = api_key)
#---- END System Parameters ----

#---- Setup Querry Domain ----

#--States to Process
state.vec <- rnassqs::nassqs_param_values(param = "state_alpha")
state.vec <- state.vec[state.vec %!in% c("OT", "GU", "MP", "AS", "US")]

#--Years to Process
year.vec <- c(1997:2020)

#--Sources
source.vec <- c("CENSUS", "SURVEY")

#---- END Querry Domain ----

#---- Query NASS Data ----

#--DEER and ELK
source.dat <- get.nass.data(
	source.vec,
	commodity.vec = c("ELK", "DEER"),
	year.vec,
	state.vec
)
saveRDS(
	source.dat,
	paste0(write.path, "NASS.Census.Elk.Deer.Data.", Sys.Date(), ".RDS")
)

#--HOGS
rm(source.dat)
source.dat <- get.nass.data(
	source.vec,
	commodity.vec = c("HOGS"),
	year.vec,
	state.vec
)
saveRDS(
	source.dat,
	paste0(write.path, "NASS.Census.Hogs.Data.", Sys.Date(), ".RDS")
)

#--CATTLE
rm(source.dat)
source.dat <- get.nass.data(
	source.vec,
	commodity.vec = c("CATTLE"),
	year.vec,
	state.vec
)
saveRDS(
	source.dat,
	paste0(write.path, "NASS.Census.Cattle.Data.", Sys.Date(), ".RDS")
)

#--SHEEP
rm(source.dat)
source.dat <- get.nass.data(
	source.vec,
	commodity.vec = c("SHEEP"),
	year.vec,
	state.vec
)
saveRDS(
	source.dat,
	paste0(write.path, "NASS.Census.Sheep.Data.", Sys.Date(), ".RDS")
)

#--GOATS
rm(source.dat)
source.dat <- get.nass.data(
	source.vec,
	commodity.vec = c("GOATS"),
	year.vec,
	state.vec
)
saveRDS(
	source.dat,
	paste0(write.path, "NASS.Census.Goats.Data.", Sys.Date(), ".RDS")
)

#---- END Extract NASS Data ----

#---- Impute missing / withheld data ----

#--Read County Adjacency File
adj.file <- read.csv("county_adjacency2010.csv", stringsAsFactors = FALSE)


#--Deer and Elk
dat <- readRDS(paste0(write.path, "NASS.Census.Elk.Deer.Data.2021-04-19.RDS"))

#Remove duplicates
dat <- unique(dat)

#Impute missing based on time series
imp.dat <- impute.missing(dat)

#Set Categories to Address
x <- unique(dat$short_desc)
short_desc.val <- x[!grepl("[$]", x)]

#Impute missing inventory
imp.dat <- impute.missing.inventory(
	imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Impute using spatial neighbors
sp.imp.dat <- spatial.impute.missing(
	imp.dat,
	adj.file = adj.file,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Fill in remaining missing using state mean
st.imp.dat <- impute.missing.using.state.mean(
	sp.imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Fill in remaining with National mean
nat.imp.dat <- impute.missing.using.national.mean(
	st.imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

saveRDS(
	nat.imp.dat,
	paste0(write.path, "NASS.Census.Elk.Deer.Data.Imputed.", Sys.Date(), ".RDS")
)


#--HOGS
dat <- readRDS(paste0(write.path, "NASS.Census.Hogs.Data.2021-04-19.RDS"))

#Remove duplicates
dat <- unique(dat)

#Impute missing based on time series
imp.dat <- impute.missing(dat)

#Set Categories to Address
x <- unique(dat$short_desc)
short_desc.val <- x[!grepl("[$]", x)]

#Impute missing inventory
imp.dat <- impute.missing.inventory(
	imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Impute using spatial neighbors
sp.imp.dat <- spatial.impute.missing(
	imp.dat,
	adj.file = adj.file,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Fill in remaining missing using state mean
st.imp.dat <- impute.missing.using.state.mean(
	sp.imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Fill in remaining with National mean
nat.imp.dat <- impute.missing.using.national.mean(
	st.imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

saveRDS(
	nat.imp.dat,
	paste0(write.path, "NASS.Census.Hogs.Imputed.", Sys.Date(), ".RDS")
)


#--CATTLE
dat <- readRDS(paste0(write.path, "NASS.Census.Cattle.Data.2021-04-19.RDS"))

#Remove duplicates
dat <- unique(dat)

#Impute missing based on time series
imp.dat <- impute.missing(dat)

#Set Categories to Address
x <- unique(dat$short_desc)
short_desc.val <- x[!grepl("[$]", x)]

#Impute missing inventory
imp.dat <- impute.missing.inventory(
	imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Impute using spatial neighbors
sp.imp.dat <- spatial.impute.missing(
	imp.dat,
	adj.file = adj.file,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Fill in remaining missing using state mean
st.imp.dat <- impute.missing.using.state.mean(
	sp.imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Fill in remaining with National mean
nat.imp.dat <- impute.missing.using.national.mean(
	st.imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

saveRDS(
	nat.imp.dat,
	paste0(write.path, "NASS.Census.Cattle.Data.Imputed.", Sys.Date(), ".RDS")
)


#--SHEEP
dat <- readRDS(paste0(write.path, "NASS.Census.Sheep.Data.2021-04-19.RDS"))

#Remove duplicates
dat <- unique(dat)

#Impute missing based on time series
imp.dat <- impute.missing(dat)

#Set Categories to Address
x <- unique(dat$short_desc)
short_desc.val <- x[!grepl("[$]", x)]

#Impute missing inventory
imp.dat <- impute.missing.inventory(
	imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Impute using spatial neighbors
sp.imp.dat <- spatial.impute.missing(
	imp.dat,
	adj.file = adj.file,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Fill in remaining missing using state mean
st.imp.dat <- impute.missing.using.state.mean(
	sp.imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Fill in remaining with National mean
nat.imp.dat <- impute.missing.using.national.mean(
	st.imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

saveRDS(
	nat.imp.dat,
	paste0(write.path, "NASS.Census.Sheep.Data.Imputed.", Sys.Date(), ".RDS")
)


#--GOATS
dat <- readRDS(paste0(write.path, "NASS.Census.Goats.Data.2021-04-19.RDS"))

#Remove duplicates
dat <- unique(dat)

#Impute missing based on time series
imp.dat <- impute.missing(dat)

#Set Categories to Address
x <- unique(dat$short_desc)
short_desc.val <- x[!grepl("[$]", x)]

#Impute missing inventory
imp.dat <- impute.missing.inventory(
	imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Impute using spatial neighbors
sp.imp.dat <- spatial.impute.missing(
	imp.dat,
	adj.file = adj.file,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Fill in remaining missing using state mean
st.imp.dat <- impute.missing.using.state.mean(
	sp.imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

#Fill in remaining with National mean
nat.imp.dat <- impute.missing.using.national.mean(
	st.imp.dat,
	statisticcat_desc.val = "INVENTORY",
	short_desc.val
)

saveRDS(
	nat.imp.dat,
	paste0(write.path, "NASS.Census.Goats.Data.Imputed.", Sys.Date(), ".RDS")
)

#---- END Impute Data ----

##---- END END ----
