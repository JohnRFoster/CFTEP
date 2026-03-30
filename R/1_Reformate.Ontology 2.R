#---------------
#
# Reshape Ontology Data
#
# By: Ryan Miller
#
#---------------

#---- Set Paths ----
data.path <- "C:/Documents/Project Documents/CattleFeverTick/Data/"
code.path <- "C:/Documents/Project Documents/CattleFeverTick/Code/"


#---- Load Libraries ----
library(tidyr)
library(stringr)
library(operators)

source(paste0(code.path, "Supporting.Functions.R"))


#---- Inspection Data ----
dat <- read.csv(
  paste0(data.path, "CFTEP+7-22+FDE+Survey+Ontology.csv"),
  stringsAsFactors = FALSE
)
dim(dat)
names(dat)


#-- Remove NA and "" from master columns

dat <- dat[dat$other_identifiers != "", ]
nrow(dat)

dat <- dat[dat$completion_date != "", ]
nrow(dat)

dat <- dat[is.na(dat$latitude) == FALSE, ]
nrow(dat)

dat <- dat[is.na(dat$longitude) == FALSE, ]
nrow(dat)

dat <- dat[is.na(dat$inspection_type) == FALSE, ]
nrow(dat)

dat <- dat[is.na(dat$completion_date) == FALSE, ]
nrow(dat)


#-- Drop problematic columns

col.vec <- c(
  "scs_state_survey_id",
  "title",
  "Remarks",
  "comments",
  "global_premises_id",
  "prem_location",
  "prem_location",
  "owner_name",
  "owner_address",
  "inspector_name",
  "Master_Restriction_Name",
  "national_identifiers",
  "Quarantine_End_Date",
  "Quarantine_Start_Date"
)

dat <- dat[, colnames(dat) %!in% col.vec]


org.col <- c(
  "other_identifiers",
  "latitude",
  "longitude",
  "county",
  "completion_date",
  "inspection_type",
  "Quarantine_Reason"
)

new.col <- c(
  "Pasture_Name",
  "Pasture_Latitude",
  "Pasture_Longitude",
  "County_Name",
  "Inspection_Date",
  "Inspection_Type",
  "Pasture_Status"
)

dat <- rename.columns(dat, org.col, new.col)


#--Remove problamatic duplicates
cnt <- plyr::count(dat)
dat <- cnt[cnt$freq == 1, ]
dat <- dat[, -ncol(dat)]
nrow(dat)


#-- Rename columns

pattern.vec <- c(
  "Inspected",
  "Infested",
  "Added",
  "Moved",
  "Dead",
  "Herds",
  "Dipped",
  "Frozen",
  "Dectomax",
  "Sprayed",
  "BM86"
)


for (i in 1:length(pattern.vec)) {
  col.names <- alter.column.names(in.dat = dat, pattern.str = pattern.vec[i])

  colnames(dat) <- col.names
} #END Loop


#Reorder Columns

col.names <- colnames(dat)

col.names <- col.names[order(col.names)]

col.names <- c(new.col, col.names[col.names %!in% new.col])


dat <- dat[, c(new.col, col.names[col.names %!in% new.col])]

#Look for duplicates (needed for reshaping below)
tmp <- plyr::count(dat[, new.col])
tmp <- tmp[order(-tmp$freq), ]
nrow(dat)


#--Reshape Data - converting wide columns to long

pattern.vec <- c(
  "Inspected",
  "Infested",
  "Added",
  "Moved",
  "Dead",
  "Herds",
  "Dipped",
  "Frozen",
  "Dectomax",
  "Sprayed",
  "BM86"
)

#Loop over
for (i in 1:length(pattern.vec)) {
  tmp <- wide.to.long(
    in.dat = dat,
    mast.col = new.col,
    pattern.str = pattern.vec[i]
  )

  col.name <- colnames(tmp)

  if (i == 1) {
    out <- tmp
  }
  if (i > 1) {
    out <- cbind.data.frame(out, tmp[, ncol(tmp)])
    colnames(out)[ncol(out)] <- col.name[length(col.name)]
  }
  #if(i>1){out<-merge(out,tmp,by=new.col)}

  print(nrow(out))
} #END Loop

nrow(out)


#--Add site level variables

#Acreage
x <- aggregate(
  acreage ~ Pasture_Name +
    Pasture_Latitude +
    Pasture_Longitude +
    County_Name +
    Inspection_Date +
    STCOUNTYFP,
  data = dat,
  FUN = mean
)

tmp <- merge(
  out,
  x,
  by = c(
    "Pasture_Name",
    "Pasture_Latitude",
    "Pasture_Longitude",
    "County_Name",
    "Inspection_Date"
  ),
  all.x = TRUE
)

tmp$County_Code <- substring(tmp$STCOUNTYFP, first = 3, last = 5)

dat <- tmp

#--Rename Columns for easy merging

org.col <- c(
  "Inspected",
  "Infested",
  "Added",
  "Moved",
  "Dead",
  "Herds",
  "Dipped",
  "Frozen",
  "Dectomax",
  "Sprayed",
  "BM86"
)
new.col <- c(
  "Inspected",
  "Infested",
  "Added",
  "Moved",
  "Dead",
  "Herds",
  "Dipped",
  "Frozen",
  "Dectomax",
  "Sprayed",
  "BM86"
)

dat <- rename.columns(dat, org.col = org.col, new.col = paste0("Qty_", new.col))


org.col <- c("acreage")
new.col <- c("Pasture_Qty_Acres")

dat <- rename.columns(dat, org.col = org.col, new.col = paste0("Qty_", new.col))


#--Read tick tracker data and add missing columns

#Read csv file
template.dat <- read.csv(
  paste0(data.path, "CFT_InspectionDataInitial_8.9.20.csv"),
  stringsAsFactors = FALSE
)
template.dat <- template.dat[FALSE, ]

#Get column names
template.cols <- colnames(template.dat)
col.names <- template.cols[template.cols %!in% colnames(dat)]

#Generate template data
template.dat <- template.dat[, colnames(template.dat) %in% col.names]

#Add NA values
template.dat[1:nrow(dat), ] <- NA

#Add to SCS data
dat <- cbind.data.frame(dat, template.dat)

#Reorder to match tick tracker
dat <- dat[, template.cols]


#--Alter numeric columns and deal with NA values
col.names <- c(
  "Qty_Herds",
  "Qty_Inspected",
  "Qty_Infested",
  "Qty_Moved",
  "Qty_Added",
  "Qty_Adults_Added",
  "Qty_Calves_Added"
)

#Convert to numeric
for (i in 1:length(col.names)) {
  dat[, col.names[i]] <- as.numeric(dat[, col.names[i]])
} #END Loop

#Convert NA to -9 and 0 to -9
for (i in 1:length(col.names)) {
  x <- dat[, col.names[i]]
  dat[is.na(x) == TRUE, col.names[i]] <- -9

  x <- dat[, col.names[i]]
  dat[(x == 0) == TRUE, col.names[i]] <- -9
  print(summary(dat[, col.names[i]]))
} #END Loop

#Convert all -9 to NA
for (i in 1:length(col.names)) {
  x <- dat[, col.names[i]]
  dat[(x == -9) == TRUE, col.names[i]] <- NA
  print(summary(dat[, col.names[i]]))
} #END Loop


#--Write Data
write.csv(
  dat,
  paste0(data.path, "scs.ontology.reshaped.", Sys.Date(), ".csv"),
  row.names = FALSE
)

##---- END END ----
