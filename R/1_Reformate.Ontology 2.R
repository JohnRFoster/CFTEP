#---------------
#
# Reshape Ontology Data
#
# By: Ryan Miller, John Foster
#
#---------------

#---- Set Paths ----
data_path <- "Data"
code_path <- "R"


#---- Load Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(operators)

# TODO: rename without spaces
source(file.path(code_path, "Supporting.Functions 2.R"))

#---- Inspection Data ----
# inspection_data <- "CFTEP+7-22+FDE+Survey+Ontology.csv" old
inspection_data <- "CFTEP 7-22 FDE Survey Ontology_19.11.24.csv"
dat <- read_csv(
  file.path(data_path, inspection_data)
)
dim(dat)
names(dat)


#-- Remove NA and "" from master columns
#-- Drop problematic columns
col_vec <- c(
  "scs_state_survey_id",
  "title",
  "Remarks",
  # "comments", not in newest data
  # "global_premises_id", not in newest data
  # "prem_location", not in newest data
  # "owner_name", not in newest data
  # "owner_address", not in newest data
  "inspector_name",
  # "Master_Restriction_Name", not in newest data
  # "national_identifiers", not in newest data
  "Quarantine_End_Date",
  "Quarantine_Start_Date"
)

dat_tidy <- dat |>
  filter(
    # other_identifiers != "", missing, in old code
    !is.na(latitude),
    !is.na(longitude),
    !is.na(inspection_type),
    !is.na(completion_date)
  ) |>
  select(-all_of(col_vec)) |>
  rename(
    # Pasture_Name = other_identifiers, not in newest data
    Pasture_Latitude = latitude,
    Pasture_Longitude = longitude,
    County_Name = county,
    Inspection_Date = completion_date,
    Inspection_Type = inspection_type,
    Pasture_Status = Quarantine_Reason
  )

#--Remove problamatic duplicates
cnt <- plyr::count(dat_tidy)
dat_cnt <- cnt[cnt$freq == 1, ]
dat_cnt <- dat_cnt[, -ncol(dat_cnt)]
nrow(dat_cnt)

#-- Rename columns

pattern_vec <- c(
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

for (i in seq_along(pattern_vec)) {
  col_names <- alter.column.names(
    in.dat = dat_cnt,
    pattern.str = pattern_vec[i]
  )

  colnames(dat_cnt) <- col_names
} #END Loop


#Reorder Columns
new_col <- c(
  "Pasture_Latitude",
  "Pasture_Longitude",
  "County_Name",
  "Inspection_Date",
  "Inspection_Type",
  "Pasture_Status"
)

col_names <- colnames(dat_cnt)
col_names <- col_names[order(col_names)]
col_names <- c(new_col, col_names[col_names %!in% new_col])
dat_cnt2 <- dat_cnt[, c(new_col, col_names[col_names %!in% new_col])]

#Look for duplicates (needed for reshaping below)
tmp <- plyr::count(dat_cnt2[, new_col])
tmp <- tmp[order(-tmp$freq), ]
nrow(dat_cnt2)

#--Reshape Data - converting wide columns to long
pattern_vec <- c(
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
for (i in seq_along(pattern_vec)) {
  tmp <- wide.to.long(
    in.dat = dat_cnt2,
    mast.col = new_col,
    pattern.str = pattern_vec[i]
  )

  col_name <- colnames(tmp)

  if (i == 1) {
    out <- tmp
  }
  if (i > 1) {
    out <- cbind.data.frame(out, tmp[, ncol(tmp)])
    colnames(out)[ncol(out)] <- col_name[length(col_name)]
  }
  #if(i>1){out<-merge(out,tmp,by=new_col)}

  print(nrow(out))
} #END Loop

nrow(out)


#--Add site level variables

#Acreage
x <- dat_cnt2 |>
  filter(!is.na(acreage)) |>
  group_by(
    # Pasture_Name,
    Pasture_Latitude,
    Pasture_Longitude,
    County_Name,
    Inspection_Date
  ) |>
  reframe(
    acreage = mean(acreage)
  )

tmp <- left_join(
  out,
  x,
  by = c(
    # "Pasture_Name",
    "Pasture_Latitude",
    "Pasture_Longitude",
    "County_Name",
    "Inspection_Date"
  )
)
# no fips in new data
# tmp$County_Code <- substring(tmp$STCOUNTYFP, first = 3, last = 5)

dat <- tmp

#--Rename Columns for easy merging
dat_qty <- dat |>
  rename(
    Qty_Inspected = Inspected,
    Qty_Infested = Infested,
    Qty_Added = Added,
    Qty_Moved = Moved,
    Qty_Dead = Dead,
    Qty_Herds = Herds,
    Qty_Dipped = Dipped,
    Qty_Frozen = Frozen,
    Qty_Dectomax = Dectomax,
    Qty_Sprayed = Sprayed,
    Qty_BM86 = BM86,
    Pasture_Qty_Acres = acreage
  )

#--Read tick tracker data and add missing columns

#Read csv file
template_dat <- read.csv(
  file.path(data_path, "CFT_InspectionDataInitial_8.9.20.csv"),
  stringsAsFactors = FALSE
)
template_dat <- template_dat[FALSE, ]

#Get column names
template_cols <- colnames(template_dat)
col_names <- template_cols[template_cols %!in% colnames(dat_qty)]

#Generate template data
template_dat <- template_dat[, colnames(template_dat) %in% col_names]

#Add NA values
template_dat[seq_len(nrow(dat_qty)), ] <- NA

#Add to SCS data
dat <- bind_cols(dat_qty, template_dat)

#Reorder to match tick tracker
dat <- dat[, template_cols]


#--Alter numeric columns and deal with NA values
col_names <- c(
  "Qty_Herds",
  "Qty_Inspected",
  "Qty_Infested",
  "Qty_Moved",
  "Qty_Added",
  "Qty_Adults_Added",
  "Qty_Calves_Added"
)

#Convert to numeric
for (i in seq_along(col_names)) {
  dat[, col_names[i]] <- as.numeric(dat[, col_names[i]])
} #END Loop

data <- dat |>
  mutate(across(starts_with("Qty"), as.numeric)) |>
  glimpse()

#Convert NA to -9 and 0 to -9
for (i in seq_along(col_names)) {
  x <- dat[, col_names[i]]
  dat[is.na(x) == TRUE, col_names[i]] <- -9

  x <- dat[, col_names[i]]
  dat[(x == 0) == TRUE, col_names[i]] <- -9
  print(summary(dat[, col_names[i]]))
} #END Loop

#Convert all -9 to NA
for (i in seq_along(col_names)) {
  x <- dat[, col_names[i]]
  dat[(x == -9) == TRUE, col_names[i]] <- NA
  print(summary(dat[, col_names[i]]))
} #END Loop

#--Write Data
write_csv(
  dat,
  file.path(data_path, "scs.ontology.reshaped.csv")
)

##---- END END ----
