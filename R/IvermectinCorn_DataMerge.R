#Merging ivermectin corn data
#December 1st, 2021

library(tidyverse)
library(runner)
library(lubridate)
ivermectindata <- read.csv("CFT_Ivermectin_Cleaned_21.9.21.csv")
## amending date variables
ivermectindata <- ivermectindata %>%
  mutate(across(
    starts_with("Date"),
    ~ replace_na(as.Date(., format = c("%d-%b-%y")), today())
  ))

inspectiondata <- read.csv("IvermectinCornCountyInspections.csv")
inspectiondata <- inspectiondata %>%
  mutate(
    Inspection_Date = replace_na(
      as.Date(Inspection_Date, format = c("%m/%d/%Y")),
      today()
    )
  )


names(ivermectindata)
# [1] "Premises_No_Original"     "Premises_No"              "Date_Visited"             "Station_No"
# [5] "Corn_Added"               "Corn_Lot_No"              "Roller_Insecticide_Added" "Inspector"
# [9] "County"                   "Area_Type"                "Quarantine_Type"          "Owner_L_Name"
# [13] "Owner_F_Name"             "Ranch_Name"               "Pasture_Name"             "Pasture_Size"
# [17] "Gate_GPS_Lat"             "Gate_GPS_Lon"             "Treatment"                "Target_Species"
# [21] "Population"               "Other_Species"            "Station_Lat"              "Station_Lon"
# [25] "Feeder_Type"              "Date_in_Use"              "Date_Retired"

names(inspectiondata)
# [1] "Pasture_Name"      "Premises_No"       "Pasture_Latitude"  "Pasture_Longitude" "Pasture_Qty_Acres" "County_Name"
# [7] "County_Code"       "Inspection_Date"   "Inspection_Type"   "Pasture_Status"    "Species"           "Qty_Herds"
# [13] "Qty_Inspected"     "Qty_Infested"      "Qty_Moved"         "Qty_Added"         "Qty_Dipped"        "Qty_Frozen"
# [19] "Qty_Sprayed"       "Qty_Dectomax"      "Qty_BM86"          "Qty_Dead"          "Source"

#OBJECTIVE: add new variables to inspectiondata, based on ivermectindata#
#Premises_No is the matching variable#

#Key variable: days between ivermectindata(Date_Visited) and inspectiondata(Inspection_Date)

#QUESTION 1: How much TOTAL corn has been added add premise since last survey
#
#Variable: Qty_Corn

## not all premise_No from the inspection data exist in the ivermectin data
premise_sub <- unique(inspectiondata$Premises_No[
  which(inspectiondata$Premises_No %in% ivermectindata$Premises_No)
])

corn_feed <- list()
for (i in 1:length(premise_sub)) {
  ## extract inspections and dates
  foo <- inspectiondata %>%
    filter(Premises_No == premise_sub[[i]]) %>%
    group_by(Premises_No) %>%
    arrange(Inspection_Date, .by_group = TRUE) %>%
    summarise(Inspection_Date = unique(Inspection_Date)) %>%
    mutate(
      Qty_Corn = 0,
      Qty_TotalCornDays = 0,
      Qty_CornFeeders = 0,
      Feeder_Density = 0
    ) ## initialize dummy columns

  if (nrow(foo) > 1) {
    ## move along number of rows in foo
    for (j in 2:nrow(foo)) {
      ## subset the data based on the date range
      goo <- ivermectindata %>%
        filter(
          Premises_No == premise_sub[[i]],
          between(
            Date_Visited, ## visitited between inspections
            foo$Inspection_Date[[j - 1]],
            foo$Inspection_Date[[j]]
          ),
          Date_in_Use < foo$Inspection_Date[[j]]
        ) ## and in use BEFORE inspection

      ## if there are additions of corn between inspections
      if (nrow(goo) > 0) {
        Corn_Count <- goo %>%
          group_by(Station_No) %>%
          ## difference between Date_in_Use & Inspection Date # Amended if retired before inspection
          mutate(
            Feeder_Corn_Days = difftime(
              data.table::fifelse(
                foo$Inspection_Date[[j]] > Date_Retired,
                Date_Retired,
                foo$Inspection_Date[[j]]
              ),
              Date_in_Use,
              units = "days"
            )
          ) %>%
          ungroup() %>%
          summarise(
            Qty_Corn = sum(Corn_Added),
            Qty_TotalCornDays = sum(as.numeric(Feeder_Corn_Days)),
            Qty_CornFeeders = n(),
            Feeder_Density = mean(Qty_CornFeeders / Pasture_Size)
          )

        ## bind to corn data to the original foo
        foo[j, "Qty_Corn"] <- Corn_Count$Qty_Corn
        foo[j, "Qty_TotalCornDays"] <- sum(Corn_Count$Qty_TotalCornDays)
        foo[j, "Qty_CornFeeders"] <- Corn_Count$Qty_CornFeeders
        foo[j, "Feeder_Density"] <- Corn_Count$Feeder_Density
      }
    }
  }
  ##
  corn_feed[[i]] <- foo
}

## generate dataframe of quantified corn to view before binding back to original
Corny_Data <- do.call(bind_rows, corn_feed)


newdat1 <- left_join(
  inspectiondata,
  Corny_Data,
  by = c("Premises_No", "Inspection_Date")
)
write.csv(newdat1, "InspectionCornDataMerge.csv")
### NOTES:
# I would go ahead and clean out all rows where Qty_Corn == before binding back to the
# inspection data (simplifies all amount bound and clears instances where feeders were
# "operational" but no corn was applied between inspections).

#QUESTION 2: Total FEEDER DAYS of ivermectin corn feeder use at premises
#for each "Station_No" compare earliest date of ivermectindata$Corn_Added for
#ivermectindata$Premises_No and count days to inspectiondata$Inspection_Date
#other option: ivermectindata$Date_In_Use to inspectiondata$Inspection_Date for that unique
#"Station_No", then add all together?
#Variable: Qty_TotalCornDays

#QUESTION 3: Number of Active Corn Feeders since last survey
#Use unique? to count the number of unique "Station_No" for that premises ("Premises_No") in that date range (since last inspection)
#Variable: Qty_CornFeeders
