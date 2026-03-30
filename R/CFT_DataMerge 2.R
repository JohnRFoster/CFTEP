setwd("path")
inspectiondata <- read.csv(
      "CFT_InspectionData_8.9.20.csv",
      stringsAsFactors = TRUE
)
premisesdata <- read.csv(
      "C:/Users/Andrew.Browne/OneDrive - USDA/CFT/TickTrackerData/TAHC+LairdTickTrackerExport-PremInspTreat-2019-09-19_Premises.csv",
      stringsAsFactors = TRUE
)
treatmentdata <- read.csv(
      "CFT_TreatmentData_8.9.20.csv",
      stringsAsFactors = TRUE
)
newdata <- read.csv(
      "C:/Users/Andrew.Browne/OneDrive - USDA/CFT/TickTrackerData/CFTEP_Premises_History.csv",
      stringsAsFactors = TRUE
)
newpremiseslceaned <- read.csv(
      "G:/CFTEP/Data/DatasetCleaning/Tick+Tracker+Prems+Cleaned.csv",
      stringsAsFactors = TRUE
)

#merge Pasture data to inspection dataset#
inspectionmerge <- merge(
      inspectiondata,
      premisesdata[, c(
            "Pasture_Name",
            "Pasture_Latitude",
            "Pasture_Longitude",
            "Premises_Latitude",
            "Premises_Longitude",
            "County_Name",
            "Area_Type",
            "Pasture_Qty_Acres",
            "Premises_Qty_Acres",
            "CFTPID"
      )],
      by = "Pasture_Name",
      all.x = TRUE
)
write.csv(inspectionmerge, "CFT_InspectionDataInitial_8.9.20.csv")
###Data needs to have ".x" columns deleted"
#Way to have it write over the columns?  Not vital#

#merge Pasture data to treatment dataset#
treatmentmerge <- merge(
      treatmentdata,
      premisesdata[, c(
            "Pasture_Name",
            "Pasture_Latitude",
            "Pasture_Longitude",
            "Premises_Latitude",
            "Premises_Longitude",
            "County_Name",
            "Area_Type",
            "Pasture_Qty_Acres",
            "Premises_Qty_Acres",
            "CFTPID"
      )],
      by = "Pasture_Name",
      all.x = TRUE
)
write.csv(treatmentmerge, "CFT_TreatmentDataInitial_8.9.20.csv")
###Data needs to have ".x" columns deleted"
#Way to have it write over the columns?  Not vital#

library(ggmap)
register_google(key = "GOOGLE_API_KEY")
texas <- c(long = -99.735462, lat = 30.638248)
texasmap <- get_map(location = texas, zoom = 5, scale = 1, maptype = "terrain")
ggmap(texasmap) + geom_point(aes(longitude, latitude), data = newdata)
ggmap(texasmap) +
      geom_point(aes(longitude, latitude), data = newpremiseslceaned)


newdata$completion_date < as.Date(newdata$completion_date)
dates <- summary(newdata$completion_date)
range(newdata$completion_date, na.rm = TRUE)
newdata <- as.date.frame(newdata$completion_date, na.rm = TRUE)
min(newdata$completion_date, na.rm = TRUE)
dates <- as.integer(diff(range(as.Date(newdata$completion_date))))
range(as.Date(newdata$completion_date))
# "2001-06-13" "2021-01-05"
summary(newdata$county)
