# Data intake and cleaning

library(foundry)

# CFTEP prem data from DIS ----
# my_dataset is set in ~/Documents/.foundry/aliases.yml
# token and hostnames are set in ~/.Renviron
df <- datasets.read_table("my_dataset")
filename <- paste0("Data/CFTEP_7-22_FDE_Survey_Ontology.csv")
readr::write_csv(df, filename)
