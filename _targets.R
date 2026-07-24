# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tibble",
               "plyr",
               "dplyr",
               "tidyr",
               "readr",
               "stringr",
               "sp",
               "raster",
               "anytime",
               "lubridate",
               "geosphere",
               "operators") # Packages that your targets need for their tasks.
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  #
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
# tar_source()
tar_source("R/1_Reformate.Ontology.R")
tar_source("R/2_Occ.Data.Prep.R")
tar_source("R/Supporting.Functions.R")

# Replace the target list below with your own:
list(

  # TODO script/target for updating DIS data
  # check if different than old data
  # if not different - do nothing

  # Cattle data from DIS - updated regularly
  tar_target(
    name = cftep_dis_data,
    command = "data/CFTEP_7-22_FDE_Survey_Ontology.csv",
    format = "file"
  ),

  # CFT inspection data
  tar_target(
    name = cft_inspection,
    command = "data/CFT_InspectionDataInitial_8.9.20.csv",
    format = "file"
  ),

  tar_target(
    name = scs_ontology_reshaped,
    command = reformat_ontology(cftep_dis_data, cft_inspection)
  ),

  tar_target(
    name = dat_occ_cov,
    command = occupancy_data_prep(cft_inspection, scs_ontology_reshaped)
  )
)
