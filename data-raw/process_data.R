library(devtools)

# read in raw data
base_path <- "data-raw"
manawatu <- read.csv(file.path(base_path, "Manawatu.csv"))

# save data within package
devtools::use_data(manawatu, overwrite=TRUE)
