library(devtools)

# read in raw data
base_path <- "data-raw"
manawatu <- read.csv(file.path(base_path, "Manawatu.csv"))

# reorder the UR2006_name factor
library(dplyr)
order <- manawatu %>% select(UR2006_num, UR2006_name) %>% unique %>% arrange(UR2006_num) %>% pull(UR2006_name)
manawatu$UR2006_name = factor(manawatu$UR2006_name, levels=order)

# save data within package
devtools::use_data(manawatu, overwrite=TRUE)
