# script to produce barplot figure...

library(meshblocknz)
library(dplyr)
# 1. Load meshblocknz and compute population by rurality by year
#    Technically want this only for the meshblocks that we're actually sampling from
mb2013 %>% filter(DHB_name == "Mid Central") %>% nrow

# 2. Once done, compute number of cases by rurality by year

# Now produce a ggplot barplot...

