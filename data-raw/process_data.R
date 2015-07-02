library(dplyr)

base_path <- "data-raw"

# read in the STs
humans <- read.csv(file.path(base_path, "human_types.csv"))
# filter out humans and years
humans <- humans %>% filter(Year >= 2005 & Year <= 2014 & Source=="Human")

# loop through and load all the data files in
files <- list.files(base_path, pattern = "phi_[0-9]+_output_[0-9]+.txt")

# now split up based on chain and iteration

data <- data.frame(file=files, iteration, chain)

# grab out the data for these chains
all <- list()
for (f in seq_along(files)) {
  d <- read.table(file.path(base_path, files[f]), header=T)
  d$Iteration <- as.numeric(gsub("phi_([0-9]+).*", "\\1", files[f]))
  d$Chain     <- as.numeric(gsub("phi_[0-9]+_output_([0-9]+).txt", "\\1", files[f]))
  d$ST <- humans$ST
  all[[f]] <- unique(d %>% select(-Isolate) %>% arrange(ST))
}
all <- do.call(rbind, all)

# transform
attr_cols <- names(all) %in% paste0("Source",0:3)
all[,attr_cols] <- exp(all[,attr_cols])
all[,attr_cols] <- all[,attr_cols] / rowSums(all[,attr_cols])

# reorder + rename
all <- all %>% select(ST, everything()) %>% rename(ASP=Loci0, GLN=Loci1, GLT=Loci2,
                                                   GLY=Loci3, PGM=Loci4, TKT=Loci5,
                                                   UNC=Loci6, Poultry=Source0, Ruminants=Source1,
                                                   Water=Source2, Other=Source3)

all$Iteration <- (as.numeric(as.factor(all$Iteration))-1) * length(unique(all$Chain)) + as.numeric(as.factor(all$Chain))

genotype_attribution <- all %>% select(-Chain)


# save data
devtools::use_data(genotype_attribution, overwrite=TRUE)
