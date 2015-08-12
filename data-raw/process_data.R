library(dplyr)

base_path <- "data-raw"

# read in the STs
humans <- read.csv(file.path(base_path, "human_types.csv"))
# filter out humans and years
humans <- humans %>% filter(Year >= 2005 & Year <= 2014 & Source=="Human")

# loop through and load all the data files in
# load any file named phi_[0-9]_output_[0-9].txt;
# [0-9] means any digit ranged 0-9; first [0-9]= number of iteration
# second [0-9]= number of chain
# make sure there is no hidden files

files <- list.files(base_path, pattern = "phi_[0-9]+_output_[0-9]+.txt")

# now split up based on chain and iteration

data <- data.frame(file=files, iteration, chain)

# grab out the data for these chains
# refer paragraph in Tracing the Source of Campylobacteriosis p.6:"
# taking each non-human group in turn, half the isolates were removed
# leaving the other half in the pool and the number of genotypes 
# unique to the removed isolates was calculated. A set of human
# isolates was drawn of equal number, and the number of unique geno-
# types calculated relative to the same pool"

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
# transform the values in sources column
# paste0(x) is short for paste(x, sep="")
# specify the column sources

attr_cols <- names(all) %in% paste0("Source",0:3)

# select values from column sources only and transform exponentially
# then store these columns in the original data frame

all[,attr_cols] <- exp(all[,attr_cols])

# modified from rowSums to colSums in order to fix p(j|i) to p(i|j)
# where j represents source and i denotes ST

all[,attr_cols] <- all[,attr_cols] / colSums(all[,attr_cols])

# reorder + rename
all <- all %>% select(ST, everything()) %>% rename(ASP=Loci0, GLN=Loci1, GLT=Loci2,
                                                   GLY=Loci3, PGM=Loci4, TKT=Loci5,
                                                   UNC=Loci6, Poultry=Source0, Ruminants=Source1,
                                                   Water=Source2, Other=Source3)

# Iteration has to be calculated as the files are combined into a
# single chain! (by interleaving) 
all$Iteration <- (as.numeric(as.factor(all$Iteration))-1) * length(unique(all$Chain)) + as.numeric(as.factor(all$Chain))

genotype_attribution <- all %>% select(-Chain)


# save data
devtools::use_data(genotype_attribution, overwrite=TRUE)
