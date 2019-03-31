library(readxl)
library(islandR)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(purrr)

# experimental work at fitting island model to whole genome MLST
samples <- 100
thin    <- 100
burnin  <- 10

#wgmlst <- read.table('wgMLST/allele_profiles.txt')
#names(wgmlst)[1] <- "genome"
wgmlst <- read.csv('wgMLST/wgmlst.csv')
wgmlst <- wgmlst %>% rename(genome=ID) %>%
  mutate(genome = sprintf('%03d.fas', genome))

# remove columns that are useless (zero variance)
v_col <- apply(wgmlst[,-1],2,var)
wgmlst <- wgmlst[, c(1,which(v_col > 0)+1)]

# add an ST column, being the pasted copy of everything else
wgmlst$ST <- as.numeric(as.factor(apply(wgmlst[,-1],1,paste,collapse='_')))

# read in isolate information
iso <- read_excel('wgMLST/st474_530isolates.xlsx') %>%
  mutate(Date = as.Date(ifelse(is.na(SampledDate), as.character(TestedDate), as.character(SampledDate)))) %>%
  select(genome, source=SA_model_source, Date) %>%
  mutate(source = ifelse(source == 'Supplier A', ifelse(Date < '2008-01-01', 'A_before', 'A_after'), source))

# join the wgmlst up to the isolate information
wgmlst <- wgmlst %>% left_join(iso) %>%
  mutate(source = fct_collapse(source,
                               Other = c("Cat_dog_pet", "Pig", "Wild_bird_other"),
                               OtherPoultry = c("Supplier_other", "Supplier B", "Spent_hen"),
                               Supplier_A = c("A_before", "A_after"),
                               EnvWater = "Environmental water"))

table(wgmlst$source, useNA='always')
table(iso$source, useNA='always')
# filter out stuff we don't want
final <- wgmlst %>% filter(!is.na(source)) # %>% filter(source %in% c('A_after', 'A_before', 'Environmental water', 'Human', 'Supplier_other', 'Cattle', 'Sheep'))
table(final$source, useNA='always')
nrow(final)

# now try and run islandR...
# we're checking here if the seed makes a difference. It shouldn't ofcourse!
seeds <- c(2,3,5)
sts <- lapply(seeds, function(x) {
             set.seed(x)
             st_fit(formula = source ~ ST,
                    non_primary = "Human",
                    data = final,
                    method="island",
                    sequences = formula(terms(~ . - Date - genome - source - ST, data=final, simplify=TRUE)),
                    samples=samples, burnin=burnin, thin=thin)
              })

bind_sampling_dists <- function(sts) {
  ltypes <- lapply(sts, function(x) { x$types })
  lsources <- lapply(sts, function(x) { x$sources })
  # check equal
  if (length(unique(ltypes)) > 1) {
    stop("Different types for each model fit")
  }
  types <- ltypes[[1]]
  if (length(unique(lsources)) > 1) {
    stop("Different sources for each model fit")
  }
  sources <- lsources[[1]]
  iters <- unlist(lapply(sts, function(x) { dim(x$sampling_distribution)[3] }))

  # construct new one
  st <- sts[[1]]
  # fill them in
  st$evolution_params <- do.call(rbind, lapply(sts, function(x) { x$evolution_params }))

  st$sampling_distribution <- array(dim = c(length(types), length(sources), sum(iters)),
                                    dimnames = list(types, sources, 1:sum(iters)))
  start <- 0
  for (i in seq_along(sts)) {
    st$sampling_distribution[,,start + 1:iters[i]] <- sts[[i]]$sampling_distribution
    start <- start + iters[i]
  }
  st
}

# combine all these together to do some trace plots
evol <- lapply(seq_along(seeds), function(x) { sts[[x]]$evolution_params %>% mutate(Seed=seeds[x])})
all <- bind_rows(evol) %>% gather(Param, Value, -Seed, -Iteration)

ggplot(all) +
  geom_line(aes(Iteration, Value, col=factor(Seed))) +
  facet_wrap(~Param, scales='free_y')

st <- bind_sampling_dists(sts)

# filter out the humans
humans <- final %>%
  filter(source == "Human") %>%
  mutate(Time = ifelse(Date < '2008-01-01', "Before", "After"))

# for a laugh, attempt to attribute (ignoring time for now)
mod = attribution(ST ~ 1, st, data=humans, iterations=10000, burnin=1000, thinning=100)
# check attribution
predict(mod, FUN=mean)

modtime = attribution(ST ~ Time, st, data=humans, iterations=10000, burnin=1000, thinning=100)

save(final, humans, st, mod, modtime, file="wgMLST/474_attribution.RData")

