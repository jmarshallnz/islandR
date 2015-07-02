# islandR

The islandR package allows source attribution using the island genomic model.

At the moment, the package encapsulates a genomic data fit to current NZ *Campylobacter jejuni* data and provides posterior
estimates so that further models can be built without the overhead of the genomic model.

The genomic model was fit using the genomic_loop branch of https://github.com/jmarshallnz/island

At a later stage, the goal is to perform all source attribution fitting directly (i.e. merge in the island model code).

## Installation

islandR is not currently available from CRAN, but you can install it from github with:

```R
# install.packages("devtools")
devtools::install_github("jmarshallnz/islandR")
```

## Usage

```R
library(islandR)

# find the number of posterior samples for each type
samples <- get_num_samples()

# retrieve the genotypes for which we have attribution results
isolates <- get_genotypes()

# get the attribution for a particular isolate
samples <- get_source_probability_sample(genotype=474)
```

