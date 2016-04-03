# islandR

The islandR package allows source attribution using the island genomic model (and in the future, other
genomic models).

This allows the attribution of cases of disease to their likely sources.

The assymmetric island model estimates the sampling distribution of genotypes on sources by using
the genetic distance between isolates to infer mutation and recombination rates in addition to
migration rates between sources. This allows improved estimation of the sampling distribution over
and above that achieved using just the prevalence of each type. Further, it allows estimating the
likely prevalence of unobserved genotypes.

## Installation

islandR is not currently available from CRAN, but you can install it from github with:

```R
# install.packages("devtools")
devtools::install_github("jmarshallnz/islandR")
```

## Usage

```R
library(islandR)

# Take a look at the expected data format using an example dataset
head(manawatu)

# Fit the sequence type distribution using the island model
st = st_fit_island(formula = Source ~ ST,
                   sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
                   non_primary = "Human",
                   data = manawatu)

# see some summaries of these
summary(st)
plot(st)

# Fit the attribution model for human cases, estimating separately by location
mod = attribution(ST ~ Location, st, data=subset(manawatu, Source == "Human"))

# Various model summaries
summary(mod)
predict(mod, FUN=mean)

# Posterior predictions including uncertainty
posterior = predict(mod, FUN=identity)
boxplot(p ~ interaction(X, Source), data=posterior, "Posterior attribution")
```

