library(islandR)

all = read.csv("~/data/sa_report/all_2015.csv")

# Fit the sequence type distribution using the island model
st = st_fit(formula = Label ~ ST,
            non_primary = "Human",
            data = all,
            method="island",
            sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC)


# Fit the attribution model for human cases, estimating separately by location
mod = attribution_ar1(formula = ST ~ 1,
                      time="YearMonth",
                      sampling_dist = st,
                      data=subset(all, Label == "Human"))

humans = manawatu %>% filter(Source == "Human")

mod = attribution(ST ~ Location, st, data=)
