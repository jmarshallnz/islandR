library(islandR)

# reproducible
set.seed(123)

# some helpful fits with which to test stuff
st = st_fit(Source ~ ST, sequences=~ASP+GLN+GLT+GLY+PGM+TKT+UNC, non_primary = "Human", method="island", data=manawatu, iters=100)
