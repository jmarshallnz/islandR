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


# OLD STUFF BELOW HERE...

library(MASS)
library(stringr)
library(dplyr)
library(tidyr)
library(lubridate)
library(islandR)



# do some plots
par(mfrow=c(3,1))
plot_traces(post, "p")
par(mfrow=c(3,1))
plot_traces(post, "theta")

par(mfrow=c(1,1))
pairs(get_var(post, "theta"), labels=paste("theta", 1:20))

post_theta = get_var(post, "theta")

#' design matrix
X = reduced.matrix

p_pred = list()
for (j in seq_along(post_theta)) {
  p_pred[[j]] = t(exp(X %*% t(post_theta[[j]])))
  colnames(p_pred[[j]]) = apply(X, 1, function(x, y) { paste(y[x == 1], collapse=":")}, colnames(X))
}
p_pred[[length(p_pred)+1]] = matrix(1, nrow(p_pred[[1]]), ncol(p_pred[[1]]))

p_sum = matrix(0, nrow(p_pred[[1]]), ncol(p_pred[[1]]))
for (j in seq_along(p_pred)) {
  p_sum = p_sum + p_pred[[j]]
}
for (j in seq_along(p_pred)) {
  p_pred[[j]] = p_pred[[j]] / p_sum
}
names(p_pred) = colnames(x$sampling_distribution[,,1])

par(mfrow=c(2,2))
for (i in 1:ncol(p_pred[[1]])) {
  for (j in 1:length(p_pred))
    plot(p_pred[[j]][,i], type="l")
}

par(mfrow=c(1,1))
for (i in 1:ncol(p_pred[[1]])) {
  plot(density(p_pred[[1]][,i]), xlim=c(0,1), type="l", main=colnames(p_pred[[1]])[i], ylim=c(0,30))
  for (j in 2:length(p_pred)) lines(density(p_pred[[j]][,i]), col=j)
  legend("topright", legend=names(p_pred), col=1:length(p_pred), lty=1)
}

# Now incorporating the random effects...
post_p = get_var(post, "p")
post_p[[length(post_p)+1]] = matrix(0, nrow(post_p[[1]]), ncol(post_p[[1]]))

for (j in seq_along(post_p)) {
  post_p[[j]] = exp(post_p[[j]])
}

p_sum = matrix(0, nrow(post_p[[1]]), ncol(post_p[[1]]))
for (j in seq_along(post_p)) {
  p_sum = p_sum + post_p[[j]]
}
for (j in seq_along(post_p)) {
  post_p[[j]] = post_p[[j]] / p_sum
}
par(mfrow=c(1,1))
plot(density(post_p[[1]]), xlim=c(0,1), type="l")
for (j in 2:length(post_p)) lines(density(post_p[[j]]), col=j)


















# TESTING!
x0      = NULL
formula = ~ 1
hdat = dat %>%
  mutate(Month = month(as.Date(Sampled.Date)), YearMonth = Month + (Year-2005)*12) %>%
  filter(Source == "Human", Year >= 2005, Year <= 2014) %>% dplyr::select(ST, UR_bool, YearMonth) %>% group_by(ST) %>% summarise(Number=n())

sts_available = get_genotypes()$ST
hdat$ST = match(hdat$ST, sts_available)

hum = list()
hum[[1]] = as.matrix(hdat)
# END TESTING!


# for each ST, create phi from posterior of previous runs
post = NULL;

nj = sj:ej #get_num_samples()
#set.seed(1)
ar = c(0,0)
for (j in nj) {
  cat("Processing prior", j, "of", sj, " to ", ej, "\n")

  phi = matrix(NA, length(sts_available), length(source_labels))
  for (st in seq_along(sts_available)) {
#    phi[st,] =  as.numeric(get_source_probability_sample(sts_available[st], j))
    phi[st,] =  as.numeric(4*colMeans(get_source_probability_sample(sts_available[st])))
  }

  # TODO: We'd ideally pass in the model matrix
  # now run our island model MCMC loop
  system.time({iter = mcmc(hum, x0, formula, phi, 20000)})
  post = c(post, iter$post)
  ar   = ar + iter$ar
}

pdf("output.pdf", width=11, height=8)
# do some plots
par(mfrow=c(3,1))
plot_traces(post, "p")

par(mfrow=c(1,1))
pairs(get_var(post, "p"), labels=paste("p", 1:3))

post_theta = get_var(post, "theta")

#' design matrix
if (is.null(x0))
  x0 = data.frame(dummy=1)
X = model.matrix(formula, data=x0)

p_pred = list()
for (j in seq_along(post_theta)) {
  p_pred[[j]] = t(exp(X %*% t(post_theta[[j]])))
  colnames(p_pred[[j]])
}
p_pred[[length(p_pred)+1]] = matrix(1, nrow(p_pred[[1]]), ncol(p_pred[[1]]))

p_sum = matrix(0, nrow(p_pred[[1]]), ncol(p_pred[[1]]))
for (j in seq_along(p_pred)) {
  p_sum = p_sum + p_pred[[j]]
}
for (j in seq_along(p_pred)) {
  p_pred[[j]] = p_pred[[j]] / p_sum
}

par(mfrow=c(2,2))
for (j in 1:4)
  plot(p_pred[[j]], type="l")

par(mfrow=c(1,1))
plot(density(p_pred[[1]]), xlim=c(0,1), type="l")
for (j in 2:4) lines(density(p_pred[[j]]), col=j)

dev.off()

# Now incorporating the random effects...
post_p = get_var(post, "p")
post_p[[length(post_p)+1]] = matrix(0, nrow(post_p[[1]]), ncol(post_p[[1]]))

for (j in seq_along(post_p)) {
  post_p[[j]] = exp(post_p[[j]])
}

p_sum = matrix(0, nrow(post_p[[1]]), ncol(post_p[[1]]))
for (j in seq_along(post_p)) {
  p_sum = p_sum + post_p[[j]]
}
for (j in seq_along(post_p)) {
  post_p[[j]] = post_p[[j]] / p_sum
}
par(mfrow=c(1,1))
plot(density(post_p[[1]]), xlim=c(0,1), type="l")
for (j in 2:4) lines(density(post_p[[j]]), col=j)

append = paste0("_",sj,"_",ej)
saveRDS(post, paste0("post",append,".rds"))

pdf(paste0("diag_traces", append, ".pdf"), width=11, height=8)

plot_traces(post, "tau")
plot_traces(post, "p")
dev.off()

pdf(paste0("diag_density", append, ".pdf"), width=11, height=8)
plot_density(post, "theta", prior_fun = function(x) { dnorm(x, 0, 1/sqrt(0.1)) })
plot_density(post, "tau", prior_fun = function(x) { dgamma(x, 0.1, 0.1) })
plot_density(post, "rho", prior_fun = function(x) { dnorm(x, 0, 1/sqrt(16)) })
dev.off()

post_theta = get_var(post, "theta")


# TESTING STUFF

post_theta = get_var(post, "p")
# END TESTING STUFF

# figure out which terms are signicant by generating P-values

pmean = NULL
for (var in post_theta) {
  pmean = rbind(pmean, colMeans(var))
}

pse = NULL
for (var in post_theta) {
  pse = rbind(pse, apply(var, 2, sd))
}

pp = NULL
for (var in post_theta) {
  pp = rbind(pp, pmin(colSums(var < 0),colSums(var > 0)) / nrow(var))
}

# compute predictions for each season/location/source

#' design matrix
if (is.null(x0))
  x0 = data.frame(dummy=1)
X = model.matrix(formula, data=x0)

p_pred = list()
for (j in seq_along(post_theta)) {
  p_pred[[j]] = t(exp(X %*% t(post_theta[[j]])))
  colnames(p_pred[[j]])
}
p_pred[[length(p_pred)+1]] = matrix(1, nrow(p_pred[[1]]), ncol(p_pred[[1]]))

p_sum = matrix(0, nrow(p_pred[[1]]), ncol(p_pred[[1]]))
for (j in seq_along(p_pred)) {
  p_sum = p_sum + p_pred[[j]]
}
for (j in seq_along(p_pred)) {
  p_pred[[j]] = p_pred[[j]] / p_sum
}

plot(density(p_pred[[1]]), xlim=c(0,1), type="l")
for (j in 2:4) lines(density(p_pred[[j]]), col=j)


# right, now plot these distributions
pdf(paste0("attribution", append, ".pdf"), width=8, height=11)
for (j in seq_along(p_pred)) {
  p_max = apply(p_pred[[j]], 2, quantile, 0.975)
  p_min = apply(p_pred[[j]], 2, quantile, 0.025)
  d = cbind(x0, p_max, p_min, t(p_pred[[j]]))
  par(mfrow=c(2,1), mar=c(4,4,2,2))
  d2 = gather(d, Iteration, Proportion, -(Time:p_min)) # TODO: ideally this wouldn't depend on other stuff

  # TODO: This needs to be updated from the laptop so we actually produce the plots we want (or maybe
  #       from the talk repo???)
  boxplot(Proportion ~ interaction(Season, Intervention), data=d2 %>% filter(Loc == 0, Proportion < p_max, Proportion > p_min), main=paste(source_labels[as.character(j)], "Urban", sep=" - "), names = rep(c("Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec"), 2), ylim=c(0,1), cex=0.5, cex.axis=0.9)
  mtext(c("2005-2007", "2008-2014"), side = 1, line = 2.5, at = c(2.5, 6.5))
  boxplot(Proportion ~ interaction(Season, Intervention), data=d2 %>% filter(Loc == 1, Proportion < p_max, Proportion > p_min), main=paste(source_labels[as.character(j)], "Rural", sep=" - "), names = rep(c("Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec"), 2), ylim=c(0,1), cex=0.5, cex.axis=0.9)
  mtext(c("2005-2007", "2008-2014"), side = 1, line = 2.5, at = c(2.5, 6.5))
}
dev.off()

# OK, now let's look at the predicted densities...
post_p = get_var(post, "p")
post_p[[length(post_p)+1]] = matrix(0, nrow(post_p[[1]]), ncol(post_p[[1]]))

for (j in seq_along(post_p)) {
  post_p[[j]] = exp(post_p[[j]])
}

p_sum = matrix(0, nrow(post_p[[1]]), ncol(post_p[[1]]))
for (j in seq_along(post_p)) {
  p_sum = p_sum + post_p[[j]]
}
for (j in seq_along(post_p)) {
  post_p[[j]] = post_p[[j]] / p_sum
}

pdf(paste0("attribution_time", append, ".pdf"), width=8, height=11)
for (j in seq_along(post_p)) {
  mu = apply(post_p[[j]], 2, mean)
  li = apply(post_p[[j]], 2, quantile, 0.025)
  ui = apply(post_p[[j]], 2, quantile, 0.975)
  par(mfrow=c(2,1), mar=c(4,4,2,2))
  times = 1:n_times
  plot(NULL, xlim=range(times), ylim=c(0,1), type="n", main=paste(source_labels[as.character(j)], "Urban", sep=" - "))
  polygon(c(times, rev(times)), c(ui[times], rev(li[times])), col="grey80", border=NA)
  lines(times, mu[times], lwd=2)

  plot(NULL, xlim=range(times), ylim=c(0,1), type="n", main=paste(source_labels[as.character(j)], "Rural", sep=" - "))
  polygon(c(times, rev(times)), c(ui[times+n_times], rev(li[times+n_times])), col="grey80", border=NA)
  lines(times, mu[times+n_times], lwd=2)
}
dev.off()

# the following is actual analysis
# pdf("urban_rural.pdf", width=11, height=8)
# labs = c("Urban", "Rural")
# for (k in 1:2) {
#
# # and on the natural scale
# e_theta = exp(cbind(post_theta[,c(1,3,5)+k-1], 0))
# e_theta = e_theta / rowSums(e_theta)
# plot(NULL, xlim=c(1,nrow(e_theta)), ylim=range(e_theta), type="n", main=labs[k])
# for (i in 1:ncol(e_theta))
#   lines(1:nrow(e_theta), e_theta[,i], col=i)
#
# dens = list()
# ylim = xlim = NA
# for (i in 1:ncol(e_theta)) {
#   dens[[i]] = density(e_theta[,i])
#   ylim = range(ylim, dens[[i]]$y, na.rm=TRUE)
#   xlim = range(xlim, dens[[i]]$x, na.rm=TRUE)
# }
# plot(NULL, xlim=xlim, ylim=ylim, type="n", main=labs[k])
# for (i in seq_along(dens)) {
#   lines(dens[[i]]$x, dens[[i]]$y, col=i)
#   # prior
#   x = seq(xlim[1], xlim[2], length.out = 100)
#   lines(x, dnorm(x, 0, 1/sqrt(0.1)), col=i, lty="dotted")
# }
# }
# dev.off()
