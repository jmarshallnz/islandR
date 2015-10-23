library(stringr)
library(dplyr)

# load in our dataset (relative to current working directory)
db_file <- "running/20150603.csv"
db <- read.csv(db_file)

# run the number of cases we want and so on

# options are:
# 1. Max number of alleles to impute
# 2. Source list
# 3. Output folder (defaults to current directory + "out_[imputed]_[sources]")

source_file <- "4_source_imputed"

sources <- read.csv(file.path("running", paste0(source_file, ".csv")), colClasses="character")
source_map <- as.numeric(sources$Number)
names(source_map) <- sources$DataSource

source_label_map <- unique(sources %>% dplyr::select(Number, Label))
source_labels <- str_replace(source_label_map$Label, "\\\\n", "\n")
names(source_labels) <- source_label_map$Number

# input parameters
alleles_to_impute <- max(c(suppressWarnings(as.numeric(sources$Imputed)), 0), na.rm=T)

human   <- "Human"

# Setup data
db <- db %>% filter(Imputed <= alleles_to_impute)

humans <- db %>% filter(Source == human, Year <= 2014)

library(islandR)
sts_available = get_genotypes()$ST

# now map our humans to rows of phi
humans$ST = match(humans$ST, sts_available)

# table up the humans to speed things up a bit
humans = as.matrix(humans %>% filter(!is.na(UR_bool)) %>% mutate(UR_bool = ifelse(UR_bool == "Rural", 1, 0)) %>% group_by(UR_bool, Quarter, ST) %>% summarise(n = n()))

hum = list()
n_times = max(humans[,2])
n_loc   = length(unique(humans[,1]))
count   = 0;
for (j in unique(humans[,1])) {
  for (i in 1:n_times) {
    count = count + 1;
    hum[[count]] = as.matrix(humans[humans[,1] == j & humans[,2] == i, 3:4])
  }
}

# for each ST, create phi from posterior of previous runs
post = NULL;

nj = sj:ej #get_num_samples()
#set.seed(1)
ar = c(0,0)
for (j in nj) {
  cat("Processing prior", j, "of", sj, " to ", ej, "\n")

  phi = matrix(NA, length(sts_available), length(source_labels))
  for (st in seq_along(sts_available)) {
    phi[st,] =  as.numeric(get_source_probability_sample(sts_available[st], j))
  }

  # TODO: We'd ideally pass in the model matrix
  # now run our island model MCMC loop
  system.time({iter = mcmc(hum, phi, 20000)})
  post = c(post, iter$post)
  ar   = ar + iter$ar
}

saveRDS(post, paste0("post_",sj,"_",ej,".rds"))

# take a look at the posteriors...

get_var = function(post, variable) {
  ncol = 1;
  if (!is.null(dim(post[[1]][[variable]]))) {
    ncol = ncol(post[[1]][[variable]])
    lapply(1:ncol, function(col) { do.call(rbind, lapply(post, function(x) { x[[variable]][,col] })) })
  } else {
    list(do.call(rbind, lapply(post, function(x) { as.numeric(x[[variable]]) })))
  }
}

# traces...
plot_traces = function(post, variable) {
  post_list = get_var(post, variable)
  for (j in seq_along(post_list)) {
    post_var = post_list[[j]]
    name = variable
    if (j > 1)
      name = paste(name, j)
    plot(NULL, xlim=c(1,nrow(post_var)), ylim=range(post_var), type="n", main=name, ylab=variable, xlab="iteration")
    for (i in 1:ncol(post_var))
      lines(1:nrow(post_var), post_var[,i], col=i)
  }
}

plot_density = function(post, variable, prior_fun) {
  post_list = get_var(post, variable)
  for (j in seq_along(post_list)) {
    post_var = post_list[[j]]
    name = variable
    if (j > 1)
      name = paste(name, j)
    dens = list()
    ylim = xlim = NA
    for (i in 1:ncol(post_var)) {
      dens[[i]] = density(post_var[,i])
      ylim = range(ylim, dens[[i]]$y, na.rm=TRUE)
      xlim = range(xlim, dens[[i]]$x, na.rm=TRUE)
    }
    plot(NULL, xlim=xlim, ylim=ylim, type="n", main=name, xlab=variable, ylab="density")
    for (i in 1:ncol(post_var)) {
      lines(dens[[i]]$x, dens[[i]]$y, col=i)
      # prior
      x = seq(xlim[1], xlim[2], length.out = 100)
      lines(x, prior_fun(x), col=i, lty="dotted")
    }
  }
}

#pdf("diag_traces.pdf", width=11, height=8)
#plot_traces(post, "theta")
#plot_traces(post, "tau")
plot_traces(post, "rho")
#dev.off()

#pdf("diag_density.pdf", width=11, height=8)
plot_density(post, "theta", prior_fun = function(x) { dnorm(x, 0, 1/sqrt(0.1)) })
plot_density(post, "tau", prior_fun = function(x) { dgamma(x, 0.1, 0.1) })
plot_density(post, "rho", prior_fun = function(x) { dnorm(x, 0, 1/sqrt(16)) })
#dev.off()

post_theta = get_var(post, "theta")

#pdf("diag_correlation.pdf", width=11, height=8)
#pairs(cbind(post_theta, post_rho, post_tau), labels=c(paste0("theta", 1:ncol(post_theta)), "rho", "tau"))
#dev.off()

#post_p   = do.call(rbind, lapply(post, function(x) { as.numeric(x[["p"]]) }))

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
x0 = expand.grid(Time = 1:n_times, Loc = 0:(n_loc-1))

x0$Season = as.factor((x0$Time - 1) %% 4 + 1)
x0$Intervention = as.factor(ifelse(x0$Time > 12, 1, 0))
x0$Loc = as.factor(x0$Loc)

x0 = expand.grid(Season = as.factor(1:4), Intervention=as.factor(0:1), Loc=as.factor(0:1))
#' design matrix
X = model.matrix(~ Season*Intervention*Loc, data=x0)

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

# right, now plot these distributions
pdf("attribution.pdf", width=8, height=11)
library(tidyr)
for (j in seq_along(p_pred)) {
  p_max = apply(p_pred[[j]], 2, quantile, 0.975)
  p_min = apply(p_pred[[j]], 2, quantile, 0.025)
  d = cbind(x0, p_max, p_min, t(p_pred[[j]]))
  par(mfrow=c(2,1), mar=c(4,4,2,2))
  d2 = gather(d, Iteration, Proportion, -(Season:p_min))
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

pdf("attribution_time.pdf", width=8, height=11)
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

for (j in seq_along(post_p)) {
  #
  #  li = apply(post_p[[j]], 2, quantile, 0.025)
  #  ui = apply(post_p[[j]], 2, quantile, 0.975)
  #  mu = post_[[[j]]][1,]
  #  mu = post_[[[j]]][1,]
  par(mfrow=c(2,1), mar=c(4,4,2,2))
  times = 1:n_times
  plot(NULL, xlim=range(times), ylim=c(0,1), type="n", main=paste(source_labels[as.character(j)], "Urban", sep=" - "))
  #  polygon(c(times, rev(times)), c(ui[times], rev(li[times])), col="grey80", border=NA)
  for (i in 1:nrow(post_p[[j]])) {
    mu = post_p[[j]][i,]
    lines(times, mu[times])
  }
  mu = apply(post_p[[j]], 2, mean)
  lines(times, mu[times], lwd=2)

  plot(NULL, xlim=range(times), ylim=c(0,1), type="n", main=paste(source_labels[as.character(j)], "Rural", sep=" - "))
  #  polygon(c(times, rev(times)), c(ui[times+n_times], rev(li[times+n_times])), col="grey80", border=NA)
  for (i in 1:nrow(post_p[[j]])) {
    mu = post_p[[j]][i,]
    lines(times, mu[times+n_times])
  }
  mu = apply(post_p[[j]], 2, mean)
  lines(times, mu[times+n_times], lwd=2)
}

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
