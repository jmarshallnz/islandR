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

get_pred = function(post, x, formula) {

  # design matrix
  X = model.matrix(formula, data=x)

  # get the prediction
  post_theta = get_var(post, "theta")

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

  # assemble the data in the middle (i.e. filter out the extremes)
  for (j in seq_along(p_pred)) {
    p_max = apply(p_pred[[j]], 2, quantile, 0.975)
    p_min = apply(p_pred[[j]], 2, quantile, 0.025)
    d = cbind(x, p_max, p_min, t(p_pred[[j]]))
    d2 = gather(d, Iteration, Proportion, matches("[0-9]+"))
    p_pred[[j]] = d2 %>% filter(Proportion > p_min, Proportion < p_max)
  }
  return(p_pred)
}
