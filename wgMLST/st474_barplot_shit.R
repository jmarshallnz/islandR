library(readxl)
library(islandR)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(purrr)

load("wgMLST/474_attribution_update.RData")

# now do prediction on each isolate. The way to do this is to combine the predictions from
# above with the predictions from st_fit. Idea is that st_fit gives P(ST | source)
# and the attribution gives P(source), so P(source | ST) = P(ST | source)*P(source) / sum_j P(ST | j) P(j)
# using Bayes Thm.
posterior <- predict(mod, FUN=identity)
post <- posterior %>% mutate(iteration = (identity-1)%/%100+1)

# Human posterior attribution per type: Combine human attribution with type
# distribution for each source via Bayes.
# i.e. P(Source | ST) = P(ST | Source) * P(Source) / sum(P(ST|Source) * P(Source))
#
# Note that Animal attribution will use the same, but with P(Source) being constant
# for each source. In this case, can just use P(Source) = 1
#
# run through for each TYPE.
human_post_types <- list()
for (i in seq_along(st$types)) {
  cat("Processing type", i, "\n")
  p_type <- st$sampling_distribution[i,,]
  colnames(p_type) <- 1:ncol(p_type)
  # for each iteration from st_fit, we have an MCMC run for the attribution, so
  # pull these out and compute via Bayes Thm
  p_type <- p_type %>% as.data.frame %>% tibble::rownames_to_column("Source") %>%
    gather(iteration, p_given_source, -Source, convert=TRUE)

  human_post_types[[i]] <- p_type %>% left_join(post, by=c("iteration", "Source")) %>%
    mutate(bayes_num = p_given_source * p) %>%
    group_by(iteration, identity) %>%
    mutate(bayes_den = sum(bayes_num), p_source_given_st = bayes_num / bayes_den) %>%
    group_by(Source) %>% summarize(v_source_given_st = var(p_source_given_st),
                                   p_source_given_st = mean(p_source_given_st),
                                   Type = st$types[i])
}
human_type_attr <- do.call(rbind, human_post_types)
human_type_attr <- human_type_attr %>% left_join(final %>% select(FILE, Type=ST, source))

# right, now generate our barplots
human_types <- human_type_attr %>% filter(source == "Human") %>%
  group_by(Type) %>% mutate(v_source_given_st = sum(v_source_given_st)) %>% ungroup %>%
  spread(Source, p_source_given_st) %>% select(FILE, Uncertainty=v_source_given_st, EnvWater:`Supplier A`)

write.csv(human_types, "wgMLST/human_barplots_with_uncertainty.csv", row.names=FALSE)
library(RColorBrewer)
col <- brewer.pal(6, "Dark2")
pdf("barplot.pdf", width=4, height=40)
par(mar=c(0.1,5,2,0.1))
b <- barplot(t(as.matrix(human_types %>% select(-Type, -genome, -source, -LabId))), beside = FALSE, border=NA,
             col=col, axes=FALSE, space=1, horiz=TRUE, offset=-0.1, las=1, cex.names = 0.5, xaxs="i", yaxs="i")
mtext(human_types$genome, side=2, line=0.5, las=1, at=b, cex=0.5)
mtext(human_types$LabId, side=2, line=2.5, las=1, at=b, cex=0.5)
par(lend=1)
legend('top', legend=names(human_types %>% select(-Type, -genome, -source, -LabId)), inset=c(-3,-0.01), col=col, xpd=TRUE,
       ncol=3, bty = 'n', cex=0.5, border=NA, lty = 1, lwd = 8)
dev.off()

# NORMALISE the sampling distribution first. i.e. divide by
st_dist <- list()
for (i in seq_along(st$types)) {
  p_type <- st$sampling_distribution[i,,]
  colnames(p_type) <- 1:ncol(p_type)
  p_type <- p_type %>% as.data.frame %>% tibble::rownames_to_column("Source") %>%
    gather(iteration, p_given_source, -Source, convert=TRUE) %>%
    mutate(Type = st$types[i])
  st_dist[[i]] <- p_type
}
st_dist <- do.call(rbind, st_dist) %>%
  group_by(Source, iteration) %>% mutate(p_given_source_norm = p_given_source / sum(p_given_source)) %>%
  group_by(Type, iteration) %>% mutate(p_source_given_st = p_given_source / sum(p_given_source),
                                       p_source_given_st_norm = p_given_source_norm / sum(p_given_source_norm))

source_type_attr <- st_dist %>% group_by(Source, Type) %>% summarize(p_source_given_st = mean(p_source_given_st),
                                                                     p_source_given_st_norm = mean(p_source_given_st_norm))

source_type_attr <- source_type_attr %>% left_join(final %>% select(genome, Type=ST, source))

source_type_attr <- source_type_attr %>% left_join(iso_info) %>% arrange(genome)

animals <- source_type_attr %>% filter(source != "Human") %>% select(-p_source_given_st) %>%
  spread(Source, p_source_given_st_norm) %>% arrange(source, genome)
pdf("animals.pdf", width=4, height=20)
par(mar=c(0.1,5,2,0.1))
b <- barplot(t(as.matrix(animals %>% select(-Type, -genome, -source, -LabId))), beside = FALSE, border=NA,
             col=col, axes=FALSE, space=1, horiz=TRUE, offset=-0.1, las=1, cex.names = 0.5, xaxs="i", yaxs="i")
mtext(animals$genome, side=2, line=0.5, las=1, at=b, cex=0.5)
mtext(animals$source, side=2, line=2.5, las=1, at=b, cex=0.5)
par(lend=1)
legend('top', legend=names(animals %>% select(-Type, -genome, -source, -LabId)), inset=c(-3,-0.01), col=col, xpd=TRUE,
       ncol=3, bty = 'n', cex=0.5, border=NA, lty = 1, lwd = 8)
dev.off()

