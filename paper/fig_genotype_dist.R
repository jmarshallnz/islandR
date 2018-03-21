library(islandR)
library(dplyr)
library(tidyr)
library(ggjoy)

# Fit the sequence type distribution using the island model
num_samples <- 10000

st_i = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = manawatu,
              method = "island",
              sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
              samples=num_samples)

st_il <- lapply(1:num_samples, function(i) { data.frame(ST=as.numeric(st_i$types), st_i$sampling_distribution[,,i], Iteration=i) })
st_id <- do.call(rbind, st_il)

# and using the Dirichlet
st_d = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = manawatu,
              method = "dirichlet",
              samples = num_samples)

st_dl <- lapply(1:num_samples, function(i) { data.frame(ST=as.numeric(st_d$types), st_d$sampling_distribution[,,i], Iteration=i) })
st_dd <- do.call(rbind, st_dl)

# Join them up
st <- bind_rows(cbind(st_id, Model='Island'), cbind(st_dd, Model='Dirichlet')) %>%
  gather(Source, P, -ST, -Iteration, -Model)

final <- st %>% group_by(Iteration, Model, ST) %>%
  mutate(P = P/sum(P)) %>% spread(Model, P)

#save(list='final', file='model_compare.Rdata')
#load('model_compare.Rdata')

# Filter out the STs to plot
sts <- c(403, 2343, 2026, 474)
plot_dat <- final %>% filter(ST %in% sts) %>%
  gather(Model, P, Dirichlet, Island) %>% ungroup %>% mutate(ST = factor(ST, levels=sts, labels = paste0("ST-", sts)),
                                                             Scale = ifelse(Model == "Dirichlet" & ST %in% c("ST-403","ST-2343"), 2, 1.2))

# and do the plot
pdf("fig_genotype_dist.pdf", width=7, height=4)
ggplot(plot_dat) + geom_joy(aes(x=P, y=Source, fill=Model, scale=Scale), alpha=0.7, size=0.1, bandwidth=0.01) +
  facet_wrap(~ST, ncol=2) +
  ylab("") +
  scale_x_continuous(name = "P(Source | ST)", limits=c(0,1), expand = c(0,0)) +
  scale_fill_manual(values = c("steelblue2", "brown")) +
  theme_bw() +
  theme(text = element_text(family="Times"),
        legend.position = c(0.92,0.895),
        legend.box.background = element_rect(),
        legend.margin = margin(1.5,3,3,3),
        legend.title = element_blank(),
        axis.text.x = element_text(hjust=c(-0.1,rep(0.5, 3), 1.1)))
dev.off()
