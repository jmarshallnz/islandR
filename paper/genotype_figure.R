library(islandR)
library(dplyr)
library(tidyr)
library(ggjoy)

# Fit the sequence type distribution using the island model
st_i = replicate(10, {
  s=st_fit(formula = Source ~ ST,
            non_primary = "Human",
            data = manawatu,
            method="island",
            sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
            samples = 1000)
  lapply(1:1000, function(i) { data.frame(ST=s$types, s$sampling_distribution[,,i], Iteration=i) })
  }
  )

# Reformat to a data.frame
st_il <- apply(st_i, 2, function(x) { do.call(rbind, x) })
st_il2 <- lapply(seq_along(st_il), function(x) { cbind(st_il[[x]], Chain=x)})
st_id <- do.call(rbind, st_il2) %>% mutate(Iteration = Iteration + (Chain-1)*1000) %>%
  select(-Chain)

# and using the Dirichlet
st_d = st_fit(formula = Source ~ ST,
            non_primary = "Human",
            data = manawatu,
            method="dirichlet",
            iters=10000)

st_dl <- lapply(1:10000, function(i) { data.frame(ST=as.numeric(st_d$types), st_d$sampling_distribution[,,i], Iteration=i) })
st_dd <- do.call(rbind, st_dl)

# Join them up
st <- bind_rows(cbind(st_id, Model='Island'), cbind(st_dd, Model='Dirichlet')) %>%
  gather(Source, P, -ST, -Iteration, -Model)

final <- st %>% group_by(Iteration, Model, ST) %>%
  mutate(P = P/sum(P)) %>% spread(Model, P)

save(list='final', file='model_compare.Rdata')
load('model_compare.Rdata')

sts <- c(403, 2343, 2026, 474)
plot_dat <- final %>% filter(ST %in% sts) %>%
  gather(Model, P, Dirichlet, Island) %>% ungroup %>% mutate(ST = factor(ST, levels=sts, labels = paste0("ST-", sts)),
                                                             Scale = ifelse(Model == "Dirichlet" & ST %in% c("ST-403","ST-2343"), 2, 1.2))

cairo_pdf("genotype_figure.pdf", width=7, height=4)
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
#        legend.background = element_rect(fill='transparent'))
#        plot.margin = unit(rep(0.5, 4), units = 'cm'))
dev.off()

#       legend.key.height = unit(30, "points"))

