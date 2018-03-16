library(islandR)
library(dplyr)
library(tidyr)

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

#final <- st %>% group_by(Model, Iteration, Source) %>%
#  mutate(P = P/sum(P)) %>% ungroup %>%
#  spread(Model, P)

final <- st %>% group_by(Iteration, Model, ST) %>%
  mutate(P = P/sum(P)) %>% spread(Model, P)

save(list='final', file='model_compare.Rdata')
load('model_compare.Rdata')

# OK, now the Kolmogorov-Smirnov statistic
kolm_smirn <- function(x, y) {
  x <- sort(x)
  y <- sort(y)
  ex <- seq_along(x)/length(x)
  ey <- seq_along(y)/length(y)
  ecdf <- c(ex, ey)[order(c(x, y))]
  max(abs(diff(ecdf)))
}

# test - seems to work
t(replicate(1000, {
x <- rnorm(1000, 0, 1)
y <- rgamma(1000, shape=10)
c(kolm_smirn(x, y), ks.test(x, y)$statistic)
})) -> foo; foo <- data.frame(foo); names(foo) <- c("x", "y")
plot(foo)
summary(lm(y ~ x, data=foo))

diff <- final %>% group_by(Source, ST) %>% summarize(D = kolm_smirn(Island, Dirichlet)) %>%
  group_by(ST) %>%
  summarize(D = sqrt(sum(D^2))) %>%
  arrange(desc(D)) %>%
  left_join(manawatu) %>%
  group_by(ST, D, ASP, GLN, GLT, GLY, PGM, TKT, UNC) %>%
  summarize(Human=sum(Source == "Human"), Source=sum(Source != "Human")) %>%
  arrange(desc(D)) %>% filter(Human > 0) %>% data.frame

sts <- manawatu %>% select(ST, ASP:UNC) %>% unique %>%
  tibble::remove_rownames() %>%
  mutate_at(.vars=vars(ASP:UNC), .funs=as.factor) %>%
  tibble::column_to_rownames(var="ST")

library(cluster)
d <- data.frame(as.matrix(daisy(sts))) %>% tibble::rownames_to_column(var="ST")
dist <- d %>% gather(ST2, D, -ST) %>%
  extract(ST2, into='ST2', regex="X([0-9]+)") %>% filter(D > 0, D < 0.3)


diff %>% left_join(dist %>%
                     mutate(ST = as.numeric(ST), ST2= as.numeric(ST2)), by="ST") %>%
  left_join(manawatu %>% select(ST2=ST, Source) %>% group_by(ST2) %>% summarize(Human=sum(Source=="Human"), Source=sum(Source != "Human")), by="ST2") %>%
  data.frame %>% select(ST, D.x, ST2, D.y, Human.x, Source.x, Human.y, Source.y)

diff %>% left_join(dist %>%
                     mutate(ST = as.numeric(ST), ST2= as.numeric(ST2)), by="ST") %>%
  left_join(manawatu %>% select(ST2=ST, Source) %>% group_by(ST2) %>% summarize(Human=sum(Source=="Human"), Source=sum(Source != "Human")), by="ST2") %>%
  data.frame %>% select(ST, D.x, ST2, D.y, Human.x, Source.x, Human.y, Source.y) %>% pull(ST) %>% unique

# try plotting 403...
plot_dat <- final %>% filter(ST %in% c(704, 10027, 403, 2370, 459, 10001, 3799, 3792, 3798,
                                        7323, 10002, 81, 2343, 2380, 2341, 618, 45, 577, 3728, 393, 3717, 227, 854, 50, 53, 474)) %>% gather(Model, P, Dirichlet, Island)

library(ggplot2)
ggplot(plot_dat) + geom_boxplot(aes(x=Model, y=P, col=Source)) + facet_wrap(~ST)


diff %>% left_join(dist %>%
                     mutate(ST = as.numeric(ST), ST2= as.numeric(ST2)), by="ST") %>%
  left_join(manawatu %>% select(ST2=ST, Source) %>% group_by(ST2) %>% summarize(Human=sum(Source=="Human"), Source=sum(Source != "Human")), by="ST2") %>%
  filter(ST %in% c(403, 2343, 2026, 474, 45))

sts <- c(403, 2343, 2026, 474)
plot_dat <- final %>% filter(ST %in% sts) %>%
  gather(Model, P, Dirichlet, Island) %>% ungroup %>% mutate(ST = factor(ST, levels=sts, labels = paste0("ST-", sts)),
                                                             Scale = ifelse(Model == "Dirichlet" & ST %in% c("ST-403","ST-2343"), 2, 1.2))

library(ggjoy)
#library(tikzDevice)
#tikz("genotype_figure.tex", standAlone = TRUE, width=5, height=6)
cairo_pdf("genotype_figure.pdf", width=7, height=4)
ggplot(plot_dat) + geom_joy(aes(x=P, y=Source, fill=Model, scale=Scale), alpha=0.7, size=0.1, bandwidth=0.01) +
  facet_wrap(~ST, ncol=2) +
  ylab("") +
  scale_x_continuous(name = "P(Source | ST)", limits=c(0,1), expand = c(0,0)) +
  scale_fill_manual(values = c("steelblue", "brown")) +
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

