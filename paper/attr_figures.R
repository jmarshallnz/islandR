library(islandR)
library(dplyr)
library(tidyr)
library(ggplot2)

# filter out those where we don't have a location, and have a factor version
# of UR2006_num
dataset <- manawatu %>% filter(Source != "Human" | !is.na(UR2006_num)) %>%
  mutate(UR2006_fact = factor(UR2006_num))

num_samples = 1

# Fit the sequence type distribution using the island model
st_i = st_fit(formula = Source ~ ST,
         non_primary = "Human",
         data = dataset,
         method = "island",
         sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
         samples = num_samples)

# and using the Dirichlet
st_d = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = dataset,
              method = "dirichlet",
              iters = num_samples)

# do the attribution
at_if <- attribution(ST ~ UR2006_fact, st_i, manawatu %>% filter(Source == "Human"))
at_il <- attribution(ST ~ UR2006_num, st_i, manawatu %>% filter(Source == "Human"))
at_df <- attribution(ST ~ UR2006_fact, st_d, manawatu %>% filter(Source == "Human"))
at_dl <- attribution(ST ~ UR2006_num, st_d, manawatu %>% filter(Source == "Human"))

df_if <- predict(at_if, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_fact([-0-9]+)=1", convert=TRUE) %>%
  replace_na(replace=list(UR_num=-3)) %>% mutate(GenotypeModel = 'Island', AttributionModel='Categorical')
df_il <- predict(at_il, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE) %>%
  mutate(GenotypeModel = 'Island', AttributionModel='Linear')
df_df <- predict(at_df, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_fact([-0-9]+)=1", convert=TRUE) %>%
  replace_na(replace=list(UR_num=-3)) %>% mutate(GenotypeModel = 'Dirichlet', AttributionModel='Categorical')
df_dl <- predict(at_dl, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE) %>%
  mutate(GenotypeModel = 'Dirichlet', AttributionModel='Linear')

df <- bind_rows(df_if, df_il, df_df, df_dl)

save(list='df', file='attribution_fits.Rdata')
load('attribution_fits.Rdata')

plot_df <- df %>%
  group_by(GenotypeModel, AttributionModel, Source, UR2006_num) %>% summarize(
    m = mean(p),
    li = quantile(p, 0.1),
    ui = quantile(p, 0.9)
  ) %>% ungroup %>%
  mutate(Source = factor(Source, levels=rev(c("Other", "Water", "Ruminants", "Poultry"))))

# plot
cairo_pdf("attr_all.pdf", width=7, height=5)
ggplot(plot_df) +
  geom_ribbon(aes(x=UR_num, ymin=li, ymax=ui, fill=Source), alpha=0.3) +
  geom_line(aes(x=UR_num, y=m, col=Source), lwd=1) +
  facet_grid(GenotypeModel~AttributionModel) +
  scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Rural", "Urban"), expand=c(0,0)) +
  scale_y_continuous(name="Percentage of cases", labels=scales::percent_format(), limits=c(0,1), expand=c(0,0)) +
  scale_color_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue2", Other="plum4", Water="green4")) +
  scale_fill_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue2", Other="plum4", Water="green4")) +
  theme_bw() +
  theme(text = element_text(family="Times"),
        legend.position = c(0.99,0.89),
        legend.justification = "right",
        legend.margin=margin(0,0,0,0),
        legend.background = element_rect(fill = 'transparent'),
        axis.text.x = element_text(hjust=c(-0.1,1.1)),
        axis.text.y = element_text(vjust=c(-0.1,rep(0.5, 3), 1.1)))
dev.off()
