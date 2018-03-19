library(islandR)
library(dplyr)
library(tidyr)
library(lubridate)
library(forcats)

manawatu <- read.csv("~/data/sa_report/data/extract_attribution.csv") %>%
  mutate(SampledDate = as.Date(SampledDate), Year = year(SampledDate)) %>%
  filter(Year >= 2005, Year <= 2014) %>%
  mutate(Source4 = fct_collapse(Source, Poultry = c("Supplier A", "Supplier B", "Supplier Other", "Otherpoultry"),
                                Ruminants = c("Cattle", "Sheep"),
                                Water = "Environmental water",
                                Other = c("Cat_dog_pet", "Duck_poultry", "Spent_hen", "Pig",
                                          "Turkey", "Water_bird_wild", "Wild_bird_other"))) %>%
  filter(Source4 != "Human" | !is.na(UR_num) ) %>%
  mutate(UR_fact = factor(UR_num), Intervention=factor(ifelse(Year < 2008, 'Before', 'After')))

num_samples <- 1000

# Fit the sequence type distribution using the island model
st_i = st_fit(formula = Source4 ~ ST,
         non_primary = "Human",
         data = manawatu,
         method="island",
         sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
         samples = num_samples)

# and using the Dirichlet
st_d = st_fit(formula = Source4 ~ ST,
              non_primary = "Human",
              data = manawatu,
              method="dirichlet",
              iters=num_samples)

# do the attribution
at_if <- attribution(ST ~ UR_num*Intervention, st_i, manawatu %>% filter(Source4 == "Human"))
at_id <- attribution(ST ~ UR_num*Intervention, st_d, manawatu %>% filter(Source4 == "Human"))

df_il <- predict(at_if, FUN=identity) %>% extract(X, into="UR_num", regex="UR_num=([-0-9]+)", convert=TRUE, remove=FALSE) %>%
  extract(X, into="Intervention", regex="Intervention([A-Za-z]+)=1") %>%
  replace_na(replace=list(Intervention="After")) %>%
  mutate(GenotypeModel = 'Island', AttributionModel='Linear')
df_dl <- predict(at_id, FUN=identity) %>% extract(X, into="UR_num", regex="UR_num=([-0-9]+)", convert=TRUE, remove=FALSE) %>%
  extract(X, into="Intervention", regex="Intervention([A-Za-z]+)=1") %>%
  replace_na(replace=list(Intervention="After")) %>%
  mutate(GenotypeModel = 'Dirichlet', AttributionModel='Linear')

df <- bind_rows(df_il, df_dl)

save(list='df', file='attribution_fits_intervention.Rdata')
load('attribution_fits_intervention.Rdata')

library(dplyr)

plot_df <- df %>%
  group_by(GenotypeModel, AttributionModel, Source, UR_num, Intervention) %>% summarize(
    m = mean(p),
    li = quantile(p, 0.1),
    ui = quantile(p, 0.9)
  ) %>% ungroup %>%
  mutate(Source = factor(Source, levels=rev(c("Other", "Water", "Ruminants", "Poultry"))),
         Intervention = factor(Intervention, levels=c("Before", "After"), labels=c("2005-2007", "2008-2014")))

# plot
library(ggplot2)

cairo_pdf("attr_intervention.pdf", width=7, height=5)
ggplot(plot_df) +
  geom_ribbon(aes(x=UR_num, ymin=li, ymax=ui, fill=Source), alpha=0.3) +
  geom_line(aes(x=UR_num, y=m, col=Source), lwd=1) +
  facet_grid(GenotypeModel~Intervention) +
  scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Highly rural", "Highly Urban"), expand=c(0,0)) +
  scale_y_continuous(name="Percentage of cases", labels=scales::percent_format(), limits=c(0,1), expand=c(0,0)) +
  scale_color_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue", Other="plum4", Water="green4")) +
  scale_fill_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue", Other="plum4", Water="green4")) +
  theme_bw() +
  theme(text = element_text(family="Times"),
        legend.position = c(0.99,0.89),
        legend.justification = "right",
        legend.margin=margin(0,0,0,0),
        legend.background = element_rect(fill = 'transparent'),
        axis.text.x = element_text(hjust=c(-0.1,1.1)),
        axis.text.y = element_text(vjust=c(-0.1,rep(0.5, 3), 1.1)))
dev.off()
