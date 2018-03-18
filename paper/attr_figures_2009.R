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
  filter(Source4 != "Human" | (!is.na(UR_num) & Year >= 2008) ) %>%
  mutate(UR_fact = factor(UR_num))

# Fit the sequence type distribution using the island model
st_i = st_fit(formula = Source4 ~ ST,
         non_primary = "Human",
         data = manawatu,
         method="island",
         sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
         samples = 1000)

# and using the Dirichlet
st_d = st_fit(formula = Source4 ~ ST,
              non_primary = "Human",
              data = manawatu,
              method="dirichlet",
              iters=1000)

# do the attribution
at_if <- attribution(ST ~ UR_fact, st_i, manawatu %>% filter(Source4 == "Human"))
at_il <- attribution(ST ~ UR_num, st_i, manawatu %>% filter(Source4 == "Human"))
at_df <- attribution(ST ~ UR_fact, st_d, manawatu %>% filter(Source4 == "Human"))
at_dl <- attribution(ST ~ UR_num, st_d, manawatu %>% filter(Source4 == "Human"))

df_if <- predict(at_if, FUN=identity) %>% extract(X, into="UR_num", regex="UR_fact([-0-9]+)=1", convert=TRUE) %>%
  replace_na(replace=list(UR_num=-3)) %>% mutate(GenotypeModel = 'Island', AttributionModel='Categorical')
df_il <- predict(at_il, FUN=identity) %>% extract(X, into="UR_num", regex="UR_num=([-0-9]+)", convert=TRUE) %>%
  mutate(GenotypeModel = 'Island', AttributionModel='Linear')
df_df <- predict(at_df, FUN=identity) %>% extract(X, into="UR_num", regex="UR_fact([-0-9]+)=1", convert=TRUE) %>%
  replace_na(replace=list(UR_num=-3)) %>% mutate(GenotypeModel = 'Dirichlet', AttributionModel='Categorical')
df_dl <- predict(at_dl, FUN=identity) %>% extract(X, into="UR_num", regex="UR_num=([-0-9]+)", convert=TRUE) %>%
  mutate(GenotypeModel = 'Dirichlet', AttributionModel='Linear')

df <- bind_rows(df_if, df_il, df_df, df_dl)

save(list='df', file='attribution_fits_2008.Rdata')
load('attribution_fits_2008.Rdata')

library(dplyr)

plot_df <- df %>%
  group_by(GenotypeModel, AttributionModel, Source, UR_num) %>% summarize(
    m = mean(p),
    li = quantile(p, 0.1),
    ui = quantile(p, 0.9)
  ) %>% ungroup %>%
  mutate(Source = factor(Source, levels=rev(c("Other", "Water", "Ruminants", "Poultry"))))

# plot
library(ggplot2)

cairo_pdf("attr_2009.pdf", width=7, height=5)
ggplot(plot_df) +
  geom_ribbon(aes(x=UR_num, ymin=li, ymax=ui, fill=Source), alpha=0.3) +
  geom_line(aes(x=UR_num, y=m, col=Source), lwd=1) +
  facet_grid(GenotypeModel~AttributionModel) +
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
