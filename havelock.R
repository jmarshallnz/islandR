library(islandR)
library(dplyr)
library(readxl)
library(forcats)
library(pubmlst)
library(tidyr)

# read in data
manawatu_sources <- read.csv("extract_for_attribution.csv")

# remap the sources
manawatu_sources <- manawatu_sources %>%
  filter(Source != "Human") %>%
#  filter(Source != "Environmental water", Source != "Wild_bird") %>%
#  mutate(SampledDate = as.Date(as.character(SampledDate))) %>%
#  filter(SampledDate >= "2014-01-01") %>%
  mutate(Source = as.character(fct_collapse(Source,
                               Poultry=c("Supplier A", "Supplier B",
                                         "Supplier Other", "Spent_hen",
                                         "Turkey", "Duck_poultry"),
                               Ruminants=c("Sheep", "Cattle"),
                               Other=c("Cat_dog_pet", "Pig"),
                               Wild_bird = c("Water_bird_wild", "Wild_bird_other"))))

# read in human data
havelock <- read_excel("For Massey 21 Aug2016a.xlsx",sheet=2)
cases <- havelock %>%
  select(Source=Group, ASP=aspA, GLN=glnA, GLT=gltA, GLY=glyA, TKT=tkt, UNC=uncA, PGM=pgm) %>%
  filter(Source == 'Clinical') %>%
  data.frame %>% impute_mlst_in_data %>%
  filter(!is.na(ST)) %>% select(-CC, -Coli)
cases$Source = "Human"

#havelock <- read_excel("~/Downloads/For Massey 21 Aug2016a.xlsx")[6:8]
#cases <- havelock %>%
#  rename(Cases = `Clinical Cases`, ST=`WGS-MLST`, Water=`Water supply`) %>%
#  filter(!is.na(ST)) %>%
#  mutate(ST = ifelse(substring(ST, 1, 3) == "ST7", "ST-7505", ST),
#         ST = as.numeric(substring(ST, 4, 10))) %>%
#  select(ST, Cases)

#cases <- cases[rep(1:nrow(cases), cases$Cases), 1] %>%
#  left_join(pubmlst) %>% mutate(Source = "Human") %>%
#  select(-Coli, -CC) %>% mutate(ST = as.integer(ST))

# expand out the rows
dat <- rbind(cases, manawatu_sources %>% select(one_of(names(cases))))
dat$Source = factor(dat$Source)
dat <- data.frame(dat)

# estimate sampling distribution on sources
st = st_fit(formula = Source ~ ST,
            non_primary = "Human",
            data = dat,
            method="island",
            sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC)

# attribute human cases
mod = attribution(ST ~ 1, st, data=subset(dat, Source == "Human"))

summary(mod)
quantiles <- predict(mod, newdata=NULL, FUN=quantile, c(0.025, 0.5, 0.975))
attribution <- quantiles %>%
  select(-X) %>% mutate(p = round(p*100, 1)) %>%
  spread(quantile, p)

restrict <- quantiles %>% spread(quantile, p) %>% rename(lc = `2.5%`, uc=`97.5%`)
posterior = predict(mod, FUN=identity)

# filter out middle 95%
filtered <- posterior %>%
  left_join(restrict) %>%
  filter(p > lc & p < uc)

pdf("havelock_attribution.pdf", width=8, height=5)
par(mar=c(4,6,4,2))
boxplot(p*100 ~ Source, data=filtered, main="Attribution of human cases", horizontal=TRUE, las=1, xlab="P(source) (%)", pch='.')
dev.off()

# experiment with ggplot2
ggplot(filtered) +
  geom_violin(aes(Source, p), fill="slate grey") + coord_flip()

st42 = t(st$sampling_dist[st$types == 42,,])
st1517 = t(st$sampling_dist[st$types == 1517,,])
st42 = st42 / rowSums(st42)
st1517 = st1517 / rowSums(st1517)

pdf("st42_1517_13sources.pdf", width=8, height=12)
par(mfrow=c(2,1), mar=c(4,6,4,2))
boxplot(st42, xlab="P(source|ST=42)", horizontal=TRUE, las=1, main="ST-42")
boxplot(st1517, xlab="P(source|ST=1517)", horizontal=TRUE, las=1, main="ST-1517")
dev.off()

st7505 = t(st$sampling_dist[st$types == 10001,,])
st7505 = st7505 / rowSums(st7505)
tab <- round(t(apply(st7505, 2, quantile, c(0.025, 0.5, 0.975)))*100, 1)

# summary table
st42tab <- round(t(apply(st42, 2, quantile, c(0.025, 0.5, 0.975)))*100, 1)
st1517tab <- round(t(apply(st1517, 2, quantile, c(0.025, 0.5, 0.975)))*100, 1)
write.csv(st42tab, "st42.csv")
write.csv(st1517tab, "st1517.csv")

