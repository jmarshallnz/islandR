# script to produce barplot figure...

library(meshblocknz)
library(dplyr)
library(tidyr)
library(sf)
library(lubridate)
library(forcats)

# 1. Load meshblocknz and compute population by rurality by year
#    Technically want this only for the meshblocks that we're actually sampling from
# load phu dataset
phu <- read_sf("~/code/epiclustR/shape_files/midcentral_phu.shp")

pop <- mb2013 %>% filter(DHB_name == "Mid Central") %>%
  filter(MB2006 %in% phu$MB06) %>% select(MB2013, MB2006, Pop2001, Pop2006, Pop2013, UR2006_num, UR2013_num) %>%
  gather(Year, Population, Pop2001:Pop2013) %>% extract(Year, into="Year", regex="([0-9]+)", convert=TRUE)

# OK, now fit a linear model per meshblock to the population
pops <- split(pop, pop$MB2013)

predict_popn <- function(dat) {
  m <- lm(Population ~ Year, data=dat)
  data.frame(Year=2005:2015, Population=predict(m, data.frame(Year=2005:2015)))
}

pop_interp <- dplyr::bind_rows(lapply(pops, predict_popn), .id='MB2013') %>% mutate(MB2013 = as.numeric(MB2013)) %>%
  left_join(mb2013 %>% select(MB2006, MB2013, UR2006_num) %>% unique)

# now rurality by year popn...

# 2. Once done, compute number of cases by rurality by year
manawatu <- read.csv("extract_attribution.csv") %>%
  mutate(SampledDate = as.Date(SampledDate), Year = year(SampledDate)) %>%
  filter(Year >= 2005, Year <= 2014) %>%
  mutate(Source4 = fct_collapse(Source, Poultry = c("Supplier A", "Supplier B", "Supplier Other", "Otherpoultry"),
                                Ruminants = c("Cattle", "Sheep"),
                                Water = "Environmental water",
                                Other = c("Cat_dog_pet", "Duck_poultry", "Spent_hen", "Pig",
                                          "Turkey", "Water_bird_wild", "Wild_bird_other"))) %>%
  filter(Source4 == "Human", !is.na(UR_num))

case_rates <- pop_interp %>% group_by(UR2006_num, Year) %>% summarize(Population=sum(Population)) %>%
  left_join(manawatu %>% group_by(UR2006_num = UR_num, Year) %>% summarize(Cases = n()))

rates <- case_rates %>% ungroup %>% mutate(Group = ifelse(UR2006_num > 0, "Urban", "Rural")) %>%
  group_by(Group, Year) %>% summarize(Population=sum(Population, na.rm=TRUE), Cases=sum(Cases, na.rm=TRUE)) %>%
  mutate(CaseRate = Cases/Population*100000) %>%
  mutate(CaseRate = ifelse(Year == 2005, CaseRate*1.33, CaseRate)) %>%
  ungroup() %>% filter(Year < 2015) %>%
  mutate(Year = factor(Year))


library(ggplot2)
cairo_pdf('case_rates.pdf', width=6, height=3.5)
ggplot(rates, aes(x=Year, y=CaseRate, fill=Group)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme(text = element_text(family="Times"),
        legend.position = c(0.07, 0.91),
        legend.margin=margin(0,0,0,0),
        legend.background = element_rect(fill = 'transparent')) +
  scale_fill_manual(NULL, values=c("grey30", "grey70")) +
  ylab("Cases per 100,000 population")
dev.off()



