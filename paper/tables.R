library(islandR)
library(meshblocknz) # from https://github.com/jmarshallnz/meshblocknz
library(dplyr)
library(tidyr)
library(knitr)

# Table 1: Allelic profiles

st1 <- c(403, 2026, 474, 2343)
manawatu %>%
  filter(ST %in% st1) %>%
  select(ST, ASP:UNC) %>%
  unique %>%
  knitr::kable(format='latex', booktabs=TRUE, linesep='', row.names=FALSE)

# Table 2: Frequency of STs
st2 <- c(42, 45, 474, 2026, 2381)
manawatu %>%
  filter(ST %in% st2) %>%
  group_by(ST, Source) %>%
  summarize(Count=n()) %>%
  ungroup %>%
  spread(Source, Count, fill=0) %>%
  select(ST, Human, Poultry, Ruminants, Water, Others=Other) %>%
  knitr::kable(format='latex', booktabs=TRUE, linesep='')

# Table 3: Rurality information
cases <- manawatu %>%
  filter(Source == "Human", !is.na(UR2006_num)) %>%
  group_by(UR2006_num, UR2006_name) %>%
  summarize(`Human cases` = n())

# Our region is that covered by Mid Central Public Health Unit, which is the intersection
# of the Mid Central District Health Board and the Manawatu-Wanganui region:
popn <- mb2013 %>% filter(DHB_name == "Mid Central", RC2013_name == "Manawatu-Wanganui Region") %>%
  select(UR2006_num, Pop2006, Pop2013) %>%
  group_by(UR2006_num) %>% summarize_all(.funs=sum)

cases %>% left_join(popn) %>%
  rename(`Rurality scale`=UR2006_num, Description=UR2006_name,
         `2006` = Pop2006, `2013` = Pop2013) %>%
  knitr::kable('latex', booktabs=TRUE, linesep='', format.args=list(big.mark=','))

# Table 4: DIC values

dic_attr <- read.csv('dic_attribution.csv')
dic_attr %>% mutate(DIC=round(DIC, 1),
                    Model = factor(AttributionModel, levels=c("Linear", "Categorical"))) %>%
  select(Model, GenotypeModel, DIC) %>%
  spread(GenotypeModel, DIC) %>%
  knitr::kable('latex', booktabs=TRUE, linesep='', format.args=list(big.mark=','))

# Additional information for the text:

# Total number of cases by human/sources
manawatu %>%
  group_by(Type=ifelse(Source == "Human", "Human", "Source")) %>%
  summarize(Count=n())

# Total number of genotypes
manawatu %>% pull(ST) %>% unique %>% length

# Proportion of genotypes observed on humans
manawatu %>%
  group_by(Type=ifelse(Source == "Human", "Human", "Source")) %>%
  summarize(Count=length(unique(ST)))
125/348

# Number of cases without location information
manawatu %>%
  filter(Source == "Human") %>%
  group_by(Missing=is.na(UR2006_num)) %>%
  summarize(Count=n())
116/(116+1344)

# Frequency of ST-403 and ST-2343
manawatu %>%
  filter(ST %in% c(403, 2343)) %>%
  group_by(ST, Source) %>%
  summarize(Count=n()) %>%
  ungroup %>%
  spread(Source, Count, fill=0)
