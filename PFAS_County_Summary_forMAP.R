library(tidyverse)

setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - PFAS + Cancer/UCMR_EA")

pfas <- read.csv(file = "./outputs/UCMR5_by_TRACT_031225.csv")
pop <- read.csv("./outputs/ACS5_UStracts2020.csv") %>% select(TRACT_ID, TOTAL_POP)
regs <- read.csv(file = "./data/StatesFIPS_CensusRegions.csv")

pfas$TRACT_ID <- as.character(str_pad(pfas$TRACT_ID, width = 11, side = "left", pad = "0"))
pfas$COUNTY_ID <- substr(pfas$TRACT_ID, 1, 5)
pop$TRACT_ID <- as.character(str_pad(pop$TRACT_ID, width = 11, side = "left", pad = "0"))
regs$FIPS <- as.character(str_pad(regs$FIPS, width = 2, side = "left", pad = "0"))

data1 <- left_join(pop, pfas, by = "TRACT_ID")

countypfas <- data1 %>%
              group_by(COUNTY_ID) %>%
              summarize(TRACTS_TOTAL = n_distinct(TRACT_ID, na.rm = T),
                POP_TOTAL = sum(TOTAL_POP, na.rm = T),
                TRACTS_TESTED = n_distinct(TRACT_ID[!is.na(PFAS_yn)], na.rm = T),
                POP_TESTED = sum(TOTAL_POP[!is.na(PFAS_yn)], na.rm = T),
                TRACTS_PFAS = n_distinct(TRACT_ID[PFAS_yn == 1], na.rm = T),
                POP_PFAS = sum(TOTAL_POP[PFAS_yn == 1], na.rm = T)) %>%
              mutate(PFASprctTotalTracts = TRACTS_PFAS/TRACTS_TOTAL*100,
                PFASprctTestedTracts = TRACTS_PFAS/TRACTS_TESTED*100,
                PFASprctTotalPop = POP_PFAS/POP_TOTAL*100,
                PFASprctTestedPop = POP_PFAS/POP_TESTED*100)

countypfas$FIPS <- substr(countypfas$COUNTY_ID, 1, 2)
countypfas1 <- left_join(countypfas, regs, by = "FIPS")
countypfas1 <- countypfas1 %>% filter(!is.na(COUNTY_ID))
write.csv(countypfas1, "./temp_ucmr5summary_countylevel_forGISmap.csv", row.names = F)


