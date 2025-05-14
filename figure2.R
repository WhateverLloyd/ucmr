#load packages
library(patchwork)
library(tidyverse)

#set working directory
setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - PFAS + Cancer/UCMR_EA")

#UCMR5 Tract-level Results
ucmr<- read.csv(file = "./outputs/UCMR5_by_TRACT_031225.csv")
ucmr$TRACT_ID <- as.character(str_pad(ucmr$TRACT_ID, width = 11, side = "left", pad = "0"))

#Guide to states, regions, divisions, etc.
regs <- read.csv(file = "./data/StatesFIPS_CensusRegions.csv")
regs$FIPS <- as.character(str_pad(regs$FIPS, width = 2, side = "left", pad = "0"))

acs <- read.csv("./outputs/ACS5_UStracts2020.csv")
acs <- acs %>% select(TRACT_ID, TOTAL_POP)
acs$TRACT_ID <- as.character(str_pad(acs$TRACT_ID, width = 11, side = "left", pad = "0"))
acs$FIPS <- substr(acs$TRACT_ID, 1, 2)

pfas <- left_join(acs, ucmr, by = "TRACT_ID")
pfas <- left_join(pfas, regs, by = "FIPS")

pfas <- pfas %>% filter(!FIPS %in% c("02","15", "60","66","69","72","78"))


pfas1 <-  pfas %>% 
          select(TRACT_ID, Region, TOTAL_POP, PFAS_yn, ChainLength_Short, 
            ChainLength_Long, ADONA, HFPO_DA, NEtFOSAA, NFDHA, NMeFOSAA, PFBA, 
            PFBS, PFDA, PFDoA, PFEESA, PFHpA, PFHpS, PFHxA, PFHxS, PFMBA, PFMPA,
            PFNA, PFOA, PFOS, PFPeA, PFPeS, PFTA, PFTrDA, PFUnA)

pfas1[4:30] <- lapply(pfas1[4:30], 
                FUN = function(x){
                  ifelse(is.na(x), "NotTested", 
                    ifelse(x==0, "NotDetected", "Detected"))})

pfas2 <-  pfas1 %>%
          pivot_longer(!c(TRACT_ID, Region, TOTAL_POP), 
            names_to = "CHEM", values_to = "RESULT") 

chemsum <-  pfas2 %>% 
            filter(!CHEM %in% c("PFAS_yn","ChainLength_Long","ChainLength_Short")) %>%
            group_by(CHEM, RESULT) %>% 
            summarize(NTRACTS = n_distinct(TRACT_ID, na.rm = T))

chemsum <- chemsum %>% pivot_wider(names_from = RESULT, values_from = NTRACTS, values_fill = 0)
chemsum$DetectionRate <- chemsum$Detected/(chemsum$Detected+chemsum$NotDetected)*100

#from this list, the 8 most commonly detected PFAS species are:
#PFPeA, PFBA, PFHxA, PFBS, PFOS, PFOA, PFHxS, PFHpA

pfas3 <-  pfas2 %>% 
          filter(CHEM %in% c("PFAS_yn","ChainLength_Short","ChainLength_Long",
            "PFPeA","PFBA","PFHxA","PFBS","PFOS","PFOA","PFHxS","PFHpA"))

pfas3$TOTAL_POPm <- pfas3$TOTAL_POP/1000000
  
#rename variables for figure
pfas3 <-  pfas3 %>%
          mutate(CHEM = case_when(
              CHEM == "PFAS_yn" ~ "Any PFAS",
              CHEM == "ChainLength_Short" ~ "Any Short Chain PFAS",
              CHEM == "ChainLength_Long" ~ "Any Long Chain PFAS",
              .default = CHEM),
            RESULT = case_when(
              RESULT == "NotTested" ~ "Not Tested",
              RESULT == "NotDetected" ~ "Not Detected",
              .default = RESULT),
            LENGTH = ifelse(CHEM == "Any PFAS", "Any", 
              ifelse(CHEM %in% c("PFBA","PFBS","PFHpA","PFHxA","PFPeA", "Any Short Chain PFAS"),"Short", "Long")))

#reorder variables for figure
pfas3$CHEM <- factor(pfas3$CHEM, levels = c("PFHxS","PFOA","PFOS","Any Long Chain PFAS","PFBS","PFHpA","PFHxA","PFBA","PFPeA","Any Short Chain PFAS","Any PFAS"))
pfas3$RESULT <- factor(pfas3$RESULT, levels = c("Not Tested","Not Detected","Detected"))
pfas3$Region <- factor(pfas3$Region, levels = c("Midwest","Northeast","West","South"))


#duplicate dataset with "CONUS"
pfasC <- pfas3 %>% mutate(Region = "CONUS")

pfas4 <- rbind(pfas3, pfasC)

pfas4$CATEGORY <- paste(pfas4$RESULT, pfas4$LENGTH, sep = "-")
pfas4$CATEGORY <- factor(pfas4$CATEGORY, levels = c("Not Tested-Any","Not Detected-Any","Detected-Any",
                                                  "Not Tested-Long","Not Detected-Long","Detected-Long",
                                                  "Not Tested-Short","Not Detected-Short","Detected-Short"))


pfas4$Region <- factor(pfas4$Region, levels = c("West","South","Midwest","Northeast","CONUS"))

fcolA <-  pfas4 %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = CATEGORY)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "#f0f0f0","#bdbdbd","#636363",
    "#fee6ce","#fdae6b","#e6550d",
    "#deebf7","#9ecae1","#3182bd")) + 
  labs(y = "", x = "Population (millions)") +
  theme_classic() +
  facet_wrap(facets = vars(Region), ncol = 2, as.table = F) +
  theme(legend.position = "none") 
fcolA




f1a <-  ggplot(data = pfas3, aes(x=TOTAL_POPm, y=CHEM, fill = RESULT, color = LENGTH)) + 
          geom_bar(stat = "identity") +
          scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
          labs(y = "", x = "") +
          theme_classic() +
          theme(legend.position = "none") +
          facet_wrap(facets = vars(Region))



f1b <-  ggplot(data = pfas3, aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
          geom_bar(stat = "identity") +
          scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
          labs(y = "", x = "Population (millions)", fill = "UCMR-5 PFAS", subtitle = "CONUS") +
          theme_classic() +
          theme(legend.position = "bottom")


f1a


f2a <-  pfas3 %>% filter(CHEM %in% c("Any PFAS", "Any Long Chain PFAS", "Any Short Chain PFAS")) %>% 
        ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
        labs(y = "", x = "") +
        theme_classic() +
        theme(legend.position = "none") +
        facet_wrap(facets = vars(Region))

f2b <-  pfas3 %>% filter(CHEM %in% c("Any PFAS", "Any Long Chain PFAS", "Any Short Chain PFAS")) %>%   
        ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
        labs(y = "", x = "Population (millions)", fill = "UCMR-5 PFAS", subtitle = "CONUS") +
        theme_classic() +
        theme(legend.position = "bottom")


f2a/f2b




f3a <-  pfas3 %>% filter(CHEM %in% c("PFHpA","PFHxS","PFOA","PFOS","PFBS","PFHxA","PFBA","PFPeA")) %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(Region))

f3b <-  pfas3 %>% filter(CHEM %in% c("PFHpA","PFHxS","PFOA","PFOS","PFBS","PFHxA","PFBA","PFPeA")) %>%   
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "Population (millions)", fill = "UCMR-5 PFAS", subtitle = "CONUS") +
  theme_classic() +
  theme(legend.position = "bottom")


f3a/f3b


frs <- pfas4 %>% filter(LENGTH == "Short", Region != "CONUS") %>% 
      ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "Short Chain PFAS") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(Region))

fcs <- pfas4 %>% filter(LENGTH == "Short", Region == "CONUS") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(Region))

frl <- pfas4 %>% filter(LENGTH == "Long", Region != "CONUS") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "Long Chain PFAS") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(Region))

fcl <- pfas4 %>% filter(LENGTH == "Long", Region == "CONUS") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(Region))

fra <- pfas4 %>% filter(LENGTH == "Any", Region != "CONUS") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "Any PFAS") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(Region))

fca <- pfas4 %>% filter(LENGTH == "Any", Region == "CONUS") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(Region))

((frs/fcs)|(frl/fcl))/(fra/fca)

fc <- pfas4 %>% filter(Region == "CONUS") %>% 
      ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
      labs(y = "", x = "", subtitle = "CONUS") +
      theme_classic() +
      theme(legend.position = "none") +
      facet_wrap(facets = vars(LENGTH), scales = "free_y", ncol = 1)

fm <- pfas4 %>% filter(Region == "Midwest") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "Midwest") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(LENGTH), scales = "free_y", ncol = 1)

fn <- pfas4 %>% filter(Region == "Northeast") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "Northeast") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(LENGTH), scales = "free_y", ncol = 1)

fs <- pfas4 %>% filter(Region == "South") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "South") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(LENGTH), scales = "free_y", ncol = 1)

fw <- pfas4 %>% filter(Region == "West") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "West") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(LENGTH), scales = "free_y", ncol = 1)



fc/((fm|fn)/(fw|fs))



fc <- pfas4 %>% filter(Region == "CONUS") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "CONUS") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(LENGTH), scales = "free_y", ncol = 1)

fm <- pfas4 %>% filter(Region == "Midwest") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "Midwest") +
  theme_classic() +
  theme(legend.position = "none") 

fn <- pfas4 %>% filter(Region == "Northeast") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "Northeast") +
  theme_classic() +
  theme(legend.position = "none") 

fs <- pfas4 %>% filter(Region == "South") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "South") +
  theme_classic() +
  theme(legend.position = "none") 

fw <- pfas4 %>% filter(Region == "West") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = RESULT)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f0f0f0","#bdbdbd","#636363")) + 
  labs(y = "", x = "", subtitle = "West") +
  theme_classic() +
  theme(legend.position = "none") 



fc/((fm|fn)/(fw|fs))


fcolR <-  pfas4 %>% 
          filter(Region != "CONUS") %>% 
          ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = CATEGORY)) + 
          geom_bar(stat = "identity") +
          scale_fill_manual(values = c(
            "#f0f0f0","#bdbdbd","#636363",
            "#fee6ce","#fdae6b","#e6550d",
            "#deebf7","#9ecae1","#3182bd")) + 
          labs(y = "", x = "") +
          theme_classic() +
          facet_wrap(facets = vars(Region)) +
          theme(legend.position = "none") 

fcolR

fcolC <-  pfas4 %>% 
  filter(Region == "CONUS") %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = CATEGORY)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "#f0f0f0","#bdbdbd","#636363",
    "#fee6ce","#fdae6b","#e6550d",
    "#deebf7","#9ecae1","#3182bd")) + 
  labs(y = "", x = "") +
  theme_classic() +
  facet_wrap(facets = vars(Region)) +
  theme(legend.position = "none") 

fcolC

fcolR/fcolC

fcolA <-  pfas4 %>% 
  ggplot(data = ., aes(x=TOTAL_POPm, y=CHEM, fill = CATEGORY)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "#f0f0f0","#bdbdbd","#636363",
    "#fee6ce","#fdae6b","#e6550d",
    "#deebf7","#9ecae1","#3182bd")) + 
  labs(y = "", x = "Population (millions)") +
  theme_classic() +
  facet_wrap(facets = vars(Region), ncol = 2, scales = "free") +
  theme(legend.position = "none") 
fcolA






pfas4 %>% filter(Region == "CONUS") %>% group_by(CHEM, CATEGORY) %>% summarize(POP=sum(TOTAL_POP, na.rm = T)) %>% View()
