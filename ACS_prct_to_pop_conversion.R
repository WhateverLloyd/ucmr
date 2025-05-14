library(tidyverse)

setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - PFAS + Cancer/UCMR_EA")

acs <- read.csv(file = "./outputs/ACS5_UStracts2020.csv")

acs2 <- acs %>%
        mutate(across(c(4:19,22:24,28,29,31:53), 
          list(~(TOTAL_POP*./100)),
          .names="{.col}_1")) %>%
        rename_with(~str_replace(.,'_prct_1', '_pop'))

write.csv(acs2, file = "./outputs/ACS5_UStracts2020wPOP.csv", row.names = F)
