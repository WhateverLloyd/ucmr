#UCMR5 Data Processing

#load packages
library(tidyverse)

#set working directory
setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/UCMR5forOPAL")

#read in most recent release of UCMR-5 data (October 2024)
ucmr <-  read_delim(file = "./data/ucmr5_oct2024/UCMR5_All.txt", delim = "\t") 
##preview data
head(ucmr)
##number of tested PWSIDS
length(unique(ucmr$PWSID))
##list of tested contaminants
unique(ucmr$Contaminant)

#remove samples for lithium (lithium is the only non-PFAS contaminant measured in UCMR5)
ucmr <- ucmr %>% filter(Contaminant != "lithium")
##remaining number of tested PWSIDS
length(unique(ucmr$PWSID))

#format date
ucmr$DATE <- as.Date(as.character(ucmr$CollectionDate), format = "%m/%d/%Y")

##date range
min(ucmr$DATE); max(ucmr$DATE)

#create unique ID for each sample by pasting together PWSID, sample date, and sample ID
ucmr$UNIQUE_ID <- paste(ucmr$PWSID, gsub("-", "", ucmr$DATE), ucmr$SampleID, sep = "_")

##check to see if there are any duplicates
ucmr %>% group_by(UNIQUE_ID, Contaminant) %>% summarize(n=n()) %>% filter(n>1)
##TX1011249_20240626_1011249ASE3 had duplicate entries for 8 PFAS, however all were non-detects
#remove duplicate entries
ucmr <- distinct(ucmr)

#convert results from ug/L to ng/L and change non-detects to zeros
ucmr$RESULT_ngL <- ifelse(is.na(ucmr$AnalyticalResultValue), 0, ucmr$AnalyticalResultValue*1000)
#also convert MRLs from ug/L to ng/L
ucmr$MRL_ngL <- ucmr$MRL*1000

#reshape data from long to wide giving each contaminant column
ucmr1 <-  ucmr %>%
          pivot_wider(id_cols = c(UNIQUE_ID, PWSID, DATE), 
            names_from = Contaminant, 
            names_glue = "{Contaminant}_Concentration_ngL",
            values_from = RESULT_ngL, 
            values_fill = NA)

#relabel columns where PFAS name has special characters
names(ucmr1)[names(ucmr1) == "11Cl-PF3OUdS_Concentration_ngL"] <- "X11ClPF3OUdS_Concentration_ngL"
names(ucmr1)[names(ucmr1) == "4:2 FTS_Concentration_ngL"] <- "X42FTS_Concentration_ngL"
names(ucmr1)[names(ucmr1) == "6:2 FTS_Concentration_ngL"] <- "X62FTS_Concentration_ngL"
names(ucmr1)[names(ucmr1) == "8:2 FTS_Concentration_ngL"] <- "X82FTS_Concentration_ngL"
names(ucmr1)[names(ucmr1) == "9Cl-PF3ONS_Concentration_ngL"] <- "X9ClPF3ONS_Concentration_ngL"
names(ucmr1)[names(ucmr1) == "HFPO-DA_Concentration_ngL"] <- "HFPODA_Concentration_ngL"

#add in MRLs
mrls <- ucmr %>% 
        select(Contaminant, MRL_ngL) %>% 
        distinct() %>%
        pivot_wider(names_from = Contaminant, 
          names_glue = "{Contaminant}_MRL_ngL",
          values_from = MRL_ngL)

#relabel columns where PFAS name has special characters
names(mrls)[names(mrls) == "11Cl-PF3OUdS_MRL_ngL"] <- "X11ClPF3OUdS_MRL_ngL"
names(mrls)[names(mrls) == "4:2 FTS_MRL_ngL"] <- "X42FTS_MRL_ngL"
names(mrls)[names(mrls) == "6:2 FTS_MRL_ngL"] <- "X62FTS_MRL_ngL"
names(mrls)[names(mrls) == "8:2 FTS_MRL_ngL"] <- "X82FTS_MRL_ngL"
names(mrls)[names(mrls) == "9Cl-PF3ONS_MRL_ngL"] <- "X9ClPF3ONS_MRL_ngL"
names(mrls)[names(mrls) == "HFPO-DA_MRL_ngL"] <- "HFPODA_MRL_ngL"

#replicate rows of MRLs to match UCMR data frame
mrls_rep <- as.data.frame(lapply(mrls, rep, nrow(ucmr1)))

#append pfas MRLs to UCMR dataframe
ucmr1 <- cbind(ucmr1, mrls_rep)

#add in PWSID lat/lon data (coordinates from centroid of PWS service area boundary)
pwsxy <-  read.csv(file = "./outputs/UCMR5allPWSIDcoordinates.csv") %>%
          select(PWSID, LAT, LON)

ucmr1 <- left_join(ucmr1, pwsxy, by = "PWSID")

ucmr1 <-  ucmr1 %>%
          select(UNIQUE_ID, PWSID, LAT, LON, DATE,
            ADONA_Concentration_ngL, ADONA_MRL_ngL,
            HFPODA_Concentration_ngL, HFPODA_MRL_ngL,
            NEtFOSAA_Concentration_ngL, NEtFOSAA_MRL_ngL,
            NFDHA_Concentration_ngL, NFDHA_MRL_ngL,
            NMeFOSAA_Concentration_ngL, NMeFOSAA_MRL_ngL,
            PFBA_Concentration_ngL, PFBA_MRL_ngL,
            PFBS_Concentration_ngL, PFBS_MRL_ngL,
            PFDA_Concentration_ngL, PFDA_MRL_ngL,
            PFDoA_Concentration_ngL, PFDoA_MRL_ngL,
            PFEESA_Concentration_ngL, PFEESA_MRL_ngL,
            PFHpA_Concentration_ngL, PFHpA_MRL_ngL,
            PFHpS_Concentration_ngL, PFHpS_MRL_ngL,
            PFHxA_Concentration_ngL, PFHxA_MRL_ngL,
            PFHxS_Concentration_ngL, PFHxS_MRL_ngL,
            PFMBA_Concentration_ngL, PFMBA_MRL_ngL,
            PFMPA_Concentration_ngL, PFMPA_MRL_ngL,
            PFNA_Concentration_ngL, PFNA_MRL_ngL,
            PFOA_Concentration_ngL, PFOA_MRL_ngL,
            PFOS_Concentration_ngL, PFOS_MRL_ngL,
            PFPeA_Concentration_ngL, PFPeA_MRL_ngL,
            PFPeS_Concentration_ngL, PFPeS_MRL_ngL,
            PFTA_Concentration_ngL, PFTA_MRL_ngL,
            PFTrDA_Concentration_ngL, PFTrDA_MRL_ngL,
            PFUnA_Concentration_ngL, PFUnA_MRL_ngL,
            X11ClPF3OUdS_Concentration_ngL, X11ClPF3OUdS_MRL_ngL,
            X42FTS_Concentration_ngL, X42FTS_MRL_ngL,
            X62FTS_Concentration_ngL, X62FTS_MRL_ngL,
            X82FTS_Concentration_ngL, X82FTS_MRL_ngL,
            X9ClPF3ONS_Concentration_ngL, X9ClPF3ONS_MRL_ngL)

#export data
write.csv(ucmr1, file = "./outputs/UCMR5_PFASSamplingData_03052025.csv", row.names = F)

