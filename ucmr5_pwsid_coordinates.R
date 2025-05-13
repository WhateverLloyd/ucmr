#PWSID Centroid Crosswalks

#this script generates a dataframe with all 7,214 UCMR5-tested public water 
#systems(PWSIDs) with latitude and longitude for their centroid. Most PWSIDs
#have a shapefile of their service area boundary which made it easy to get 
#centroid coordinates. However, some PWSIDs did not have a shapefile and we
#used ZIP code data to estimate the location of these PWSIDs.


#load packages
library(tidyverse)

#set working directory
setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/UCMR5forOPAL")


#import most recent release of UCMR-5 data (October 2024)
ucmr <-  read_delim(file = "./data/ucmr5_oct2024/UCMR5_All.txt", delim = "\t") 
#remove samples for lithium (lithium is the only non-PFAS contaminant measured in UCMR5)
ucmr <- ucmr %>% filter(Contaminant != "lithium")


#create dataset of UCMR5-tested PWSID coordinates by combining data from EPA's
#shapefile of community water system service area boundaries, and other PWSID-ZIP code 
#crosswalks.

#import PWSID service area boundary shapefile centroid coordinates
pwsXY1 <- read.csv(file = "./data/PWSID_Centroids.csv") %>%
          select(PWSID, LAT, LON)

#365 UCMR tested PWSIDs are missing from this shapefile. For the remaning
#PWSIDs, we used the PWSID to ZIP code crosswalk provided in the UCMR5 data 
#release. We found this crosswalk to be less accurate than the shapefile, so we
#will only use this location data when shapefile is absent. Further, this 
#crosswalk provides ZIP codes which do not have spatial boundaries. We must 
#convert the ZIP codes to their respective ZIP code tabulated areas (ZCTAs) and
#use the centroid coordinates calculated for those ZCTAs.

#import PWS-ZIP crosswalk
ucmr_zip <- read_delim(file = "./data/ucmr5_oct2024/UCMR5_ZIPCodes.txt", delim = "\t") 
ucmr_zip$ZIPCODE <- as.character(str_pad(ucmr_zip$ZIPCODE, width = 5, side = "left", pad = "0"))

#import ZIP-ZCTA crosswalk
zip_zcta <- read.csv(file = "./data/zip_to_zcta_2018.csv")
zip_zcta$ZIPCODE <- as.character(str_pad(zip_zcta$ZIP_CODE, width = 5, side = "left", pad = "0"))
zip_zcta$ZCTA <- as.character(str_pad(zip_zcta$ZCTA, width = 5, side = "left", pad = "0"))
zip_zcta <- zip_zcta %>% select(ZIPCODE, ZCTA)

#join PWS-ZIP to ZIP-ZCTA
pws_zcta <- left_join(ucmr_zip, zip_zcta, by = "ZIPCODE") %>%
            select(PWSID, ZCTA) %>%
            distinct() %>%
            drop_na()

#import ZCTA coordinates
zctaXY <- read.csv(file = "./data/ZCTA_2020_Centroids.csv")
zctaXY$ZCTA <- as.character(str_pad(zctaXY$ZCTA, width = 5, side = "left", pad = "0"))

#add coordinates to new PWS-ZCTA crosswalk
pws_zctaXY <- inner_join(pws_zcta, zctaXY, by = "ZCTA")

#UCMR's PWS-ZIP crosswalk includes PWSIDs which map to multiple ZIPs/ZCTAs. 
#since we are only interested in a single coordinate per PWSID, take the
#mean lat/lon values of all ZCTAs associated with PWSIDs.
pwsXY2 <- pws_zctaXY %>%
          group_by(PWSID) %>%
          summarize(LAT = mean(LAT, na.rm = T),
            LON = mean(LON, na.rm = T))

#15 UCMR-tested PWSIDs were missing from both the PWS shapefile and 
#UCMR-provided crosswalk. I searched these PWSIDs in EPA's database
#(https://enviro.epa.gov/envirofacts/sdwis/search) and created our own crosswalk
#joining these remaining PWSIDs to their respective ZCTAs
pwsmiss <- read.csv(file = "./data/UCMRunmatchedPWSIDSwZCTAS.csv")
pwsmiss$ZCTA <- as.character(str_pad(pwsmiss$ZCTA, width = 5, side = "left", pad = "0"))
pwsXY3 <-  left_join(pwsmiss, zctaXY, by = "ZCTA") %>%
            select(PWSID, LAT, LON)

#add "SOURCE" column to each set of PWSID and coordinates indicating the method
#used to get lat/lon. Then, join PWSID coordinates from different sources into 
#same dataframe for export.

#PWS shapefile centroid (only keep PWSIDS in UCMR5)
pwsXY1 <- pwsXY1 %>%
          filter(PWSID %in% ucmr$PWSID) %>%
          mutate(SOURCE = "PWS Shapefile Centroid")

#UCMR ZIP code crosswalk (only keep PWSIDS in UCMR5 that are missing from shapefile)
pwsXY2 <- pwsXY2 %>%
          filter(PWSID %in% ucmr$PWSID) %>%
          filter(!PWSID %in% pwsXY1$PWSID) %>%
          mutate(SOURCE = "UCMR ZIP Crosswalk")

#PWSID look up 
pwsXY3 <- pwsXY3 %>%
          filter(PWSID %in% ucmr$PWSID) %>%
          filter(!PWSID %in% pwsXY1$PWSID) %>%
          filter(!PWSID %in% pwsXY2$PWSID) %>%
          mutate(SOURCE = "SDWIDS Search")

#combine 3 PWSID coordinates into single dataframe
pwsXYfinal <- rbind(pwsXY1, pwsXY2, pwsXY3)

#check for duplicate PWSIDS
length(unique(pwsXYfinal$PWSID)) == nrow(pwsXYfinal)

#export final PWSID centroid dataset
write.csv(pwsXYfinal, "./outputs/UCMR5allPWSIDcoordinates.csv",row.names = F)