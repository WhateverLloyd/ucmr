#PWS source data 
library(sf)
library(tidyverse)
library(readxl)

#set working directory
setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - PFAS + Cancer")

#IMPORT AND PROCESS UCMR5 DATA BY PWSID ----

#read in most recent release of UCMR-5 data (Jan 2025)
ucmr5 <-  read_delim(file = "./data_sources/ucmr5_jan2025/UCMR5_All.txt", delim = "\t") 

#calculate max and average PFAS values at each PWSID (regardless of PFAS species)
ucmr5a <- ucmr5 %>%
  #remove lithium samples (only non-PFAS contaminant in UCMR5)
  filter(Contaminant != "lithium") %>%
  #convert results from ug/L to ng/L and change non-detects (NAs) to 0
  mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
                             AnalyticalResultValue*1000)) %>%
  group_by(PWSID) %>%
  summarize(MAX_PFASngL = max(RESULT_ngL),
            AVG_PFASngL = mean(RESULT_ngL)) %>%
  mutate(PFAS_yn = ifelse(MAX_PFASngL>0,1,0))

#repeat but calculating max/avg of individual PFAS species at PWSID 
ucmr5b <- ucmr5 %>%
  filter(Contaminant != "lithium") %>%
  mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
                             AnalyticalResultValue*1000)) %>%
  group_by(PWSID, Contaminant) %>%
  summarize(max = max(RESULT_ngL)) %>%
  pivot_wider(names_from = Contaminant, 
              values_from = max, values_fill = NA)

#relabel columns where PFAS name has special characters
names(ucmr5b)[names(ucmr5b) == "11Cl-PF3OUdS"] <- "11Cl_PF3OUdS"
names(ucmr5b)[names(ucmr5b) == "4:2 FTS"] <- "4_2FTS"
names(ucmr5b)[names(ucmr5b) == "6:2 FTS"] <- "6_2FTS"
names(ucmr5b)[names(ucmr5b) == "8:2 FTS"] <- "8_2FTS"
names(ucmr5b)[names(ucmr5b) == "9Cl-PF3ONS"] <- "9Cl_PF3ONS"
names(ucmr5b)[names(ucmr5b) == "HFPO-DA"] <- "HFPO_DA"

#rerun PWSID summaries by PFAS groupings (eg, short vs. long)
#import PFAS taxonomy document
chem <- read_xlsx(path = "./data_sources/UCMR5_PFAS_Taxonomy.xlsx")

ucmr5grp <- left_join(ucmr5, chem, by = c("Contaminant" = "Abbreviation"))

#PFAS categories
ucmr5c <- ucmr5grp %>%
  filter(Contaminant != "lithium") %>%
  mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
                             AnalyticalResultValue*1000)) %>%
  group_by(PWSID, PFAS_CAT) %>%
  summarize(max = max(RESULT_ngL)) %>%
  pivot_wider(names_from = PFAS_CAT, 
              names_glue = "PFASgroup_{PFAS_CAT}",
              values_from = max, values_fill = NA)

#PFAS chain length
ucmr5d <- ucmr5grp %>%
  filter(Contaminant != "lithium") %>%
  mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
                             AnalyticalResultValue*1000)) %>%
  group_by(PWSID, Length) %>%
  summarize(max = max(RESULT_ngL)) %>%
  pivot_wider(names_from = Length,
              names_glue = "ChainLength_{Length}",
              values_from = max, values_fill = NA)


#calculate EPA's hazard index (HI) from UCMR5 data
#HI = (HFPO-DA/10) + (PFBS/2000) + (PFNA/10) + (PFHxS/10)
haz <-  ucmr5 %>%
  #keep only the 4 contaminants used in hazard index
  filter(Contaminant %in% c("HFPO-DA","PFBS","PFNA","PFHxS")) %>%
  #convert non-detects to 0's and results from ug/L (ppm) to ng/L (ppt)
  mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
                             AnalyticalResultValue*1000)) %>%
  #find the maximum concentration for each compound by PWSID
  group_by(PWSID, Contaminant) %>%
  summarize(MAX_ngL = max(RESULT_ngL, na.rm = T)) %>%
  #adjust concentrations following hazard index formula
  mutate(MAX_PFAS_adj = case_when(
    Contaminant == "HFPO-DA" ~ MAX_ngL/10,
    Contaminant == "PFBS" ~ MAX_ngL/2000,
    Contaminant == "PFNA" ~ MAX_ngL/10,
    Contaminant == "PFHxS" ~ MAX_ngL/10)) %>%
  #sum adjusted values to get HI and return number of HI compounds sampled
  group_by(PWSID) %>%
  summarize(PFAS_HI = sum(MAX_PFAS_adj, na.rm = T), n_spp = n()) %>%
  #return NA if HI is <1 and not all HI PFAS species were sampled 
  mutate(PFAS_HI = ifelse(n_spp == 4 | PFAS_HI >= 1, PFAS_HI, NA)) %>%
  select(-n_spp)

#Maximum Contaminant Levels (MCL)

#determine if PFAS are above or below MCL for each PWSID
mcl <-  ucmr5 %>%
  #filter results to just contaminants with MCL
  filter(Contaminant %in% c("PFOA", "PFOS", "PFHxS", "HFPO-DA", "PFNA")) %>%
  #convert non-detects to 0's and results from ug/L (ppm) to ng/L (ppt)
  mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
                             AnalyticalResultValue*1000)) %>%
  #calculate maximum result for each compound by PWSID
  group_by(PWSID, Contaminant) %>%
  summarize(MAX_ngL = max(RESULT_ngL, na.rm = T)) %>%
  #create variable indicating if concentration is above or below MCL
  mutate(PFAS_MCL = case_when(
    Contaminant %in% c("PFOA","PFOS") & MAX_ngL >= 4 ~ 1,
    Contaminant %in% c("PFHxS","HFPO-DA","PFNA") & MAX_ngL >= 10 ~ 1,
    .default = 0)) %>%
  #Indicate if PWSID has PFAS above its MCL. Also test to see if all 5 
  #species were tested and return NA if species were missing, unless one 
  #species was above MCL. 
  group_by(PWSID) %>%
  summarize(PFAS_MCL = max(PFAS_MCL), n_spp = n()) %>%
  mutate(PFAS_MCL = ifelse(n_spp == 5 | PFAS_MCL == 1, PFAS_MCL, NA)) %>%
  select(-n_spp)

#join all datasets together
ucmr5final <- left_join(ucmr5a, ucmr5b, by = "PWSID") %>%
  left_join(., ucmr5c, by = "PWSID") %>%
  left_join(., ucmr5d, by = "PWSID") %>%
  left_join(., haz, by = "PWSID") %>%
  left_join(., mcl, by = "PWSID")

#add binary variable indicating if PWSID was out of compliance (either HI >= 1 or
#any of the 5 regulated compounds exceeded their MCL)
ucmr5final <- ucmr5final %>%
  mutate(PFAS_EPA = case_when(is.na(PFAS_MCL) | is.na(PFAS_HI) ~ NA,
                              PFAS_MCL == 1 | PFAS_HI >= 1 ~ 1, .default = 0))

#ADD PWS-REPORTED SOURCE DATA FROM UCMR ----

#read in PWS additional info file
pwsinfo <- read_delim(file = "./data_sources/ucmr5_jan2025/UCMR5_AddtlDataElem.txt", delim = "\t")

#add variable for PWS-reported PFAS source
pwssrcA <-  pwsinfo %>% 
  filter(AdditionalDataElement == "PotentialPFASSources") %>%
  mutate(UCMRSOURCE_Any = ifelse(Response == "DK", NA, ifelse(Response == "Yes", 1, 0)))  %>%
  group_by(PWSID) %>%
  summarize(UCMRSOURCE_Any = max(UCMRSOURCE_Any, na.rm = T))

pwssrcA[pwssrcA==-Inf] <- NA
pwssrcA %>% group_by(UCMRSOURCE_Any) %>% summarize(n=n())
#5,815 PWS reported no, 773 reported yes, and 1,357 reported dont know

#add industry categories reported with PFAS source
pwssrcB <-  pwsinfo %>% 
  filter(AdditionalDataElement == "PotentialPFASSourcesDetail")  %>%
  mutate(FLAG = 1) %>%
  group_by(PWSID, Response) %>%
  summarize(FLAG = max(FLAG, na.rm = T))

pwssrcC <-  pwssrcB %>% 
  pivot_wider(id_cols = PWSID, 
              names_from = Response,
              names_glue = "UCMRSOURCE_{Response}",
              values_from = FLAG, 
              values_fill = 0)

pwssrcD <- full_join(pwssrcA, pwssrcC, by = "PWSID")

#replace NA's with zeros (unless PWS reported "don't know", in which case, keep NAs)
src_nona <- pwssrcD %>% filter(!is.na(UCMRSOURCE_Any)) 
src_nona[is.na(src_nona)]<-0
src_isna <- pwssrcD %>% filter(is.na(UCMRSOURCE_Any))
pwssrcE <- rbind(src_nona, src_isna)

#join to ucmr PFAS data
ucmr5final <- left_join(ucmr5final, pwssrcE, by = "PWSID") 

#ADD INDUSTRY SOURCE DATA FROM PFAS ANALYTICS TOOL ----

#import database of industry sources and create spatial file from lat/lon 
inds <- read_xlsx(path = "./UCMR_EA/data/IndustrySourcesPFAS_110424.xlsx")
inds <- inds %>% 
        filter(Industry != "-", Latitude != "-", Longitude != "-") %>%
        dplyr::select(Industry, Latitude, Longitude) %>% 
        distinct() %>%
        st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326)



#import PWS shapefile
pwssf <- st_read("./UCMR_EA/data/Water_System_Boundaries_012825/Final.shp")
pwssf <- st_transform(pwssf, crs = 4326)

#use 'st_join' to join PWS shapefile with intersecting industry sites
sf_use_s2(FALSE)
indpws <- st_join(pwssf, inds)
indpws <- st_drop_geometry(indpws) %>% as.data.frame()

#count number of industry sites by industry categorty for each PWSID
indpws1 <-  indpws %>% 
            group_by(PWSID, Industry) %>% 
            summarize(Count=n()) %>% 
            drop_na()

#create dataframe to replace industry category types with shorter labels
indlab <- data.frame(Industry = c("Airports (Part 139)","Airports", "Cement Mfg", 
              "Chemical Mfg", "Cleaning Product Mfg", "Consumer Products", 
              "National Defense", "Electronics Industry", "Fire Protection", 
              "Fire Training", "Furniture and Carpet", "Industrial Gas", 
              "Glass Products", "Metal Machinery Mfg", "Metal Coating", 
              "Mining and Refining", "Oil and Gas", "Paints and Coatings", 
              "Paper Mills and Products", "Petroleum", "Plastics and Resins", 
              "Printing", "Textiles and Leather", "Waste Management"),
            Label = c("AIP", "AIR", "CEM", "CHE", "CLE", "CON", "DEF", "ELE", 
              "FIR", "FIT", "FUR", "GAS", "GLA", "MET", "MTC", "MIN", "OIL", 
              "PAI", "PAP", "PET", "PLA", "PRI", "TEX", "WAM"))
  
indpws1 <- left_join(indpws1, indlab, by = "Industry")

#reshape data from long to wide giving each industry sector its own column with
#count data. also create new variable of total industry counts regardless of type
indpws2 <-  indpws1 %>% 
            select(!Industry) %>%
            pivot_wider(names_from = Label, 
              names_glue = "PATSOURCE_{Label}",
              values_from = Count, 
              values_fill = 0) %>%
            mutate(PATSOURCE_Any = rowSums(across(PATSOURCE_AIR:PATSOURCE_GAS))) 

#join to ucmr5 results

ucmr5final <- left_join(ucmr5final, indpws2, by = "PWSID")

#replace NAs in pfas analytics tool source columns with zeros
narp <- function(x){ifelse(is.na(x),0,x)}
ucmr5final[72:96] <- lapply(ucmr5final[72:96], FUN = narp)

#COMBINE DATA WITH PWSID GIS OUTPTUS (coordiantes, area, population, etc.)
pwsgis <- read.csv(file = "./data_working/UCMR5PWSID_GISoutputs.csv")

ucmr5final <- left_join(pwsgis, ucmr5final, by = "PWSID")

#export final dataset 
write.csv(ucmr5final, file = "./data_working/UCMR5byPWSwSOURCES_Jan2025Update.csv",row.names = F)

#GET PWSID centroid coordiantes ----
# 
# ##Public Water Systems shapefile
# pwssf <- st_read("./UCMR_EA/data/Water_System_Boundaries_012825/Final.shp")
# pwssf <- st_transform(pwssf, crs = 4326)
# 
# #get PWSID coordinates for UCMR5 tested PWSID in PWS shapefile
# pwsxy1 <- pwssf %>% 
#   filter(PWSID %in% ucmr5final$PWSID) %>% 
#   st_centroid() %>%
#   mutate(LON = unlist(map(geometry,1)),
#          LAT = unlist(map(geometry,2))) %>%
#   st_drop_geometry() %>%
#   dplyr::select(PWSID, LAT, LON)
# 
# #ZCTA shapefile and coordinates for each ZCTA centroid
# zctasf <- st_read("../General - Medicare Cohort/CCHealth/Data/ZIP_ZCTA_CENSUS Data/ZCTA Spatial Files/tl_2010_us_zcta510/tl_2010_us_zcta510.shp")
# zctasf <- st_transform(zctasf, crs = 4326)
# zctasf$ZCTA <- as.character(str_pad(zctasf$ZCTA5CE10, width = 5, side = "left", pad = "0"))
# sf_use_s2(FALSE)
# zctaxy <- st_centroid(zctasf) %>% 
#           st_set_crs(st_crs(pwssf))
# zctaxy$ZCTA <- as.character(str_pad(zctasf$ZCTA5CE10, width = 5, side = "left", pad = "0"))
# zctaxy <- zctaxy %>%
#           mutate(LON = unlist(map(geometry,1)),
#             LAT = unlist(map(geometry,2))) %>%
#           st_drop_geometry() %>%
#           dplyr::select(ZCTA, LAT, LON)
# 
# #PWSID-ZIP/ZCTA Crosswalks
# 
# #EPA's UCMR5 PWS-to-ZIP Crosswalk
# ucmr5xw <- read_delim(file = "./data_sources/ucmr5_oct2024/UCMR5_ZIPCodes.txt", delim = "\t") 
# ucmr5xw$ZIPCODE <- as.character(str_pad(ucmr5xw$ZIPCODE, width = 5, side = "left", pad = "0"))
# 
# #PWS-ZCTA crosswalk for 10 PWSIDS missing from PWS shapefile and EPA UCMR crosswalks
# pwsxw <- read.csv(file = "./data_sources/UCMRunmatchedPWSIDSwZCTAS.csv")
# pwsxw$ZCTA <- as.character(str_pad(pwsxw$ZCTA, width = 5, side = "left", pad = "0"))
# pwsxw <- pwsxw %>% select(PWSID, ZCTA)
# 
# #ZIP to ZCTA crosswalk
# zxw <- read.csv(file = "../General - CRB Exposure/Census Geography/zip_to_zcta_2018.csv")
# zxw$ZCTA <- as.character(str_pad(zxw$ZCTA, width = 5, side = "left", pad = "0"))
# zxw$ZIPCODE <- as.character(str_pad(zxw$ZIP_CODE, width = 5, side = "left", pad = "0"))
# zxw <- zxw %>% select(ZIPCODE, ZCTA) 
# 
# #convert EPA crosswalks from ZIP code to ZCTA and remove PWSIDs present in shapefile
# ucmr5xw <-  left_join(ucmr5xw, zxw, by = "ZIPCODE") %>% 
#             filter(!PWSID %in% pwssf$PWSID) %>%
#             select(!ZIPCODE) %>%
#             drop_na(.) %>%
#             distinct(.)
# 
# #join EPA crosswalk and missing PWSID crosswalk
# pwsxy2 <- rbind(ucmr5xw, pwsxw)
# 
# #remove PWSIDS that weren't UCMR tested
# pwsxy2 <- pwsxy2  %>% 
#           filter(PWSID %in% ucmr5final$PWSID)
# 
# #join lat/lon coordinates of ZCTAs matched to PWSIDs
# pwsxy2 <- left_join(pwsxy2 , zctaxy, by = "ZCTA")
# 
# #get mean coordinates for PWSIDS which matched to multiple ZCTAs
# pwsxy3 <- pwsxy2  %>%
#           group_by(PWSID) %>%
#           summarize(LAT = mean(LAT), LON = mean(LON))
# 
# pwsxy4 <- rbind(pwsxy1,pwsxy3) 

