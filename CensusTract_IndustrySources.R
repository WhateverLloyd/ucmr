#load packages
library(readxl)
library(sf)
library(tidyverse)
library(janitor)

#set working directory
setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - PFAS + Cancer/UCMR_EA/")

#import 2020 census tract population centroids
trxy20 <- read_delim(file = "./data/CenPop2020_Mean_TR.txt", delim = ",")

trxy20 <- trxy20 %>%
          mutate(TRACT_ID = as.character(paste(
            str_pad(STATEFP, width = 2, side = "left", pad = "0"),
            str_pad(COUNTYFP, width = 3, side = "left", pad = "0"),
            str_pad(TRACTCE, width = 6, side = "left", pad = "0"),
            sep = ""))) %>%
          filter(!STATEFP %in% c("02","15","72")) %>%
          select(TRACT_ID, LONGITUDE, LATITUDE)

trxy20 <- st_as_sf(trxy20, coords = c("LONGITUDE","LATITUDE"), crs = 4326)

st_crs(trxy20)

trxy20p <- st_transform(trxy20, crs = 5070)

#import PFAS industry sources dataset from PFAS analytics tool
indsrc <- read_xlsx(path = "./data/IndustrySourcesPFAS_110424.xlsx")

indsrc <- indsrc %>%
          select(Industry, Longitude, Latitude) %>%
          na.omit()


indsrc  <- st_as_sf(indsrc , coords = c("Longitude","Latitude"), crs = 4326)

indsrcp <- st_transform(indsrc, crs = 5070)

#draw circular buffers around census tract population centroids
tr01km <- st_buffer(trxy20p, dist = 1000)
tr05km <- st_buffer(trxy20p, dist = 5000)

#intersect with industry source points

tr1src <- st_intersects(tr01km, indsrcp)
trxy <- st_drop_geometry(trxy20p) %>% rowid_to_column(var = "TRACT_rowid")
inds <- st_drop_geometry(indsrcp) %>% rowid_to_column(var = "SRC_rowid")

tr1src_df <- as.data.frame(tr1src)
names(tr1src_df) <- c("TRACT_rowid","SRC_rowid")
tr1src_df <- left_join(tr1src_df, trxy, by = "TRACT_rowid")
tr1src_df <- left_join(tr1src_df, inds, by = "SRC_rowid")

tr1src_out <- tr1src_df %>% 
              group_by(TRACT_ID, Industry) %>% 
              summarize(Count = n()) %>% 
              pivot_wider(names_from = Industry,
                values_from = Count, 
                values_fill = 0)

tr5src <- st_intersects(tr05km, indsrcp)
tr5src_df <- as.data.frame(tr5src)
names(tr5src_df) <- c("TRACT_rowid","SRC_rowid")
tr5src_df <- left_join(tr5src_df, trxy, by = "TRACT_rowid")
tr5src_df <- left_join(tr5src_df, inds, by = "SRC_rowid")

tr5src_out <- tr5src_df %>% 
  group_by(TRACT_ID, Industry) %>% 
  summarize(Count = n()) %>% 
  pivot_wider(names_from = Industry,
              values_from = Count, 
              values_fill = 0)


tr1src_out <- clean_names(tr1src_out, case = "all_caps")
colnames(tr1src_out)[-1] <- paste(colnames(tr1src_out)[-1], "1kmBuff", sep = "_")

tr5src_out <- clean_names(tr5src_out, case = "all_caps")
colnames(tr5src_out)[-1] <- paste(colnames(tr5src_out)[-1], "5kmBuff", sep = "_")

tr1src_out$ALL_1kmBuff <- rowSums(tr1src_out[-1])
tr5src_out$ALL_5kmBuff <- rowSums(tr5src_out[-1])

final <- left_join(trxy, tr1src_out, by = "TRACT_ID") %>%
          left_join(., tr5src_out, by = "TRACT_ID")

final[is.na(final)] <- 0

final <- final %>% select(!TRACT_rowid)


write.csv(final, "./outputs/IndustrySourceCounts_CensusTracts_1and5kmBuff.csv",row.names = F)
