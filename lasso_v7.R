#set up workspace ----
library(glmnet)
library(ggpubr)
library(tidyverse)
library(patchwork)

setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - PFAS + Cancer/UCMR_EA")


#prepare datasets ----

#UCMR5 Tract-level Results
ucmr <- read.csv(file = "./outputs/UCMR5_by_TRACT_031225.csv")
##format 11-digit tract ID
ucmr$TRACT_ID <- as.character(str_pad(ucmr$TRACT_ID, width = 11, side = "left", pad = "0"))
##create new binary variables indicating if tract is within PWS receiving PFAS treatment
##or PWS with known PFAS source
ucmr$TREATMENT_yn <- ifelse(is.na(ucmr$PFAStreatment), NA, ifelse(ucmr$PFAStreatment == 0, 0, 1))
ucmr$SOURCE_yn <- ifelse(is.na(ucmr$KnownSource), NA, ifelse(ucmr$KnownSource == 0, 0, 1))

#American Community Survey
acs <- read.csv("./outputs/ACS5_UStracts2020.csv")
##format 11-digit tract ID
acs$TRACT_ID <- as.character(str_pad(acs$TRACT_ID, width = 11, side = "left", pad = "0"))
##select variables of interest
acs <-  acs %>% 
  select(TRACT_ID, TOTAL_POP, MEDIAN_AGE, UNDER5_prct, OVER65_prct, 
         OVER85_prct, FEMALE_prct, WHITE_prct, BLACK_prct, AMIND_prct, 
         ASIAN_prct, PACISL_prct, HISPLAT_prct, FOREIGN_BORN_prct, LNGISO_prct,
         MED_INCOME_usd, GINI, DISABILITY_prct, SNAP_prct, POVERTY_prct, 
         AVG_HOUSEHOLD_SIZE, MED_RENT_usd, MED_HOME_VAL_usd, RENTING_prct,
         HOME_AGE, NOHS_prct, BACHELORS_prct, UNEMPLOYMENT_prct, 
         INSURANCE_NONE_prct)

#fix poverty percentage (from population not in poverty to population in poverty)
acs$POVERTY_prct <- 100 - acs$POVERTY_prct

#NHGIS
nhg <- read.csv(file = "./data/nhgis0002_csv/nhgis0002_ds258_2020_tract.csv")
nhg$TRACT_ID <- as.character(str_pad(nhg$GEOCODE, width = 11, side = "left", pad = "0"))

nhg <- nhg %>% 
  rename("POP_TOTAL" = U7I001, 
         "POP_URBAN" = U7I002, 
         "POP_RURAL" = U7I003) %>%
  mutate(POP_DENS = POP_TOTAL/AREALAND, 
         POP_URBANprct = POP_URBAN/POP_TOTAL*100,
         POP_RURALprct = POP_RURAL/POP_TOTAL*100) %>%
  select(TRACT_ID, AREALAND, POP_TOTAL, POP_DENS, POP_URBAN, 
         POP_URBANprct, POP_RURAL, POP_RURALprct)

nhg[is.na(nhg)] <- 0

#join datasets
ucmrall <-  left_join(acs, ucmr, by = "TRACT_ID") %>%
  left_join(., nhg, by = "TRACT_ID")

#add state/region/division to data
##import state guide to census regions
regs <- read.csv(file = "./data/StatesFIPS_CensusRegions.csv")
##format 2-digit state FIPS codes
regs$FIPS <- as.character(str_pad(regs$FIPS, width = 2, side = "left", pad = "0"))
ucmrall$FIPS <- substr(ucmrall$TRACT_ID, start = 1, stop = 2)
##join datasets
ucmrall <- left_join(ucmrall, regs, by = "FIPS")
##keep only 48 states in CONUS plus DC
ucmrall <- ucmrall %>% filter(!FIPS %in% c("60","66","69","72","78","02","15")) 

#organize final dataset for lasso modelling
final1 <- ucmrall %>%
  select("TRACT_ID", "Postal", "Region", "Division", "EPA", "PFAS_yn", 
         "PFAS_EPA", "TREATMENT_yn","SOURCE_yn", "POP_DENS", 
         "POP_URBANprct", "POP_RURALprct", "MEDIAN_AGE", "UNDER5_prct",
         "OVER65_prct", "OVER85_prct", "FEMALE_prct", "WHITE_prct", 
         "BLACK_prct", "AMIND_prct", "ASIAN_prct", "PACISL_prct", 
         "HISPLAT_prct", "FOREIGN_BORN_prct", "LNGISO_prct", 
         "MED_INCOME_usd", "POVERTY_prct", "GINI", "MED_RENT_usd", 
         "MED_HOME_VAL_usd", "RENTING_prct", "HOME_AGE", 
         "AVG_HOUSEHOLD_SIZE", "UNEMPLOYMENT_prct", "DISABILITY_prct",
         "SNAP_prct", "NOHS_prct", "BACHELORS_prct", "INSURANCE_NONE_prct")


#run models----

#1. scale dataset
final1 <- final1 %>% mutate(across(c(10:39), ~ scale(.)[, 1]))

#2. stratify and scale datasets by region (CONUS, Midwest, Northeast, South, West)

##CONUS
finalC <- final1 %>%
  select(c(1,6,10:39)) %>%
  drop_na() %>%
  mutate(across(-c(TRACT_ID,PFAS_yn), 
                ~ scale(.)[,1])) %>%
  mutate(STUDYAREA = "CONUS")

##Midwest
finalM <- final1 %>%
  filter(Region == "Midwest") %>%
  select(c(1,6,10:39)) %>%
  drop_na() %>%
  mutate(across(-c(TRACT_ID,PFAS_yn), ~ scale(.)[,1])) %>%
  mutate(STUDYAREA = "Midwest")

##Northeast
finalN <- final1 %>%
  filter(Region == "Northeast") %>%
  select(c(1,6,10:39)) %>%
  drop_na() %>%
  mutate(across(-c(TRACT_ID,PFAS_yn), ~ scale(.)[,1])) %>%
  mutate(STUDYAREA = "Northeast")

##South
finalS <- final1 %>%
  filter(Region == "South") %>%
  select(c(1,6,10:39)) %>%
  drop_na() %>%
  mutate(across(-c(TRACT_ID,PFAS_yn), ~ scale(.)[,1])) %>%
  mutate(STUDYAREA = "South")

##West
finalW <- final1 %>%
  filter(Region == "West") %>%
  select(c(1,6,10:39)) %>%
  drop_na() %>%
  mutate(across(-c(TRACT_ID,PFAS_yn), ~ scale(.)[,1])) %>%
  mutate(STUDYAREA = "West")


#3. run models on stratified and scaled datasets

set.seed(seed = 1001)
rnorm(7,0,1)

#CONUS
data <- finalC

#create X and Y matrices for use in 'glmnet'
xmat <- as.matrix(data[,3:32])
ymat <- as.matrix(data[,2])

#use cross validation to select lambda (1 se)
cvmod <- cv.glmnet(x = xmat, y = ymat,  family = "binomial", alpha=1)
lmda <- cvmod$lambda.1se
#run lasso again with selected lambda value
mod <- glmnet(x = xmat, y = ymat, family = "binomial", alpha = 1, lambda = lmda)

#save model outputs
C_out <- stack(mod$beta[,])
C_out$outcome <- "PFASyn"
C_out$geography <- "CONUS"
C_out$lambda_val <- lmda
C_out$n_tracts <- nrow(data)

#Midwest
data <- finalM

xmat <- as.matrix(data[,3:32])
ymat <- as.matrix(data[,2])

cvmod <- cv.glmnet(x = xmat, y = ymat,  family = "binomial", alpha=1)
lmda <- cvmod$lambda.1se
mod <- glmnet(x = xmat, y = ymat, family = "binomial", alpha = 1, lambda = lmda)

M_out <- stack(mod$beta[,])
M_out$outcome <- "PFASyn"
M_out$geography <- "Midwest"
M_out$lambda_val <- lmda
M_out$n_tracts <- nrow(data)

#Northeast
data <- finalN

xmat <- as.matrix(data[,3:32])
ymat <- as.matrix(data[,2])

cvmod <- cv.glmnet(x = xmat, y = ymat,  family = "binomial", alpha=1)
lmda <- cvmod$lambda.1se
mod <- glmnet(x = xmat, y = ymat, family = "binomial", alpha = 1, lambda = lmda)

N_out <- stack(mod$beta[,])
N_out$outcome <- "PFASyn"
N_out$geography <- "Northeast"
N_out$lambda_val <- lmda
N_out$n_tracts <- nrow(data)

#South
data <- finalS

xmat <- as.matrix(data[,3:32])
ymat <- as.matrix(data[,2])

cvmod <- cv.glmnet(x = xmat, y = ymat,  family = "binomial", alpha=1)
lmda <- cvmod$lambda.1se
mod <- glmnet(x = xmat, y = ymat, family = "binomial", alpha = 1, lambda = lmda)

S_out <- stack(mod$beta[,])
S_out$outcome <- "PFASyn"
S_out$geography <- "South"
S_out$lambda_val <- lmda
S_out$n_tracts <- nrow(data)

#West
data <- finalW

xmat <- as.matrix(data[,3:32])
ymat <- as.matrix(data[,2])

cvmod <- cv.glmnet(x = xmat, y = ymat,  family = "binomial", alpha=1)
lmda <- cvmod$lambda.1se
mod <- glmnet(x = xmat, y = ymat, family = "binomial", alpha = 1, lambda = lmda)

W_out <- stack(mod$beta[,])
W_out$outcome <- "PFASyn"
W_out$geography <- "West"
W_out$lambda_val <- lmda
W_out$n_tracts <- nrow(data)

#join all outputs into single dataset
ALL_out <- rbind(C_out, M_out, N_out, S_out, W_out)

#plot results ---- 

#add absolute value and direction of coefficients
ALL_out$absval <- abs(ALL_out$values)
ALL_out$dir <- ifelse(ALL_out$values == 0,"Null", 
                      ifelse(ALL_out$values > 0, "Positive", "Negative"))


unique(ALL_out$ind)
varguide <- data.frame("ind" = c(
  "POP_DENS",
  "POP_URBANprct",
  "POP_RURALprct",
  "MEDIAN_AGE",        
  "UNDER5_prct",
  "OVER65_prct",
  "OVER85_prct",
  "FEMALE_prct",
  "WHITE_prct",
  "BLACK_prct",
  "AMIND_prct",
  "ASIAN_prct",
  "PACISL_prct",
  "HISPLAT_prct",
  "FOREIGN_BORN_prct",
  "LNGISO_prct",
  "MED_INCOME_usd",
  "POVERTY_prct",
  "GINI",
  "MED_RENT_usd",
  "MED_HOME_VAL_usd",
  "RENTING_prct",
  "HOME_AGE",
  "AVG_HOUSEHOLD_SIZE",
  "UNEMPLOYMENT_prct",
  "DISABILITY_prct",
  "SNAP_prct",
  "NOHS_prct",
  "BACHELORS_prct",
  "INSURANCE_NONE_prct"),
  "label" = c(
    "Population Density",
    "Urban Population",
    "Rural Population",
    "Median Age",        
    "Age Under 5",
    "Age Over 65",
    "Age Over 85",
    "Female",
    "White",
    "Black",
    "American Indian",
    "Asian",
    "Pacific Islander",
    "Hispanic or Latino",
    "Foreign Born",
    "Linguistic Isolation",
    "Median Income",
    "Poverty",
    "Gini Index",
    "Median Rent",
    "Median Home Value",
    "Renting",
    "Home Age",
    "Household Size",
    "Unemployment",
    "Disability",
    "SNAP",
    "No High School",
    "Bachelors or Higher",
    "Uninsured"))

ALL_out <- left_join(ALL_out, varguide, by = "ind")

#sort variables from largest to smallest


#filter data by outcome and coefficient value
dat1 <- ALL_out %>% filter(absval >= 0.10)

#sort variables from largest to smallest
dat1$label <- factor(dat1$label, levels = c("Median Income","Hispanic or Latino", "Black","Uninsured","White","Median Rent","Household Size","Urban Population","Asian","Renting","Foreign Born","Poverty","Home Age","Population Density","Median Home Value"))

ggballoonplot(data = dat1, x="geography",y="label",size="absval", fill = "dir") + 
  theme(legend.position = "right") + 
  scale_fill_manual(values = c("white","black")) +
  labs(subtitle="",size = "Coefficient Size", fill = "Coefficient Sign")
# 
# dat1 <- ALL_out %>% filter(outcome == "PFASyn", absval >= 0.10)
# 
# ggballoonplot(data = dat1, x="geography",y="label",size="absval", fill = "dir") + 
#   theme(legend.position = "right") + 
#   scale_fill_manual(values = c("white","black")) +
#   labs(subtitle="")
# 
# dat1 <- ALL_out %>% filter(outcome == "PFASyn", absval >= 0.10)
# 
# ggballoonplot(data = dat1, x="geography",y="label",size="absval", fill = "dir") + 
#   theme(legend.position = "right") + 
#   scale_fill_manual(values = c("white","black")) +
#   labs(subtitle="")
# 
# dat1 <- ALL_out %>% filter(outcome == "PFASyn", absval >= 0.10)
# 


