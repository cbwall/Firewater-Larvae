########################################################################
### This calculates the Chlorophyll from a 96 well plate in corals from the SpectraMax iD5
### 
### Created by Dr. Nyssa Silbiger, edit by Chris Wall
### Edited on 03/05/2026
#########################################################################

## create a folder for your day of sampling (foldername below) and put your plate.csv files in there. All your data will be exported to that folder

#libraries

# load packages and use pacman to load all the packages you are missing!
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library
pacman::p_load("here", "tidyverse")


## File names -------------------
foldername<-'ChlaSpec' # folder of the day
filename<-'chla_raw.csv' # data
sampleID<-'chla.map.csv' # template of sample IDs
platename<-'oki_chla_plate4' # this will be the name of your file
metadata_file<-"Chl_metadata_plate4.csv" # file with slurry vols and SA

## What is the path length?
PL<-0.57 # pathlength for Silbiger lab plate with 200ul 


# DONT CHANGE ANYTHING BELOW HERE ----------------------------------

# Read in the Sample IDs
sampleIDNames<-read_csv(here('data',foldername,sampleID))

# Extract the letter and numbers from well, using sub, then use sprintf to format the output with leading zeros for numbers
letters <- sub("([A-Z]+)([0-9]+)", "\\1", sampleIDNames$well)
numbers <- as.numeric(sub("([A-Z]+)([0-9]+)", "\\2", sampleIDNames$well))
sampleIDNames$well <- sprintf("%s%02d", letters, numbers)

# Result: "A01" "A02" "A10" "A11"
sampleIDNames<-sampleIDNames %>% na.omit %>%
  select(well, Sample.ID) # pull out the location and the IDs

# read in metadata
metadata<-read_csv(here("data", foldername, metadata_file))

### format the 96 well plate data ################
#read in the rows w/o the dye
ChlPlate<-read.csv(here('data',foldername,filename), row.names = c("A","B","C","D","E","F","G","H"),nrows = 8, fileEncoding="latin1", skip = 17)
ChlPlate<-ChlPlate[, -c(1:2)]
# make the rownames a column
ChlPlate$Rows<-rownames(ChlPlate)

## Pull out each of the wavelengths to make own column
#630 wavelength
Abs_630<-ChlPlate %>%
  select(Rows,X1:X12)%>%
  pivot_longer(cols = X1:X12, values_to =  "Abs_630", names_to = "column")

#663 wavelength
Abs_663<-ChlPlate %>%
  select(Rows,X1.1:X12.1)%>% # pull out the 663 columns
  pivot_longer(cols = X1.1:X12.1, values_to =  "Abs_663", names_to = "column") %>%
  separate(col = column, extra = "drop", into = 'column') # remove the .1

#750 wavelength
Abs_750<-ChlPlate %>%
  select(Rows,X1.2:X12.2)%>% # pull out the 750 columns
  pivot_longer(cols = X1.2:X12.2, values_to =  "Abs_750", names_to = "column") %>%
  separate(col = column, extra = "drop", into = 'column') # remove the .2


## bring everything to one dataframe (reduce allows me to use the left_join function multiple times at once)
AllData<-Reduce(function(...) left_join(...), list(Abs_630,Abs_663,Abs_750))
#AllData$column<-gsub("X", "", AllData$column) # get rid of the X

# remove the X in the column name and add a 0 in front of single digits to keep the well name similar to how the instrument exports
AllData<-AllData %>%
  mutate(column = as.numeric(str_remove(AllData$column, pattern = "X")), # remove the X
  column = ifelse(column<10, gsub("(\\d)+", "0\\1", column),column), # only add a 0 in front of 1-9
  well = paste0(Rows,column)) # paste with row name

AllData.df <- AllData %>%
  left_join(sampleIDNames, by = "well") %>%
  na.omit %>%
  dplyr::mutate(is_blank = str_detect(tolower(coalesce(Sample.ID, "")), "blank"),
                sample.type = ifelse(str_detect(Sample.ID, "blank"), "blank", "samples"))

### Run chla analysis ###
###MAYA UPDATE THIS TO USE CORRECTED VALUES!!!
# chl from Jeffry and Humphreys
ChlData_raw<-AllData.df %>%
  mutate(chla = 11.43*((Abs_663 - Abs_750) / PL) - 0.64*((Abs_630 - Abs_750)/PL),
         chlc = 27.09*((Abs_630 - Abs_750) / PL) - 3.63*((Abs_663 - Abs_750)/PL),
         Totalchl = chla+chlc,
           daterun = Sys.Date()) %>%
  group_by(Sample.ID) %>%
  summarise(mean.Totalchl = mean(Totalchl), sd = sd(Totalchl), CV=(sd(Totalchl)/mean(Totalchl)*100))

write.csv(ChlData_raw, "output/chla.calc_indiv.csv")

## below to merge with sample data

  left_join(sampleIDNames) %>% # join with the sampleIDs
  drop_na(Sample.ID) %>%
  left_join(metadata) %>%
  mutate(chla_ug_cm2 = (chla*Slurry_vol_ml)/SA_cm2,  # normalize to coral slurry and SA
         chlc_ug_cm2 = (chlc*Slurry_vol_ml)/SA_cm2,
         Totalchl_ug_cm2 = (Totalchl*Slurry_vol_ml)/SA_cm2) %>%
  filter(!is_blank)
         

# averaged data by sample ID
ChlData_ave<-ChlData_raw %>%
  group_by(Sample.Name) %>%
  summarise_at(vars(chla:Totalchl, chla_ug_cm2:Totalchl_ug_cm2), .funs = list(function(x){mean(x, na.rm = TRUE)},
                                                 function(x){sd(x, na.rm = TRUE)}/sqrt(n()))) %>%
  rename(chla_mean = chla_fn1, chlc_mean=chlc_fn1, TotalChl_mean = Totalchl_fn1,
         chla_SE = chla_fn2, chlc_SE = chlc_fn2, Totalchl_SE = Totalchl_fn2,
         chla_ug_cm2_mean = chla_ug_cm2_fn1, chlc_ug_cm2_mean = chlc_ug_cm2_fn1,
         Totalchl_ug_cm2_mean = Totalchl_ug_cm2_fn1, 
         chla_ug_cm2_SE = chla_ug_cm2_fn2, chlc_ug_cm2_SE = chlc_ug_cm2_fn2,
         Totalchl_ug_cm2_SE = Totalchl_ug_cm2_fn2) %>%
  select(Sample.Name,chla_mean,chla_SE,chlc_mean,chlc_SE,TotalChl_mean,Totalchl_SE,
         chla_ug_cm2_mean,chla_ug_cm2_SE, chlc_ug_cm2_mean,chlc_ug_cm2_SE,Totalchl_ug_cm2_mean,Totalchl_ug_cm2_SE ) # put in a better order
  
### Export the data
## all the info
write_csv(ChlData_ave, here("Data", foldername, paste0(platename,"_summary.csv"))) # print the summary data
write_csv(ChlData_raw, here("Data", foldername, paste0(platename,"_raw.csv"))) # print the raw data





######### NOT USED since 750 is internal corrections
blank_vals <- AllData.df %>%
  filter(is_blank) %>%
  summarise(
    blank_630 = mean(Abs_630, na.rm = TRUE),
    blank_663 = mean(Abs_663, na.rm = TRUE),
    blank_750 = mean(Abs_750, na.rm = TRUE)
  )

AllData.cor <- AllData.df %>%
  mutate(
    Abs_630_corr = Abs_630 - blank_vals$blank_630,
    Abs_663_corr = Abs_663 - blank_vals$blank_663,
    Abs_750_corr = Abs_750 - blank_vals$blank_750
  )

