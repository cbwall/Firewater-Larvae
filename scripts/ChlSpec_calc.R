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
pacman::p_load("here", "tidyverse", "ggplot2", "ggh4x")


## File names -------------------
foldername<-'ChlaSpec' # folder of the day
filename<-'chla_pyrolarv.csv' # data
sampleID<-'chla.map_pyrolarv.csv' # template of sample IDs
metadata_file<-"firewater larvae_metadata.csv" # file with slurry vols and SA

## What is the path length?
PL<-0.484 # pathlength NEW glass plate Silbiger lab plate with 200ul 


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
metadata<-read_csv(here("data", metadata_file))

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

# remove the X in the column name and add a 0 in front of single digits to keep the well name similar to how the instrument exports
AllData<-AllData %>%
  mutate(column = as.numeric(str_remove(AllData$column, pattern = "X")), # remove the X
         column = ifelse(column<10, gsub("(\\d)+", "0\\1", column),column), # only add a 0 in front of 1-9
         well = paste0(Rows,column)) # paste with row name

AllData.df <- AllData %>%
  left_join(sampleIDNames, by = "well") %>% # merge with sampleIDnames
  filter(!Sample.ID %in% "EMPTY") %>% # remove empty wells
  na.omit %>% #omit.NAs
  dplyr::mutate(is_blank = str_detect(tolower(coalesce(Sample.ID, "acetone.blank")), "blank"),
                sample.type = ifelse(str_detect(Sample.ID, "blank"), "blank", "samples"))

### blank correct
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

### Run chla analysis ### 
# chl from Jeffry and Humphreys
# show per well values (technical replicates)
AllData.stat.summary<-AllData.cor %>%
  mutate(chla = 11.43*((Abs_663_corr - Abs_750_corr) / PL) - 0.64*((Abs_630_corr - Abs_750_corr)/PL),
         chlc = 27.09*((Abs_630 - Abs_750_corr) / PL) - 3.63*((Abs_663_corr - Abs_750_corr)/PL),
         totalchl = chla+chlc,
         daterun = Sys.Date()) %>%
  group_by(Sample.ID) %>%
  left_join(metadata, by="Sample.ID")%>%
  filter(!is_blank) 
  
write.csv(AllData.stat.summary, "output/chl.calc_summary.csv")

# calculate for the normalized data
# show as grouped by sample.ID (mean of technical replicates)
AllData.normalized <- AllData.stat.summary %>%
  group_by(Sample.ID) %>%
  summarise(mean.totalchl.ug.ml = mean(totalchl), 
            sd = sd(totalchl), 
            CV=(sd(totalchl)/mean(totalchl)*100)) %>%
left_join(metadata) %>%
  mutate(number.of.larvae = 50,
         extract.vol.mL = 0.5) %>%
  mutate(totalchl.ug.larv = (mean.totalchl.ug.ml*extract.vol.mL/number.of.larvae))

write.csv(AllData.normalized, "output/chl.norm.csv")

###### fix structure
make.fac <- c("Sample.ID", "Treatment", "Level", "Rep")
AllData.normalized[make.fac] <- lapply(AllData.normalized[make.fac], factor)
str(AllData.normalized)

AllData.normalized$Treatment<- factor(AllData.normalized$Treatment, levels=c("control", "burned", "inorganic"))
AllData.normalized$Level<- factor(AllData.normalized$Level, levels=c("control", "low", "medium", "high"))


# plot
tot.chl.plot<- ggplot(AllData.normalized, aes(interaction(Treatment, Level, sep = "\n "), totalchl.ug.larv, fill=Treatment))+
  geom_jitter(width=0.1, aes(color=Treatment), alpha=0.7)+
  geom_boxplot(outlier.shape = NA, alpha=0.6)+
  scale_fill_manual(values=c("dodgerblue", "coral", "goldenrod")) +
  scale_color_manual(values=c("dodgerblue", "coral", "goldenrod")) +
  ylab(expression(paste("total chl a+c2 (\u00B5g larvae"^{-1},")")))+
  xlab("Treatment State")+ theme_bw() + guides(fill="none", color="none")

tot.chl.plot
dev.print(pdf, "figures/pyrolarv_chl.pdf", width = 7, height = 5)
dev.off()


