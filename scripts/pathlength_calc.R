########################################################################
### This calculates the pathlength (200ul) of water in a glass 96-well plate using the SpectraMax iD5
### 
### Created by Chris Wall from Silbiger/Powell script
### Edited on 11/02/2026
#########################################################################

# create a folder for your data

#libraries
# load packages and use pacman to load all the packages you are missing!
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library
pacman::p_load("here", "tidyverse", "ggplot2")


## File names ------------------- # can just read in, but you may use this format for Chla plates too
foldername<-'ChlaSpec' # folder of the day
filename<-'path length_full plate.csv' # data


### format the 96 well plate data ################
#read in the rows
pathPlate<-read.csv(here('data',foldername,filename), row.names = c("A","B","C","D","E","F","G","H"),nrows = 8, fileEncoding="latin1", skip = 2)
pathPlate<-pathPlate[, -c(1:2)]

# make the rownames a column and filter
pathPlate$Rows<-rownames(pathPlate)
pathPlate <- pathPlate %>% select(-"X.1")

## Pull out each of the wavelengths to make own column
#900 wavelength
Abs_900<-pathPlate %>%
  select(Rows,X1:X12)%>%
  pivot_longer(cols = X1:X12, values_to =  "Abs_900", names_to = "column")

#977 wavelength
Abs_977<-pathPlate %>%
  select(Rows,X1.1:X12.1)%>% # pull out the 663 columns
  pivot_longer(cols = X1.1:X12.1, values_to =  "Abs_977", names_to = "column") %>%
  separate(col = column, extra = "drop", into = 'column') # remove the .1

## bring everything to one dataframe (reduce allows me to use the left_join function multiple times at once)
path.dat<-Reduce(function(...) left_join(...), list(Abs_900,Abs_977))
#AllData$column<-gsub("X", "", AllData$column) # get rid of the X

# remove the X in the column name and add a 0 in front of single digits to keep the well name similar to how the instrument exports
path.dat.df<-path.dat %>%
  mutate(column = as.numeric(str_remove(AllData$column, pattern = "X")), # remove the X
         column = ifelse(column<10, gsub("(\\d)+", "0\\1", column),column), # only add a 0 in front of 1-9
         well = paste0(Rows,column)) # paste with row name

K.fac<- 0.18 # for dividing Abs 977 - 900.

# pathlength calculation, inspect the output see if anything is off here
Path.calc<-path.dat.df %>%
  mutate(Abs977.900 = (Abs_977 - Abs_900) / K.fac)

# visualize
Path.calc.plot<-ggplot(data = Path.calc, aes(x = Abs977.900)) +
  geom_histogram(fill="coral", alpha=0.6)+
  geom_vline(aes(xintercept = mean(Abs977.900)), color = "red", linetype = "dashed", size = 1)

# remove outlier
Path.calc.corrected<-Path.calc %>% filter(!Abs977.900 > 0.58)

# confirm revised plot
Path.calc.plot.revise<-ggplot(data = Path.calc.corrected, aes(x = Abs977.900)) +
  geom_histogram(fill="coral", alpha=0.6)+
  geom_vline(aes(xintercept = mean(Abs977.900)), color = "red", linetype = "dashed", size = 1)

# summarise
Path.calc.summary<- Path.calc.corrected %>%
  summarise(mean.path = mean(Abs977.900), sd = sd(Abs977.900), CV=(sd(Abs977.900)/mean(Abs977.900)*100))

#export
write.csv(Path.calc, "output/glass.plate.path.csv")
