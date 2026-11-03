########################################################################
### This calculates the soluble protein from a 96 well plate in corals from the SpectraMax iD5
### 
### Created by Chris Wall
### modified from Dr.Maya Powell, Dr. Nyssa Silbiger, and Hollie Putnam scripts
###
### Edited on 02/17/2026
#########################################################################

# load packages
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library

# use pacman to load all the packages you are missing!
pacman::p_load("here", "tidyverse", "broom", "dplyr", "readr", "ggpmisc", "cowplot")

############################### 
### The script is formatted to extract data and perform standard curve for multiple plates all from a single spec read-out .csv
### It uses R projects for directories: you should have a folder names "ProtSpec" in your project folder for data and output

####### read in the spec data
plates<-read.csv("data/ProtSpec/prot_raw.csv", skip=17)
plates<-plates[, -c(1:2, 15:21)]
plates[]<-lapply(plates, as.numeric)
colnames(plates)<- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")

# plate 2, plates come in as "last is first" based on spec output
plate2<-plates[(1:8),]
row.names(plate2) = c("A","B","C","D","E","F","G","H")
plate2$rows<-as.character(row.names(plate2))

# plate 1
plate1<-plates[(13:20),] # 4 extra trash rows here and some at end too
row.names(plate1) = c("A","B","C","D","E","F","G","H")
plate1$rows<-as.character(row.names(plate1))

### read in sampleIDs, in long form. Should have "sample.ID" (standards or sample names), "run" number, alphanumeric well
sampleIDs<-read.csv("data/ProtSpec/prot_map.csv") #metadata for spec
sampleIDs<- sampleIDs %>%
  dplyr::mutate(sample.type = ifelse(str_detect(sample.ID, "Stand"), "std", "samples")) %>% 
  select(sample.ID, sample.type, run, well)

# separate runs
sampleIDs.pt1<- sampleIDs[(sampleIDs$run=="run1"),]
sampleIDs.pt2<-sampleIDs[(sampleIDs$run=="run2"),]

# make into long format data
protPlate1 <- plate1 %>%
  dplyr::select(rows,1:12) %>%
  pivot_longer(cols = -rows, values_to = "abs", names_to = "columns") %>%
  mutate(columns = str_pad(as.character(columns), width = 1), # if want a leading 0 then: width = 2, pad = "0"
                          plate="plate1", 
                          well = paste0(rows, columns)) %>%  #make "wells" with paste of row and column
  select(plate, well, abs)

protPlate2 <- plate2 %>%
  dplyr::select(rows,1:12) %>%
  pivot_longer(cols = -rows, values_to = "abs", names_to = "columns") %>%
  mutate(column = str_pad(as.character(columns), width = 1), # only add a 0 in front of 1-9
                          plate="plate2", 
                          well = paste0(rows, columns)) %>%  #make "wells" with paste of row and column
 select(plate, well, abs)

#####################################
## combine with sampleID
#####################################
protPlate1.IDs<-merge.data.frame(protPlate1, sampleIDs.pt1, by="well", all=TRUE)
protPlate2.IDs<-merge.data.frame(protPlate2, sampleIDs.pt2, by="well", all=TRUE)

AllData<-rbind(protPlate1.IDs, protPlate2.IDs) #combine
AllData<-na.omit(AllData) # remove NAs

# rearrane
AllData<- AllData %>% dplyr::select(sample.ID, sample.type, run, plate, well, abs)

#####################################
# now you have 2 cleaned plate outputs
#####################################


#####################################
# standard curve: we will fit via linear regression and polynomial
#####################################
standards <- tribble(
  ~std, ~BSA_ug.mL,
  "StandardA",        2000, 
  "StandardB",        1500,
  "StandardC",        1000,
  "StandardD",         750,
  "StandardE",         500,
  "StandardF",         250,
  "StandardG",         125,
  "StandardH",          25,
  "StandardI",           0
)

# standard for plate1
std_curve_p1 <- AllData %>%
  filter(run == "run1", sample.type == "std") %>%
  select(run, well, sample.ID, abs) %>%
  rename(std = sample.ID) %>%
  group_by(std) %>%
  summarise(mean_abs = mean(abs)) %>%  # calculate mean of standard triplicates
  mutate(mean_abs_cor = mean_abs - mean_abs[std == "StandardI"]) %>%   # subtract blank absorbance value from all
  left_join(standards)

# standard for plate2
std_curve_p2 <- AllData %>%
  filter(run == "run2", sample.type == "std") %>%
  select(run, well, sample.ID, abs) %>%
  rename(std = sample.ID) %>%
  group_by(std) %>%
  summarise(mean_abs = mean(abs)) %>%  # calculate mean of standard triplicates
  mutate(mean_abs_cor = mean_abs - mean_abs[std == "StandardI"]) %>%   # subtract blank absorbance value from all
  left_join(standards)

# extract your average blank value to subtract from your wells
blank.plate1.cor<-as.numeric(std_curve_p1[9,2]) #0.0707
blank.plate2.cor<-as.numeric(std_curve_p2[9,2]) #0.0714


#####################################
# Fit model for standard curve and plot it
#####################################

######### try with linear fits ######### 
######### plate 1
mod1 <- lm(BSA_ug.mL ~ mean_abs_cor, data = std_curve_p1)
fitted.p1 <- mod1 %>% broom::augment() #add columns to the original dataset such as predictions, residuals and cluster assignments

# extract slope and intercept of linear model
intercept.p1<-mod1$coefficients[1]
slope.p1<-mod1$coefficients[2]

std_curve_p1_plot <- std_curve_p1 %>%
  ggplot(aes(x = mean_abs_cor, y = BSA_ug.mL)) +
  geom_point(color = "red", size = 3) +
  geom_line(data = fitted.p1, aes(x = mean_abs_cor, y = .fitted)) +
  labs(title = "Standard curve plate1")

######### plate 2
mod2 <- lm(BSA_ug.mL ~ mean_abs_cor, data = std_curve_p2)
fitted.p2 <- mod2 %>% broom::augment()

# extract slope and intercept of linear model
intercept.p2<-mod2$coefficients[1]
slope.p2<-mod2$coefficients[2]

std_curve_p2_plot <- std_curve_p2 %>%
  ggplot(aes(x = mean_abs_cor, y = BSA_ug.mL)) +
  geom_point(color = "red", size = 3) +
  geom_line(data = fitted.p2, aes(x = mean_abs_cor, y = .fitted)) +
  labs(title = "Standard curve plate2")

################## ################## ################## 
######### try with 3rd order polynomial fits, advised since daga isn't wholly linear 
################## ################## ################## 
######### plate 1
std_poly_p1<-ggplot(data = std_curve_p1, aes(x = mean_abs_cor, y = BSA_ug.mL)) +
  stat_poly_line(formula = y ~ poly(x, 3, raw = TRUE)) +
  stat_poly_eq(formula = y ~ poly(x, 3, raw = TRUE), use_label("eq", "R2")) +
  geom_point(color = "red", size = 3)

######### plate 2
std_poly_p2<-ggplot(data = std_curve_p2, aes(x = mean_abs_cor, y = BSA_ug.mL)) +
  stat_poly_line(formula = y ~ poly(x, 3, raw = TRUE)) +
  stat_poly_eq(formula = y ~ poly(x, 3, raw = TRUE), use_label("eq", "R2")) +
  geom_point(color = "red", size = 3)

# Build the plot and extract the equation
std_poly_p1_build <- ggplot_build(std_poly_p1)
equation_label.p1 <- std_poly_p1_build$data[[2]]$label
print(equation_label.p1)

std_poly_p2_build <- ggplot_build(std_poly_p2)
equation_label.p2 <- std_poly_p2_build$data[[2]]$label
print(equation_label.p2)

# write or the equations, solve for y = protein (x = mean corrected absorbance)
# "x" will be your blank corrected absorbance
y1<- prot1 ~ 1.28 + 664*(x) + 213*(x)^2 - 40.9*(x)^3 # 3rd order polynomial for plate1
y2<- prot2~ -3.60 + 690*(x) + 112*(x)^2 - 3.10*(x)^3 # 3rd order polynomial for plate2

#####################################
#### subtract control on plate 1 and 2
#####################################
#just have samples, average across replicates
AllData.samples<-AllData %>% 
  dplyr::select(run, sample.type, sample.ID, abs) %>%
  filter(grepl("samples", sample.type)) %>%  # only sample names
  group_by(run, sample.ID) %>% # group by sample IDs
  summarise(mean_abs = mean(abs)) %>% # get means by sample
  dplyr::mutate(mean_abs_cor = if_else(run=="run1", #blank corrected
         mean_abs - blank.plate1.cor, 
         mean_abs - blank.plate2.cor),
         prot.ug.mL = if_else(run=="run1", #mutliply by slope for standard curve 
                              (1.28 + 664*(mean_abs_cor) + 213*(mean_abs_cor)^2 - 40.9*(mean_abs_cor)^3),  
                              (-3.60 + 690*(mean_abs_cor) + 112*(mean_abs_cor)^2 - 3.10*(mean_abs_cor)^3)))
  
# visually check that all points fall within the standard curve

##### plate 1
std_curve_p1_plot + 
  geom_point(data = AllData.samples[(AllData.samples$run=="run1"),], 
             aes(x = mean_abs_cor, y = prot.ug.mL), pch = "X", cex = 5, alpha = 0.3) +
  labs(title = "Plate 1: projected on standard curve")

std.data.plate1<- std_poly_p1 + 
  geom_point(data = AllData.samples[(AllData.samples$run=="run1"),], 
             aes(x = mean_abs_cor, y = prot.ug.mL), pch = "X", cex = 5, alpha = 0.3) +
    labs(title = "Plate 1: projected on std polynomial")

##### plate 2
std_curve_p2_plot + 
  geom_point(data = AllData.samples[(AllData.samples$run=="run2"),], 
             aes(x = mean_abs_cor, y = prot.ug.mL), pch = "X", cex = 5, alpha = 0.3) +
  labs(title = "Plate 2:  projected on standard curve")

std.data.plate2<- std_poly_p2 + 
  geom_point(data = AllData.samples[(AllData.samples$run=="run2"),], 
             aes(x = mean_abs_cor, y = prot.ug.mL), pch = "X", cex = 5, alpha = 0.3) +
  labs(title = "Plate 2: projected on std polynomial")
  
plot_grid(std.data.plate1, std.data.plate2)
dev.print(pdf, "data/ProtSpec/Std.plots.pdf", width=9, height=5)
dev.off()

##### final dataframe
prot<-AllData.samples

### Export the data
write_csv(prot, "output/prot_summary_p1.2.csv") # print the summary data

################ ################ ################ 
################ PAU!
################ ################ ################ 

