#*housekeeping ----

setwd("C:/wd")
library(ggplot2)
library(tidyr)
library(dplyr)


#*Import 16S data----
#Import normalised, relative abundance data with all results (no "Othres" group)
Rel_abund_16S <- read.csv("/Normalisoidut_16S_maarat_analyysiin_rel_abund.csv",
                          header = TRUE)

#Function to replace "-" in tetracycline concentrations with "."
modify_strings <- function(df) {
  for (i in 1:nrow(df)) {
    if (startsWith(df$Samples_real[i], "t")) {
      df$Samples_real[i] <- sub("^(tc-0[^-]*)-(.*)", "\\1.\\2", df$Samples_real[i])
    }
  }
  return(df)
}

#Call the function
Rel_abund_16S2 <- modify_strings(Rel_abund_16S)



#Combine all data from same species
Rel_abund_16S2.2 <- Rel_abund_16S2 %>%
  group_by(Samples_real, Species, ab, concentration, temperature, replicate,) %>%
  summarise(normalisoidut = sum(normalisoidut), Rel_abund = sum(Rel_abund))


#*Import epicPCR data----
#Import epicPCR data with zero results included
epic_data_all <-
  read.csv2(
    "/Epicdata_with_zero_results.csv",
    header = T
  )

#Rearrange data to long form to make it easier to work with it in R
epic_all_long <- gather(epic_data_all, species, presence, Escherichia.coli:Microvirga.lotononidis...., factor_key=TRUE)

#Convert dot in species name to space
epic_all_long[,"species"]<-sub("[.]", " ", epic_all_long[,"species"])

#Function to remove extra dots in Microvirga species name
remove_dots <- function(df) {
  df$species <- gsub("\\.{4}$", "", df$species)
  return(df)
}

#Call function
epic_all_long <- remove_dots(epic_all_long)

#Create an empty column for sample name
epic_all_long[ , 'Samples_real'] = NA

#Combine treatment information to sample name
for (i in 1:nrow(epic_all_long)) {
  epic_all_long[i,]$Samples_real <- paste(epic_all_long[i,]$ab,
                                          epic_all_long[i,]$concentration,
                                          epic_all_long[i,]$temperature,
                                          epic_all_long[i,]$replicate,
                                          sep = "-")
}

#*Combine abundance data to epicdata----
#Create empty column
epic_all_long[ , 'Rel_abund'] = NA

#Capitalise column name to match between data sets so it can be used for reference
colnames(epic_all_long)[colnames(epic_all_long) == "species"] <- "Species"

epic_all_long2 <- epic_all_long

#Correct typo
epic_all_long2 <- epic_all_long2 %>%
  mutate(Species = ifelse(Species == 'Azospirillum brasiliense', 'Azospirillum brasilense', Species))


#Combine 16S data with epic data based on sample name and species name
# Initialize relative abundance column in data frame with 0
epic_all_long2$Rel_abund <- 0

epic_all_long3 <- epic_all_long2

#Loop through each row in data frame
for (i in 1:nrow(epic_all_long3)) {
  # Find the matching row in 16S data based on sample name and species
  match_row <- Rel_abund_16S2.2[Rel_abund_16S2.2$Samples_real == epic_all_long3$Samples_real[i] & 
                                  Rel_abund_16S2.2$Species == epic_all_long3$Species[i], ]
  
  # If a match is found, update epic-data with the relative abundance information
  if (nrow(match_row) > 0) {
    epic_all_long3$Rel_abund[i] <- match_row$Rel_abund
  }
}



#Add the species which weren't detected with epicPCR to the data frame along with their relative abundance data
epic_all_long4 <- epic_all_long3

#Match column names so data can be moved
colnames(Rel_abund_16S2.2)[colnames(Rel_abund_16S2.2) == "relative_abundance"] <- "Rel_abund"
colnames(Rel_abund_16S2.2)[colnames(Rel_abund_16S2.2) == "AB"] <- "ab"
colnames(Rel_abund_16S2.2)[colnames(Rel_abund_16S2.2) == "Concentration"] <- "concentration"
colnames(Rel_abund_16S2.2)[colnames(Rel_abund_16S2.2) == "Temperature"] <- "temperature"
colnames(Rel_abund_16S2.2)[colnames(Rel_abund_16S2.2) == "Replicate"] <- "replicate"


Rel_abund_16S3 <- subset(Rel_abund_16S2.2, select = -normalisoidut)


# Identify rows in in 16S data that are not in epic-data based on columns: species, sample name and relative abundance
new_rows <- anti_join(Rel_abund_16S3, epic_all_long4, by = c("Species", "Samples_real", "Rel_abund"))

# Add missing rows to epic-data with presence as 0 and sample number, origin and total number of species as 'NA'
if (nrow(new_rows) > 0) {
  new_rows$presence <- 0
  new_rows$sample.number.for.sequencing <- NA
  new_rows$origin <- NA
  new_rows$no_of_species <- NA
  epic_all_long4 <- rbind(epic_all_long4, new_rows)
}

epic_all_long4 <- epic_all_long4[!(epic_all_long4$Species %in% "NA NA"),]


epic_all_long4$presence <- as.factor(epic_all_long4$presence)


#Convert each origin to its own column
epic_all_long5 <- epic_all_long4 %>%
  mutate(pucL = ifelse(presence == 0, 0, as.integer(origin == "pucL")),
         pucD = ifelse(presence == 0, 0, as.integer(origin == "pucD")),
         intL = ifelse(presence == 0, 0, as.integer(origin == "intL")),
         intD = ifelse(presence == 0, 0, as.integer(origin == "intLD")),
         tnpL = ifelse(presence == 0, 0, as.integer(origin == "tnpL")),
         tnpD = ifelse(presence == 0, 0, as.integer(origin == "tnpD")))



#Combine the rows from same samples and species that were created when separating origins to their own columns
# Group by antibiotic, concentration, temperature, replicate and species
#and summarize to combine rows with the same values in other columns
epic_all_long6 <- epic_all_long5 %>%
  group_by(ab, concentration, temperature, replicate, Species,
           Samples_real, Rel_abund) %>%
  summarize(pucL = max(pucL),
            pucD = max(pucD),
            intL = max(intL),
            intD = max(intD),
            tnpL = max(tnpL),
            tnpD = max(tnpD)) %>%
  ungroup()


#How many species we have?
num_unique_values <- epic_all_long6 %>% 
  pull(Species) %>% 
  n_distinct

print(num_unique_values)
#Too many for plot!

#*Remove uninteresting species from plot----

#Remove species that don't carry the gene and are rare (rel.ab. < 0.02)

#Save a list of all species in data frame
unique_species <- unique(epic_all_long6$Species)
print(unique_species)

#Save a list of gene carrier species
epic_hosts <- unique(epic_all_long$Species)
print(epic_hosts)

#These are the species to keep
Species_keep <- c("Acidovorax defluvii", "Acinetobacter baumannii", "Aeromonas caviae", "Azospirillum brasilense",
                    "Citrobacter koseri", "Elizabethkingia meningoseptica", "Enterococcus faecalis",
                    "Escherichia coli", "Hafnia alvei", "Kluyvera intermedia", "Listeria innocua",
                    "Paracoccus denitrificans", "Pseudomonas chlororaphis", "Pseudomonas putida",
                    "Sphingobacterium multivorum", "Sphingobacterium spiritivorum", "Others")


#List of species that are rare and don't carry the gene
Species_extra <- unique_species[unique_species %in% epic_hosts | unique_species %in% Species_keep]
print(Species_extra)

#Create a new data frame by excluding rows with rare species 
epic_some_long <- subset(epic_all_long6, Species %in% Species_extra)


#There should be 19 species left
num_unique_values <- epic_some_long %>% 
  pull(Species) %>% 
  n_distinct()

#Check that this gives 19
print(num_unique_values)




#*Data for plot----
#Create a new column with presence-absence data for whether species was detected in 16S results
#(This is to create "not observed" bubbles)
epic_some_long2$Species_presence <- ifelse(epic_some_long2$Rel_abund > 0, 1, 0)

#Change it to a factor
epic_some_long2$Species_presence <- as.factor(epic_some_long2$Species_presence)

#Save unique species names in reverse order for creating y coordinates
lajien_laskeminen_kuvaajaan <- rev(unique(epic_some_long2$Species))

#Create y coordinates
epic_some_long3 <- epic_some_long2
epic_some_long3$y <- NA
epic_some_long3$y <- match(epic_some_long3$Species, lajien_laskeminen_kuvaajaan)
epic_some_long3 <- epic_some_long3[, !(names(epic_some_long3) %in% c("x"))]

#*Add missing samples that weren't observed in epic nor 16S data----
#There are still some data points missing that were absent in both 16S and epicPCR results
# Check for missing species names in each sample
all_combinations <- expand.grid(unique(epic_some_long3$Samples_real), unique(epic_some_long3$Species))
names(all_combinations)[names(all_combinations) == "Var1"] <- "Samples_real"
names(all_combinations)[names(all_combinations) == "Var2"] <- "Species"
existing_combinations <- as.data.frame(epic_some_long3[, c("Samples_real", "Species")])
missing_combinations <- setdiff(as.data.frame(all_combinations), existing_combinations)


# Create new rows for missing combinations
if (nrow(missing_combinations) > 0) {
  missing_rows <- data.frame(Samples_real = missing_combinations$Samples_real,
                             ab = sapply(missing_combinations$Samples_real,
                                         function(sample_name) epic_some_long3$ab[epic_some_long3$Samples_real == sample_name][1]),
                             replicate = sapply(missing_combinations$Samples_real,
                                                function(sample_name) epic_some_long3$replicate[epic_some_long3$Samples_real == sample_name][1]),
                             concentration = sapply(missing_combinations$Samples_real,
                                                    function(sample_name) epic_some_long3$concentration[epic_some_long3$Samples_real == sample_name][1]),
                             temperature = sapply(missing_combinations$Samples_real,
                                                  function(sample_name) epic_some_long3$temperature[epic_some_long3$Samples_real == sample_name][1]),
                             y = sapply(missing_combinations$Species,
                                        function(sample_name) epic_some_long3$y[epic_some_long3$Species == sample_name][1]),
                             Species_presence = rep(0, nrow(missing_combinations)),
                             Species = missing_combinations$Species,
                             Rel_abund = rep(0, nrow(missing_combinations)),
                             pucL = rep(0, nrow(missing_combinations)),
                             pucD = rep(0, nrow(missing_combinations)),
                             tnpL = rep(0, nrow(missing_combinations)),
                             intL = rep(0, nrow(missing_combinations)),
                             tnpD = rep(0, nrow(missing_combinations)),
                             intD = rep(0, nrow(missing_combinations))
  )
  
  # Add missing rows to the original data frame
  epic_some_long4 <- rbind(epic_some_long3, missing_rows)
}

#*Create plots----
#Divide data based on temperature to three plots
epic_some_t15 <- epic_some_long4[epic_some_long4$temperature == 15, ]
epic_some_t25 <- epic_some_long4[epic_some_long4$temperature == 25, ]
epic_some_t37 <- epic_some_long4[epic_some_long4$temperature == 37, ]

#Save distance for arrow from the center of the relative abundance marker
nuolen_pituus <- 0.35

#*Plot for 15 degree treatment results----
BUBBLEPLOT_15_grey <- ggplot(epic_some_t15, aes(x = replicate, y = Species,
                                                 shape = ifelse(Species_presence == 1, "22", "21"),
                                                 stroke = ifelse(Species_presence == 0, 1.5, 0),
                                                 colour = Species,
                                                 fill = ifelse(Species_presence == 1, Species, 0),
                                                 size = ifelse(Rel_abund > 0, Rel_abund, 0.15))) +
  facet_grid(. ~ concentration)+
  scale_size_continuous(breaks = c(0, 0.01, 0.1, 0.25, 0.5))+
  geom_point() +
  geom_text(label = "\u25B2", size = 3, colour = "red", inherit.aes = F,
            aes(x = replicate, y = ifelse(pucL == 1, y+nuolen_pituus, NA), family="Courier New"))+ #▲
  geom_text(label = "\u25BA", size = 3, colour = "black", inherit.aes = F,
            aes(x = replicate + nuolen_pituus, y = ifelse(tnpL == 1, y, NA), family = "Courier New"))+ # outer arrow or line of the right arrow ►
  geom_text(label = "\u25BA", size = 1.75, colour = "white", inherit.aes = F,
            aes(x = replicate + nuolen_pituus, y = ifelse(tnpL == 1, y, NA), family = "Courier New"))+ #inner arrow ►
  geom_text(label = "\u25BC", size = 3, colour = "blue", inherit.aes = F,
            aes(x = replicate, y = ifelse(intL == 1, y-nuolen_pituus, NA), family="Courier New"))+ #▼
  geom_text(label = "\u25C4", size = 3, colour = "#00c6ff", inherit.aes = F,
            aes(x = replicate - nuolen_pituus, y = ifelse(pucD == 1, y, NA), family="Courier New"))+ #◄
  scale_fill_manual(values=c("Acidovorax defluvii" = "#515151", "Acinetobacter baumannii" = "darkgrey", 
                             "Aeromonas caviae" = "#515151", "Agrobacterium tumefaciens" = "darkgrey", 
                             "Azospirillum brasilense" = "#515151", "Citrobacter koseri" = "darkgrey",
                             "Elizabethkingia meningoseptica" = "#515151", "Enterococcus faecalis" = "darkgrey",
                             "Escherichia coli" = "#59b0b5", "Hafnia alvei" = "darkgrey",
                             "Kluyvera intermedia" = "#515151", "Listeria innocua" = "darkgrey",
                             "Microvirga lotononidis" = "#515151", "Paracoccus denitrificans" = "darkgrey",
                             "Pseudomonas chlororaphis" = "#515151", "Pseudomonas putida" =  "darkgrey",
                             "Sphingobacterium multivorum" = "#515151", "Sphingobacterium spiritivorum" = "darkgrey",
                             "Tolumonas auensis" = "#515151", "0" = "white")) +
  scale_colour_manual(values=c("Acidovorax defluvii" = "#515151", "Acinetobacter baumannii" = "darkgrey", 
                               "Aeromonas caviae" = "#515151", "Agrobacterium tumefaciens" = "darkgrey", 
                               "Azospirillum brasilense" = "#515151", "Citrobacter koseri" = "darkgrey",
                               "Elizabethkingia meningoseptica" = "#515151", "Enterococcus faecalis" = "darkgrey",
                               "Escherichia coli" = "#59b0b5", "Hafnia alvei" = "darkgrey",
                               "Kluyvera intermedia" = "#515151", "Listeria innocua" = "darkgrey",
                               "Microvirga lotononidis" = "#515151", "Paracoccus denitrificans" = "darkgrey",
                               "Pseudomonas chlororaphis" = "#515151", "Pseudomonas putida" =  "darkgrey",
                               "Sphingobacterium multivorum" = "#515151", "Sphingobacterium spiritivorum" = "darkgrey",
                               "Tolumonas auensis" = "#515151", "0" = "white")) +
  scale_y_discrete(limits = rev(unique(epic_some_t15$Species)))+
  scale_shape_manual(values = c(21, 22)) +
  labs(x = "Replicate",
       size = "Relative\nabundance (16S)")+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "#ececec", fill = NA, linewidth = 2),
        axis.text.y=element_text(face = "italic"),
        panel.grid.minor.x = element_line(colour = "#ececec", linewidth = 0.4),
        strip.background = element_rect(fill = "#c1e1e9"),
        text = element_text(size = 15))+
  guides(colour = "none", fill = "none", shape = "none",
         size = guide_legend(override.aes = list(shape = 22, fill = "black")))

BUBBLEPLOT_15_grey



#*Plot for 25 degree treatment results----
BUBBLEPLOT_25_grey <- ggplot(epic_some_t25, aes(x = replicate, y = Species,
                                                 shape = ifelse(Species_presence == 1, "22", "21"),
                                                 colour = Species,
                                                 stroke = ifelse(Species_presence == 0, 1.5, 0),
                                                 fill = ifelse(Species_presence == 1, Species, 0),
                                                 size = ifelse(Rel_abund > 0, Rel_abund, 0.15))) +
  facet_grid(. ~ concentration)+
  scale_size_continuous(breaks = c(0, 0.01, 0.1, 0.25, 0.5))+
  geom_point() +
  geom_text(label = "\u25B2", size = 3, colour = "red", inherit.aes = F,
            aes(x = replicate, y = ifelse(pucL == 1, y+nuolen_pituus, NA), family="Courier New"))+ #▲
  geom_text(label = "\u25BA", size = 3, colour = "black", inherit.aes = F,
            aes(x = replicate + nuolen_pituus, y = ifelse(tnpL == 1, y, NA), family = "Courier New"))+ #►
  geom_text(label = "\u25BA", size = 1.75, colour = "white", inherit.aes = F,
            aes(x = replicate + nuolen_pituus, y = ifelse(tnpL == 1, y, NA), family = "Courier New"))+ #►
  geom_text(label = "\u25BC", size = 3, colour = "blue", inherit.aes = F,
            aes(x = replicate, y = ifelse(intL == 1, y-nuolen_pituus, NA), family="Courier New"))+ #▼
  geom_text(label = "\u25C4", size = 3, colour = "#00c6ff", inherit.aes = F,
            aes(x = replicate - nuolen_pituus, y = ifelse(pucD == 1, y, NA), family="Courier New"))+ #◄
  scale_fill_manual(values=c("Acidovorax defluvii" = "#515151", "Acinetobacter baumannii" = "darkgrey", 
                             "Aeromonas caviae" = "#515151", "Agrobacterium tumefaciens" = "darkgrey", 
                             "Azospirillum brasilense" = "#515151", "Citrobacter koseri" = "darkgrey",
                             "Elizabethkingia meningoseptica" = "#515151", "Enterococcus faecalis" = "darkgrey",
                             "Escherichia coli" = "#ffb400", "Hafnia alvei" = "darkgrey",
                             "Kluyvera intermedia" = "#515151", "Listeria innocua" = "darkgrey",
                             "Microvirga lotononidis" = "#515151", "Paracoccus denitrificans" = "darkgrey",
                             "Pseudomonas chlororaphis" = "#515151", "Pseudomonas putida" =  "darkgrey",
                             "Sphingobacterium multivorum" = "#515151", "Sphingobacterium spiritivorum" = "darkgrey",
                             "Tolumonas auensis" = "#515151", "0" = "white")) +
  scale_colour_manual(values=c("Acidovorax defluvii" = "#515151", "Acinetobacter baumannii" = "darkgrey", 
                               "Aeromonas caviae" = "#515151", "Agrobacterium tumefaciens" = "darkgrey", 
                               "Azospirillum brasilense" = "#515151", "Citrobacter koseri" = "darkgrey",
                               "Elizabethkingia meningoseptica" = "#515151", "Enterococcus faecalis" = "darkgrey",
                               "Escherichia coli" = "#ffb400", "Hafnia alvei" = "darkgrey",
                               "Kluyvera intermedia" = "#515151", "Listeria innocua" = "darkgrey",
                               "Microvirga lotononidis" = "#515151", "Paracoccus denitrificans" = "darkgrey",
                               "Pseudomonas chlororaphis" = "#515151", "Pseudomonas putida" =  "darkgrey",
                               "Sphingobacterium multivorum" = "#515151", "Sphingobacterium spiritivorum" = "darkgrey",
                               "Tolumonas auensis" = "#515151", "0" = "white")) +
  scale_y_discrete(limits = rev(unique(epic_some_t25$Species)))+
  scale_shape_manual(values = c(21, 22)) +
  labs(x = "Replicate",
       size = "Relative\nabundance (16S)")+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "#ececec", fill = NA, linewidth = 2),
        axis.text.y=element_text(face = "italic"),
        panel.grid.minor.x = element_line(colour = "#ececec", linewidth = 0.4),
        strip.background = element_rect(fill = "#f1ec3e"),
        text = element_text(size = 15))+
  guides(colour = "none", fill = "none", shape = "none",
         size = guide_legend(override.aes = list(shape = 22, fill = "black")))

BUBBLEPLOT_25_grey



#*Plot for 37 degree treatment results----
BUBBLEPLOT_37_grey <- ggplot(epic_some_t37, aes(x = replicate, y = Species,
                                                 shape = ifelse(Species_presence == 1, "22", "21"),
                                                 colour = Species,
                                                 stroke = ifelse(Species_presence == 0, 1.5, 0),
                                                 fill = ifelse(Species_presence == 1, Species, 0),
                                                 size = ifelse(Rel_abund > 0, Rel_abund, 0.15))) +
  facet_grid(. ~ concentration)+
  scale_size_continuous(breaks = c(0, 0.01, 0.1, 0.25, 0.5))+
  geom_point() +
  geom_text(label = "\u25B2", size = 3, colour = "red", inherit.aes = F,
            aes(x = replicate, y = ifelse(pucL == 1, y+nuolen_pituus, NA), family="Courier New"))+ #▲
  geom_text(label = "\u25BA", size = 3, colour = "black", inherit.aes = F,
            aes(x = replicate + nuolen_pituus, y = ifelse(tnpL == 1, y, NA), family = "Courier New"))+ #►
  geom_text(label = "\u25BA", size = 1.75, colour = "white", inherit.aes = F,
            aes(x = replicate + nuolen_pituus, y = ifelse(tnpL == 1, y, NA), family = "Courier New"))+ #►
  geom_text(label = "\u25BC", size = 3, colour = "blue", inherit.aes = F,
            aes(x = replicate, y = ifelse(intL == 1, y-nuolen_pituus, NA), family="Courier New"))+ #▼
  geom_text(label = "\u25C4", size = 3, colour = "#00c6ff", inherit.aes = F,
            aes(x = replicate - nuolen_pituus, y = ifelse(pucD == 1, y, NA), family="Courier New"))+ #◄
  scale_fill_manual(values=c("Acidovorax defluvii" = "#515151", "Acinetobacter baumannii" = "darkgrey", 
                             "Aeromonas caviae" = "#515151", "Agrobacterium tumefaciens" = "darkgrey", 
                             "Azospirillum brasilense" = "#515151", "Citrobacter koseri" = "darkgrey",
                             "Elizabethkingia meningoseptica" = "#515151", "Enterococcus faecalis" = "darkgrey",
                             "Escherichia coli" = "#bf0202", "Hafnia alvei" = "darkgrey",
                             "Kluyvera intermedia" = "#515151", "Listeria innocua" = "darkgrey",
                             "Microvirga lotononidis" = "#515151", "Paracoccus denitrificans" = "darkgrey",
                             "Pseudomonas chlororaphis" = "#515151", "Pseudomonas putida" =  "darkgrey",
                             "Sphingobacterium multivorum" = "#515151", "Sphingobacterium spiritivorum" = "darkgrey",
                             "Tolumonas auensis" = "#515151", "0" = "white")) +
  scale_colour_manual(values=c("Acidovorax defluvii" = "#515151", "Acinetobacter baumannii" = "darkgrey", 
                               "Aeromonas caviae" = "#515151", "Agrobacterium tumefaciens" = "darkgrey", 
                               "Azospirillum brasilense" = "#515151", "Citrobacter koseri" = "darkgrey",
                               "Elizabethkingia meningoseptica" = "#515151", "Enterococcus faecalis" = "darkgrey",
                               "Escherichia coli" = "#bf0202", "Hafnia alvei" = "darkgrey",
                               "Kluyvera intermedia" = "#515151", "Listeria innocua" = "darkgrey",
                               "Microvirga lotononidis" = "#515151", "Paracoccus denitrificans" = "darkgrey",
                               "Pseudomonas chlororaphis" = "#515151", "Pseudomonas putida" =  "darkgrey",
                               "Sphingobacterium multivorum" = "#515151", "Sphingobacterium spiritivorum" = "darkgrey",
                               "Tolumonas auensis" = "#515151", "0" = "white")) +
  scale_y_discrete(limits = rev(unique(epic_some_t37$Species)))+
  scale_shape_manual(values = c(21, 22)) +
  labs(x = "Replicate",
       size = "Relative\nabundance (16S)")+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "#ececec", fill = NA, linewidth = 2),
        axis.text.y=element_text(face = "italic"),
        panel.grid.minor.x = element_line(colour = "#ececec", linewidth = 0.4),
        strip.background = element_rect(fill = "#efcaca"),
        text = element_text(size = 15))+
  guides(colour = "none", fill = "none", shape = "none",
         size = guide_legend(override.aes = list(shape = 22, fill = "black")))

BUBBLEPLOT_37_grey

