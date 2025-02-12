#*housekeeping----

setwd("C:/wd")
library(phyloseq)
library(tidyverse)
library(microViz)
library(knitr)
library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)

#Read phyloseq object
physeq3 <- readRDS("physeq3.rds")

sample_data(physeq3) <- read_delim("metadata_16S.csv", delim=";") %>%
  column_to_rownames("Samples")

#*Remove negative controls from data ----
#Save negative controls
Samples_toRemove <- c("PCR1-neg-ctrl", "PCR2-neg-ctrl", "PCR3-neg-ctrl")

#To see what samples get removed, this will return a ps object that contains the samples you want to remove
subset_samples(physeq3, Samples_real %in% Samples_toRemove)

#To remove those from your phyloseq object, this will return a new ps object with the samples removed
physeq32 <- subset_samples(physeq3, !(Samples_real %in% Samples_toRemove))


#Get data out from phyloseq objest by converting it to a plot
plot_bar(physeq32, "Samples_real","Abundance", "Species")


#change it to ggplot form for easier editing
p <- plot_bar(physeq32, "Samples_real", fill="Species")
pd <- p$data
pd

head(pd)
dim(pd)
str(pd)
pd$Temperature <- as.character(pd$Temperature)

#Remove rows with Abundance = 0
pd2 <- filter(pd, Abundance > 0)
head(pd2)
str(pd2)

#Change species name maltophila2 -> maltophila (maltophila had two versions of 16S rRNA gene)
pd2$Species[pd2$Species == 'maltophila2'] <- 'maltophila'


#Add genus name to have full species name in Species column
for (i in 1:nrow(pd2)) {
  pd2[i,]$Species <- paste(pd2[i,]$Genus, pd2[i,]$Species,  sep = " ")
}

#There is no other Stenotrophomonas species in our community, so we can change "Stenotrophomonas NA" to
#Stenotrophomonas maltophila
pd2 <- pd2 %>%
  mutate(Species = ifelse(Species == 'Stenotrophomonas NA', 'Stenotrophomonas maltophila', Species))



#*Normalise 16S results based on 16S rRNA gene copy count----

#Table with information about how many 16S rRNA gene copies each strain/species has
kuustoistaSmäärät <- read.csv("/lajien_16S_maarat_Rmuoto.csv",
                              sep = ";")

#Function to normalise the 16S results based on how many copies each strain/species has
normalisoi_16S <- function(taulukko1, taulukko2) {
  # Create a new column in dataframe "one" to store the divided values
  taulukko1$normalisoidut <- NA
  
  # Loop through each row in dataframe "one"
  for (i in 1:nrow(taulukko1)) {
    # Get the value from column "Species" in dataframe "one"
    info <- taulukko1$Species[i]
    
    # Find the corresponding row in dataframe "two" where species matches
    matching_row <- taulukko2[taulukko2$Species == info, ]
    
    # If a matching row is found, divide the abundance in "one" with the values in column "x16S-copies" in "two"
    if (nrow(matching_row) > 0) {
      taulukko1$normalisoidut[i] <- taulukko1$Abundance[i] / matching_row$X16S_copies
    }
  }
  
  # Return the updated dataframe "one"
  return(taulukko1)
}

#Run the function to normalise the 16S results
pd3 <- normalisoi_16S(pd2, kuustoistaSmäärät)


#*Convert 16S data to relative abundance----

#Calculate combined abundance for each sample to be able to use in calculating relative abundance of each species
readsums <- aggregate(normalisoidut ~ Samples_real, data=pd3, FUN=sum)

#Create column for relative abundance to data frame
pd3[ , 'relative_abundance'] = NA

pd3.2 <- pd3

#Function for calculating relative abundance
calculate_relative_abundance <- function(df1, df2) {
  for (i in 1:nrow(df1)) {
    matching_row <- df2[df2$Samples_real == df1$Samples_real[i], ]
    if (nrow(matching_row) > 0) {
      df1$relative_abundance[i] <- df1$normalisoidut[i] / matching_row$normalisoidut
    }
  }
  
  return(df1)
}

# Call the function to calculate relative abundance
pd3.2 <- calculate_relative_abundance(pd3.2, readsums)


#*Group the low abundance species under "Others"----

#Define species with relative abundance under 0.02
Species_others <- c("Acidovorax defluvii", "Acinetobacter baumannii", "Aeromonas caviae", "Azospirillum brasilense",
                    "Citrobacter koseri", "Elizabethkingia meningoseptica", "Enterococcus faecalis",
                    "Escherichia coli", "Hafnia alvei", "Kluyvera intermedia", "Listeria innocua",
                    "Paracoccus denitrificans", "Pseudomonas chlororaphis", "Pseudomonas putida",
                    "Sphingobacterium multivorum", "Sphingobacterium spiritivorum", "Others")

#Function to combine low abundance species to "Others"
muuta_lajit_others <- function(datataulukko) {
  for (i in 1:nrow(datataulukko)) {
    if(!(datataulukko$Species[i] %in% Species_others)){
      datataulukko$Species[i] <- "Others"
    }
  }
  return(datataulukko)
}

#Call the function
pd4.2 <- muuta_lajit_others(pd3.2) 


#*16S Plot----

ABefekti_kuvaaja_RelAbund <- ggplot(data = pd4.2, aes(x = Replicate, y = relative_abundance, fill = factor(Species,
                                                                                                           levels = Species_others))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("Acidovorax defluvii" = "#ff0000", "Acinetobacter baumannii" = "#bf0202", 
                             "Aeromonas caviae" = "darkblue",
                             "Citrobacter koseri" = "lightgreen", "Elizabethkingia meningoseptica" = "#049e78",
                             "Enterococcus faecalis" = "#f1ec3e",
                             "Escherichia coli" = "#ff0000", "Kluyvera intermedia" = "#afcf88",
                             "Paracoccus denitrificans" = "#00c6ff", "Pseudomonas chlororaphis" = "#b8d9db",
                             "Pseudomonas putida" =  "#59b0b5", "Sphingobacterium multivorum" = "#ffb400",
                             "Sphingobacterium spiritivorum" = "#864835", "Hafnia alvei" = "pink",
                             "Listeria innocua" = "black", 
                             "Others"= "gray"),
                    name = "Species")+
  ylab("Normalized relative abundance") +
  theme(legend.position = 'bottom',
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = "gray",
                                          linewidth = 0.4),
        panel.grid.minor.y = element_line(color = "lightgray",
                                          linewidth = 0.2),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(face = "italic"))+
  facet_grid(~Temperature + Concentration)
ABefekti_kuvaaja_RelAbund #muista tallentaa tyyliin 2000 ja 1100, muuten huono laatu!

#Lämpötilavärit
g <- ggplot_gtable(ggplot_build(ABefekti_kuvaaja_RelAbund))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c("#c1e1e9", "#c1e1e9", "#c1e1e9", "#c1e1e9", "#c1e1e9", "#fdf9eb", "#fdf9eb", "#fdf9eb", "#fdf9eb", 
           "#fdf9eb", "#efcaca", "#efcaca", "#efcaca", "#efcaca", "#efcaca")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  g$grobs[[i]]$grobs[[2]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)
