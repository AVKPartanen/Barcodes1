#*housekeeping----

setwd("C:/wd")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(lattice)
library(vegan)
source(file = "HighstatLibV14.R")

# Data in
data1 <-
  read.csv(
    "/Koe-tulokset_epic.csv",
    sep = ";",
    header = T
  )
str(data1)

data1[, 7:22][data1[, 7:22] > 0] <-
  1 #Change all values above 0 to 1 as epicPCR is not a quantitative method. Only presence absence
data1 #Check

#*Number of species per sample, new column----
data3 <- data1
str(data3)

data3$no_of_species <-
  {no_of_species = apply(
    data1 %>% select('Escherichia.coli':'Microvirga.lotononidis....'),
    1,
    sum
  )}


#*remove samples with unrecognised barcode----
data4 <- subset(data3, origin != "unknown")
data4

data4 <-
  data4 %>% arrange(sample.number.for.sequencing, origin)
tail(data4)


#Samples with zero results are missing, lets add those


#*Parameters used in the function----

k <- 1 #parameter for counting rows

#Different origins
all_origins <-
  rep(c("intD", "intL", "pucD", "pucL", "tnpD", "tnpL"), times = 75)

#each replicate bottle
all_replicates <-
  rep(c(
    rep(1, times = 6),
    rep(2, times = 6),
    rep(3, times = 6),
    rep(4, times = 6),
    rep(5, times = 6)
  ),
  times = 15)

#75 bottles in the experiment, running bottle number
all_smpl_no <- rep(1:75, each = 6)

#Each antibiotic
all_ab <-
  rep(c(
    rep('tc', times = 10 * 6),
    rep('su', times = 10 * 6),
    rep('no', times = 5 * 6)
  ), times = 3)

#All different concentrations
all_concentrations <-
  rep(c(
    rep(0.02, times = 5 * 6),
    rep(0.2, times = 5 * 6),
    rep(2, times = 5 * 6),
    rep(20, times = 5 * 6),
    rep(0, times = 5 * 6)
  ), times = 3)

#Three temperature treatments
all_temperatures <-
  c(rep(15, times = 5 * 6 * 5),
    rep(25, times = 5 * 6 * 5),
    rep(37, times = 5 * 6 * 5))

#Create new parameter which to add the missing rows to 
data6 <- data4

head(data6)
tail(data4)

#*Function to add the missing rows----
lisaa_puuttuvat_rivit <- function(datataulukko) {
  k <- 1
  while (k <= 450) {
    smpl_no1 <- datataulukko[,"sample.number.for.sequencing"]
    origin1 <- datataulukko[, "origin"]
    if (smpl_no1[k] == all_smpl_no[k] &&
        origin1[k] == all_origins[k]) {
      k <- k + 1
    }
    else {
      datataulukko <- dplyr::add_row(
        datataulukko,
        origin = all_origins[k],
        sample.number.for.sequencing = ceiling(k / 6),
        ab = all_ab[k],
        concentration = all_concentrations[k],
        temperature = all_temperatures[k],
        replicate = all_replicates[k],
        .before = k
      )
      k <- k + 1
    }
  }
  #changes NAs to 0
  for (j in 1:ncol(datataulukko)) {
    datataulukko[[j]][is.na(datataulukko[[j]])] = 0
  }
  return(datataulukko)
}

#Call for function
data6 <- lisaa_puuttuvat_rivit(data6)

head(data6)
str(data6) #Check that right amount of rows and numbers in right order


#Change sample number, antibiotic, concentration etc. to factors
data6$origin <- factor(data6$origin)
data6$ab <- factor(data6$ab)
data6$replicate <- factor(data6$replicate)
data6$sample.number.for.sequencing <- factor(data6$sample.number.for.sequencing)

str(data6) #Check



##########################################################################################
#################################_____STATISTICS_____#####################################
##########################################################################################


#*Create a metadata file, add column and create rownames ----
metadata1 <- subset(data6[,1:6])
metadata1$name <- NA

for (i in 1:nrow(metadata1)) {
  metadata1[i,]$name <- paste(metadata1[i,]$sample.number.for.sequencing,
                              metadata1[i,]$ab,
                              metadata1[i,]$concentration,
                              metadata1[i,]$temperature,
                              metadata1[i,]$replicate,
                              metadata1[i,]$origin,
                              "counts.txt",
                              sep = "_")
}

rownames(metadata1)<- metadata1$name

#Add same rownames for main data frame
rownames(data6)<- metadata1$name

#*Remove E. coli results from living origins, as they don't indicate HGT (E. coli being the original host) ----
data62 <- data6

Ec_elossa <- c("pucL", "intL", "tnpL")

for (i in 1:nrow(data62)) {
  if (data62$origin[i] %in% Ec_elossa) {
    data62$Escherichia.coli[i] <- 0
  }
}

#Recalculate the total number of species
data62$no_of_species <-
  {no_of_species = apply(
    data62 %>% select('Escherichia.coli':'Microvirga.lotononidis....'),
    1,
    sum
  )}




##########################################################################################################
####################################### ____DATA EXPLORATION____ ################################################
##########################################################################################################

#This data exploration has been executed according to
#Zuur et al. (2010): A protocol for data exploration to avoid common statistical problems.
#Methods in Ecology and Evolution

#*Missing values?----
colSums(is.na(data62))
#No

#*Outliers ----
outliers_T <- ggplot(data62, aes(x=temperature, y=sample.number.for.sequencing)) +
  geom_point() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
outliers_species <- ggplot(data62, aes(x=no_of_species, y=sample.number.for.sequencing)) +
  geom_point() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
outliers_c <- ggplot(data62, aes(x=concentration, y=sample.number.for.sequencing)) +
  geom_point() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggarrange(outliers_T, outliers_species, outliers_c + rremove("x.text"),
          ncol = 2, nrow = 2)
#No outliers

#* How many observations per----
#...origin
table(data62$origin)

#...bottle
table(data62$sample.number.for.sequencing)

#...bacterial species (unnecessary?)
table(data62$Escherichia.coli)
#15/450
table(data62$Kluyvera.intermedia)
#7/450
table(data62$Pseudomonas.putida)
#8/450
table(data62$Sphingobacterium.spiritivorum)
#12/450
table(data62$Acidovorax.defluvii)
#10/450
table(data62$Aeromonas.caviae)
#11/450
table(data62$Elizabethkingia.meningoseptica)
#1/450
table(data62$Sphingobacterium.multivorum)
#5/450
table(data62$Azospirillum.brasiliense)
#5/450
table(data62$Pseudomonas.chlororaphis)
#5/450
table(data62$Enterococcus.faecalis)
#7/450
table(data62$Agrobacterium.tumefaciens)
#1/450
table(data62$Acinetobacter.baumannii)
#21/450
table(data62$Paracoccus.denitrificans)
#8/450
table(data62$Tolumonas.auensis)
#1/450
table(data62$Microvirga.lotononidis....)
#1/450

#How many times a certain number of total new species have been observed for each treatment-origin combination
table(data62$no_of_species)
#  0   1   2   3   4   5 
#379  49   7   6   8   1


#collinearity between categorical covariates and others has been taken into account in experiment design
#antibiotic and concentration are not independent from each others, others are fine

#Continuous covariates concentration and temperature have been taken into account in experiment design

#*Relations between covariates and explanatory parameters----
#Correlation between temperature and total number of new species?

xyplot(
  jitter(no_of_species, factor = 0.5) ~ jitter(temperature),
  data = data62,
  xlab = "Temperature",
  ylab = "N:o of species",
  main = "all samples",
  strip = function(bg = 'white', ...)
    strip.default(bg = 'white', ...),
  scales = list(
    alternating = TRUE,
    x = list(relation = "free"),
    y = list(relation = "same")
  ),
  panel = function(x, y) {
    panel.grid(h = -1, v = 2)
    panel.points(x, y, col = 1)
    panel.abline(lm(y ~ x))        #Add regression line
  }
)
#weak positive correlation


#What about different origins?
#Integron from living host
xyplot(
  jitter(no_of_species, factor = 0.5) ~ jitter(temperature),
  data = subset(data62, origin == "intL"),
  xlab = "Temperature",
  ylab = "N:o of species",
  main = "integron in living host",
  strip = function(bg = 'white', ...)
    strip.default(bg = 'white', ...),
  scales = list(
    alternating = TRUE,
    x = list(relation = "free"),
    y = list(relation = "same")
  ),
  panel = function(x, y) {
    panel.grid(h = -1, v = 2)
    panel.points(x, y, col = 1)
    panel.abline(lm(y ~ x))        #Add regression line
  }
)
#No change

#plasmid from living host
xyplot(
  jitter(no_of_species, factor = 0.5) ~ jitter(temperature),
  data = subset(data62, origin == "pucL"),
  xlab = "Temperature",
  ylab = "N:o of species",
  main = "plasmid in living host",
  strip = function(bg = 'white', ...)
    strip.default(bg = 'white', ...),
  scales = list(
    alternating = TRUE,
    x = list(relation = "free"),
    y = list(relation = "same")
  ),
  panel = function(x, y) {
    panel.grid(h = -1, v = 2)
    panel.points(x, y, col = 1)
    panel.abline(lm(y ~ x))        #Add regression line
  }
)
#Increase

#Plasmid from dead host
xyplot(
  jitter(no_of_species, factor = 0.5) ~ jitter(temperature),
  data = subset(data62, origin == "pucD"),
  xlab = "Temperature",
  ylab = "N:o of species",
  main = "extracted plasmid",
  strip = function(bg = 'white', ...)
    strip.default(bg = 'white', ...),
  scales = list(
    alternating = TRUE,
    x = list(relation = "free"),
    y = list(relation = "same")
  ),
  panel = function(x, y) {
    panel.grid(h = -1, v = 2)
    panel.points(x, y, col = 1)
    panel.abline(lm(y ~ x))        #Add regression line
  }
)
#Increase


#Transposon from living host
xyplot(
  jitter(no_of_species, factor = 0.5) ~ jitter(temperature),
  data = subset(data62, origin == "tnpL"),
  xlab = "Temperature",
  ylab = "N:o of species",
  main = "transposon in living host",
  strip = function(bg = 'white', ...)
    strip.default(bg = 'white', ...),
  scales = list(
    alternating = TRUE,
    x = list(relation = "free"),
    y = list(relation = "same")
  ),
  panel = function(x, y) {
    panel.grid(h = -1, v = 2)
    panel.points(x, y, col = 1)
    panel.abline(lm(y ~ x))        #Add regression line
  }
)
#Increase

#Some have positive correlation

#*Origin vs n:o of species----
ggplot(data62, aes(x=origin, y=no_of_species)) + 
  geom_boxplot()


#*antibiotic vs n:o of species----
ggplot(data62, aes(x=ab, y=no_of_species)) + 
  geom_boxplot()

#*concentration vs. N.o of species----
#tetracycline
xyplot(
  jitter(no_of_species, factor = 0.5) ~ jitter(concentration, factor = 0.2),
  data = subset(data62, ab == "tc" | ab == "no"),
  xlab = "Concentration",
  ylab = "N:o of hosts",
  strip = function(bg = 'white', ...)
    strip.default(bg = 'white', ...),
  scales = list(
    alternating = TRUE,
    x = list(relation = "free"),
    y = list(relation = "same")
  ),
  panel = function(x, y) {
    panel.grid(h = -1, v = 2)
    panel.points(x, y, col = 1)
    panel.abline(lm(y ~ x))        #Add regression line
  }
)
#Hardly any signal, weak decrease


#*sulfamethazine----
xyplot(
  jitter(no_of_species, factor = 0.5) ~ jitter(concentration, factor = 0.2),
  data = subset(data62, ab == "su" | ab == "no"),
  xlab = "Concentration",
  ylab = "N:o of hosts",
  strip = function(bg = 'white', ...)
    strip.default(bg = 'white', ...),
  scales = list(
    alternating = TRUE,
    x = list(relation = "free"),
    y = list(relation = "same")
  ),
  panel = function(x, y) {
    panel.grid(h = -1, v = 2)
    panel.points(x, y, col = 1)
    panel.abline(lm(y ~ x))        #Add regression line
  }
)
#Hardly any signal, weak increase

#*replicate vs. n:o of species, boxplot----
ggplot(data = data62, aes(y = no_of_species, x = replicate)) +
  geom_boxplot() + 
  xlab("replicate") + ylab("N:o of species") + 
  theme(text = element_text(size=15))
#no effect


#*sample n:o vs n:o of species----
ggplot(data = data62, aes(y = no_of_species, x = sample.number.for.sequencing)) +
  geom_boxplot() + 
  xlab("Sample number for sequencing") + ylab("N:o of species") + 
  theme(text = element_text(size=15))
#Variation differs between bottles, add to model as random effect?


#*How big percentage of results are zeroes?----
sum(data62$no_of_species == 0) / nrow(data62)
#0.8422222, a lot of zero results!
#zero-inflated model?



#######################################################################################################
################################## PERMANOVA ---- #########################################
###########################################################################################

#*convert data to matrix, metadata removed----
str(data62)
data62[7:22] <- lapply(data62[7:22],as.numeric)
str(data62)
data62_mat<-as.matrix(data62[7:22])
str(data62_mat)

#Metadata for PERMANOVA
metadata2.2 <- subset(data62[,1:6])
str(metadata2.2)

identical(rownames(metadata2.2), rownames(data62_mat))
#Check, should give TRUE

#Remove rows with no observations
data62_mat2 <- data62_mat[rowSums(data62_mat[])>0,]

#*New metadata with matching rows ----
matching_rows <- intersect(rownames(metadata2.2), rownames(data62_mat2))
metadata2.3 <- metadata2.2[matching_rows, ]

identical(rownames(metadata2.3), rownames(data62_mat2))
#Check, should give TRUE

#*Lets investigate the variance ----

data63 <- data62[1:22]
data63[7:22] <- lapply(data63[7:22],as.numeric)
str(data63)
data63 <- data63[rowSums(data63[7:22])>0,]

identical(rownames(data63), rownames(data62_mat2))
#Should be TRUE

data63_long <-
  gather(data63,
         species,
         presence,
         Escherichia.coli:Microvirga.lotononidis....,
         factor_key = TRUE)
head(data63_long)


AB_mean <- aggregate(data63_long$presence, list(data63_long$ab), FUN=mean)
#  Group.1          x
#1      no 0.09583333
#2      su 0.10879630
#3      tc 0.10344828

c_mean <- aggregate(data63_long$presence, list(data63_long$concentration), FUN=mean)
#  Group.1          x
#1    0.00 0.09583333
#2    0.02 0.10546875
#3    0.20 0.10096154
#4    2.00 0.08333333
#5   20.00 0.12916667

T_mean <- aggregate(data63_long$presence, list(data63_long$temperature), FUN=mean)
#  Group.1          x
#1      15 0.09166667
#2      25 0.13541667
#3      37 0.09943182

rep_mean <- aggregate(data63_long$presence, list(data63_long$replicate), FUN=mean)
#  Group.1          x
#1       1 0.06944444
#2       2 0.07083333
#3       3 0.14166667
#4       4 0.10576923
#5       5 0.11513158

orig_mean <- aggregate(data63_long$presence, list(data63_long$origin), FUN=mean)
#  Group.1          x
#1    intL 0.1250000
#2    pucD 0.0625000
#3    pucL 0.1125000
#4    tnpL 0.1136364

var(AB_mean$x)
#4.243785e-05
var(c_mean$x)
#0.0002832681
var(T_mean$x)
#0.0005448782
var(rep_mean$x)
#0.0009452693
var(orig_mean$x)
#0.0007756543

#AB < c < T < origin < replicate


adonis2(data62_mat2~ ab*concentration + temperature + origin,
        strata = metadata2.3$sample.number.for.sequencing,
        data=metadata2.3,
        permutations = 9999,
        method = 'jaccard')

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  strata 
#Permutation: free
#Number of permutations: 9999

#adonis2(formula = data62_mat2 ~ ab * concentration + temperature + origin, data = metadata2.3, permutations = 9999,
#method = "jaccard", strata = metadata2.3$sample.number.for.sequencing)
#                  Df SumOfSqs      R2       F Pr(>F)    
# ab                2   1.2794 0.04314  2.3370 0.0001 ***
# concentration     1   0.7052 0.02378  2.5763 0.0001 ***
# temperature       1   2.9619 0.09986 10.8203 0.0001 ***
# origin            3   7.3952 0.24933  9.0053 0.0001 ***
# ab:concentration  1   0.3468 0.01169  1.2668 0.0037 ** 
# Residual         62  16.9715 0.57220                   
# Total            70  29.6599 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



###################################___MANN-WHITNEY ANALYSIS
#*The effect of different factors to the total number of species----
#Significant comparisons (you can run all, but I've removed those that give N.S. results for prettier picture)
no_of_species_comparisons <- list(c('tnpL', 'intL'), c('tnpL','pucD'), c('tnpL', 'pucL'))

#Violin plot of differences between different origins
violinplot3<-ggplot(subset(data62, ! origin %in% c('intD', 'tnpD')), aes(x=origin, y=no_of_species, colour=origin,
                                                                         fill=origin)) +
  geom_violin(alpha=0.5, colour=NA) +
  geom_jitter(width = 0.25, height = 0) +
  scale_color_manual(values=c('#59b0b5', 'black', '#864835', '#049e78'))+
  scale_fill_manual(values=c('#b8d9db', 'grey', '#ffb400','lightgreen'))+
  scale_x_discrete(labels = c("intL" = "Integron, \n alive", "pucD" = "Plasmid,\n dead",
                              "pucL" = "Plasmid,\n alive", "tnpL" = "Transposon,\n alive"))+
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = 'black', fill = 'white'))+
  theme(panel.border = element_rect(colour = 'black'))+
  ylab('Number of new hosts')+
  xlab('Origin')+
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6))+
  theme(legend.position='none', panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = no_of_species_comparisons,method = 'wilcox.test',label = 'p.signif', hide.ns = T)
violinplot3

#According to Mann-Whitney test (Wilcoxon rank sum test), only transposon originating from the living host differed
#from the other origins



#Significant comparisons (you can run all, but I've removed those that give N.S. results for prettier picture)
no_of_species_comparisons_T <- list(c('15','37'), c('25','37'))

#Violin plot of differences between different temperatures
violinplot4<-ggplot(data62, aes(x=as.factor(temperature), y=no_of_species, colour=as.factor(temperature),
                                fill=as.factor(temperature))) +
  geom_violin(alpha=0.5, colour=NA) +
  geom_jitter(width = 0.25, height = 0) +
  scale_color_manual(values=c('#59b0b5', '#ffb400', '#bf0202'))+
  scale_fill_manual(values=c('#b8d9db', '#f1ec3e', '#ff0000'))+
  scale_x_discrete(labels = c("15" = "15 °C", "25" = "25 °C", "37" = "37 °C"))+
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6))+
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = 'black', fill = 'white'))+
  theme(panel.border = element_rect(colour = 'black'))+
  ylab('Number of new hosts')+
  xlab('Temperature treatment')+
  theme(legend.position='none', panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = no_of_species_comparisons_T,method = 'wilcox.test',label = 'p.signif')
violinplot4

#According to Mann-Whitney test (Wilcoxon rank sum test), 37 degree treatment differed statistically significantly
#from the other two. 15 ja 25 degree treatments had no statistically significant difference.