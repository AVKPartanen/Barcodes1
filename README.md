# Barcode project
__This project contains code for a scientific article that is under submission process. Abstract of the project and link to article will be added later.__

## About code dependencies
The data frame created in Relative abundance data analysis is further used in bubble plot combining relative abundance data and epicPCR carrier data  
[Relative abundance data](Abundance_data.R) -> 
[EpicPCR and 16S plot](Bubble_plot_with_epic_16S_and_BC_data.R)

## The plot combining epicPCR and 16S data has been finalised outside R
The code for bubble plots with relative abundance data and epicPCR gene carrier data draws most of the plots, but they have been finalised for more polished look outside R. There the “not observed” bubble was added to relative abundance legend and the legend for the origin markers was introduced. I also changed the antibiotic concentration on top of each group to text including antibiotic and written information whether it was the treatment with lower or higher concentration, and I added the temperature on top of the plot.
[EpicPCR and 16S plot](Bubble_plot_with_epic_16S_and_BC_data.R)

## Where to find the DADA pipeline used and changes compared to that
The DADA pipeline to create the phyloseq object that is used in analysis of the 16S results can be found from (https://benjjneb.github.io/dada2/tutorial.html) .  The filtering was done with settings:
-	trimming length was 220 and 200
-	maximum expected error was 2 and 2  
The reads were assigned taxonomy based on 16S data base from the strains and species instead of SILVA database.  Full genomes of some of the HAMBI strains used has been published by  __Hogle et al. (2024):__ Complete genome sequences of 30 bacterial species from a synthetic community. _Microbiology Resource Announcements_, __13__, (6).


# Otsikko 
## Alaotsikko

Tekstia samalla rivilla. 
Jatkuu

Eri rivilla.  
Jartkuu.

_kursiivi_
__lihavoitu__


```r
koodiblokki
```

[Statistical data analysis](Statistical_analysis.R)
