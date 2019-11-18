# infantMilkConsumption


## Table of contents
* [Overview](##Overview)
* [Contents](#Contents)
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Demo](#Demo)
* [License](#License)


## Overview
Code used to create statistical tables and plots used in the manuscript "Formula with Added Sugars Significantly Alters Gut Microbiota and Cognitive Development in Hispanic Infants". The aim of the research is to determine differences in the gut microbiota who consume lactose-reduced formula (a formula that replaces lactose with added sugar) compared to traditional formula and breastfeeding.

## Repo Contents
* [R](#R)
* [man](#man)
* [data](#data)

## System Requirements
This program was tested and developed on the macOS Catalina 10.15 system. The following packages in R are necessary for this package:
* R (>= 3.5.0)
* data.table
* ggplot2
* ggsignif
* phyloseq
* dplyr 
* stringr
* readr

RStudio Version 1.2.5019 is optional.

## Installation

To install type the following lines in R:

```
install.packages("devtools")
library("devtools")
```

You will then need to download the R package "infantMilkConsumption" using the following lines:

```
install.github("rbarner/infantMilkConsumption")
library("infantMilkConsumption")
```

## Demo 
(The runtime of each of these functions should be <30 seconds each.)

1. First we want to determine if there are any differences in the general characteristics of the infants who are either fed solely breastmilk (BF), traditional formula (TF) or lactose-reduced formula with added sugar (ASF). To do this run the following command:
```
compute_possibleCovariates_milkType()
```
This command will give us the means of the three groups and an p-value from an ANOVA testing for differences in means (or p-value from chi-square test if categorical).  We see that there are significant differences in age in days of baby ("baby_age"), BMI of mother ("mom_current_BMI") and consumption of fruit/juice of infants at 6 months ("fruit_incj_Inf_1m6m") so we will include these variables in the statistical model.

2. Next, we determine the differences in taxonomic groups at the phylum level in the infants at 6 months of age based on the three milk consumption types. 
```
compute_taxa_milkType_association()
```
This will give us the means, differences between TF and BF groups, differences between ASF and BF groups and p-values associated with those differences. The phylum Actinobacteria has the greatest statistically significant difference so we will plot this graph.
```
plot_taxa_milkType_association()
```

We will also compute differences in alpha- and beta- diversities
```
compute_alphaDiversity_milkType_association()
```

```
compute_betaDiversity_milkType_association()
```

3. We have determined predicted KEGG modules so we will determine differences in the abundance of these predicted modules by milk consumption type.
```
compute_keggModule_milkType_association()
```
This will give us the means, differences between TF and BF groups, differences between ASF and BF groups and p-values associated with those differences. The KEGG module M00221 is "putative simple sugar transport system" and with the increase in glucose in the ASF group, it is possible that this is high in the ASF group.

```
plot_keggModule_milkType_association()
```

4. Lastly, we want to determine if the type of milk consumption at 6 months has any impact on the development and growth of the child at 24 months. To test for development we will use the scaled scores from the Bayley's Scale of Infant Development.
```
compute_bayley_milkType_association()
```
From this we see that the difference in cognitive score ("bsid_cog_ss") in ASF and BF groups is statistically significant so we will create a boxplot to visualize these differences

```
plot_bayley_milkType_association()
```

We also see that the difference in the score of fine motor skills ("bsid_mot_fm_ss") between ASF and BF groups is also statistically significant so we will create a boxplot to visualize this difference as well.
```
plot_bayley_milkType_association(variable = "bsid_mot_fm_ss",variableName = "Fine Motor (FM) Scaled Score")
```

We also test for differences in growth measures among the groups using the following line:
```
compute_somaticGrowth_milkType_association()
```

## License
This package is licensed under GNU General Public License (GPL-3)
