# Author: Alex Ollhoff
# Description: Plot genetic disimilarity between BRIDG6 diverse parents and Rasmusson
##################################################################################################
library(ggplot2)
library(dplyr)
library(reshape2)

# Pairwise genetic distance between Ras and donor parents
#Distance <- read.csv("~/Documents/PhD/NAM/Population_description/Genetic_distance/NAM_parents_with_genetic_distance.csv", header = T)
Distance <- read.csv("~/Desktop/data_for_thesis/NAM_parents_with_genetic_distance.csv", header = T)
head(Distance)

Distance$Pop_location_sorted = factor(Distance$Pop_location, levels = c('Central European', 'Coastal Mediterranean', 'East African', 'Asian', 'Admixed', 'Un'))

# Sort by proportion disimilarity
Distance_by_diff <- Distance[order(Distance[,11]),]

# Add new sort column
Distance_by_diff$Sort_by_diff <- c(1:88)

# Plot
ggplot() +
  labs(x = "Donor Parents", y = "Proportion of Variants\nDiffering from Rasmusson") +
  theme(strip.background = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(), legend.position = "none", axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +
  geom_bar(data = Distance_by_diff, aes(x = as.numeric(Sort_by_diff), y = as.numeric(Proportion_disimilarity_toRas), fill = factor(Pop_location)), na.rm = T, stat = "identity") +
  scale_fill_manual(name = "NAM Parent\nSubpopulations", values = c("goldenrod1", "violet", "olivedrab3", "turquoise1", "red", "grey55")) +
  #facet_grid( Pop_location_sorted ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_continuous(breaks = c(0, 0.010, 0.020, 0.030, 0.040), labels = c("0", "0.01", "0.02", "0.03", "0.04"), limits = c(0, 0.042), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), labels = Distance_by_diff$NAM_names) +
  coord_flip() 
# Save as portrait 3x7.5 device size in Population_description as "Pairwise_distance_by either genetic/pheno/overalldist"

# Quantify trend
Distance_few <- select(Distance, c(line_name, Proportion_disimilarity_toRas, Pop_location))

gendist <- lm(Proportion_disimilarity_toRas ~ factor(Pop_location), data = Distance_few)
anova(gendist)
summary(gendist)