#### IUCN 
# Grouping fish observations by location, depth zone, and habitat and computing:
#how many are primary/secondary/not fishery targeted 
#Targeted or not by aquarium industry 
#The type of primary and secondary threats
#Representing data as stacked barplots.
#Group 10 to 40 m as shallow; group 50 to 60m in Pigadi as deep; group 80 to 90m in Fourni as deep; show separately deep animal forest inside vs. outside 
#total of 3 stacked barplots per location, 6 stacked barplots total (can all be shown in same graph).
library(tidyr)
library(dplyr)
library(stringr)


IUCN_data <- read.csv("IUCN_compilation.csv")
head(IUCN_data)

 ############## Random data analysis 
##how many are primary/secondary/not fishery targeted
# Summarize the number of species by fishery targeting
fishery_summary <- IUCN_data %>%
  group_by(fishery) %>%
  summarize(count = n(), .groups = "drop")

print(fishery_summary)

species_by_fishery_flat <- IUCN_data %>%
  select(species, fishery) %>%
  distinct() %>%
  arrange(fishery, species)

print(species_by_fishery_flat)

## Targeted or not by aquarium industry 
aquarium_summary <- IUCN_data %>%
  group_by(aquarium) %>%
  summarize(count = n(), .groups = "drop")

print(aquarium_summary)

species_by_aquarium_flat <- IUCN_data %>% 
  select(species, aquarium) %>% 
  distinct() %>% 
  arrange(aquarium, species)

print(species_by_aquarium_flat)

####################################################
####################################################
############ Joining of data sets

fish_transects <- read.csv("fish-transects.csv") %>%
  filter(species != "fishing line") %>%
  left_join(read.csv("IUCN_compilation.csv") %>%
              select(species,IUCN),by="species")

# data set for IUCN was NA uh oh
#########
# Clean to match species 
clean_species <- function(df) {
  df %>%
    mutate(species = str_trim(species) %>% str_to_lower())
}

fish_transects <- clean_species(fish_transects)
IUCN_data <- clean_species(IUCN_data)

# Perform the join second time 
Newdataset2 <- fish_transects %>%
  left_join(IUCN_data, by = "species")

# Check for unmatched species
unmatched_species <- setdiff(fish_transects$species, IUCN_data$species)
print(unmatched_species)

length(unmatched_species)
unique_fish_transects_sp <- unique(fish_transects$species)
length(unique_fish_transects_sp)

# The length within the unmatched species is the same as the unique species within the fish 
# transect data indicating the species are not merging 
#
########################TESTING######################

unique_species <- fish_transects %>%
  select(species) %>%
  distinct()

# Test join with unique species only
test_join <- unique_species %>%
  left_join(IUCN_data, by = "species")

print(test_join)
#### This test still led to the IUCN data being NA #############
###############################################################
#
#
###### Try again 
# Inspect unique species in each dataset
print("Unique species in fish_transects:")
print(unique(fish_transects$species))

print("Unique species in IUCN_data:")
print(unique(IUCN_data$species))

# Find species in fish_transects not in IUCN_data
unmatched_species <- setdiff(fish_transects$species, IUCN_data$species)
print("Unmatched species:")
print(unmatched_species)


#### third try
clean_species <- function(df) {
  df %>%
    mutate(species = str_trim(species) %>%        # Remove spaces
             str_to_lower() %>%             # Convert to lowercase
             str_replace_all("[^a-z ]", ""))  # Remove non-alphanumeric
}

fish_transects <- clean_species(fish_transects)
IUCN_data <- clean_species(IUCN_data)

Newdataset3 <- fish_transects %>%
  left_join(IUCN_data, by = "species")

#####
# Create manual corrections

species_corrections <- tibble(
  original = c("gobie sp", "ind"),
  corrected = c("gobius sp", NA)  # Correct or drop species if appropriate
)

# Apply corrections
fish_transects <- fish_transects %>%
  left_join(species_corrections, by = c("species" = "original")) %>%
  mutate(species = ifelse(!is.na(corrected), corrected, species)) %>%
  select(-corrected)

# Reattempt the join
Newdataset4 <- fish_transects %>%
  left_join(IUCN_data, by = "species")
### did not work :(

