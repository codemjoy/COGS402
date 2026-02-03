library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(stringr)
library(exactextractr)

bc_das <- st_read("lda_000b21a_e.shp") %>%
  filter(PRUID == "59") %>%
  st_transform(bc_das, crs = 3005)  # BC Albers Projection

census <-read.csv("age_dist_only.csv")

# Pivot total counts
census_total <- census %>%
  select(DGUID, CHARACTERISTIC_NAME, C1_COUNT_TOTAL) %>%
  pivot_wider(
    names_from = CHARACTERISTIC_NAME,
    values_from = C1_COUNT_TOTAL, 
  ) %>%
  rename_with(~paste0("total_", make.names(.)), -DGUID)

# Pivot female counts
census_female <- census %>%
  select(DGUID, CHARACTERISTIC_NAME, C3_COUNT_WOMEN.) %>%
  pivot_wider(
    names_from = CHARACTERISTIC_NAME,
    values_from = C3_COUNT_WOMEN.,
    values_fill = 0
  ) %>%
  rename_with(~paste0("female_", make.names(.)), -DGUID)

# Pivot male counts 
census_male <- census %>%
  select(DGUID, CHARACTERISTIC_NAME, C2_COUNT_MEN.) %>%
  pivot_wider(
    names_from = CHARACTERISTIC_NAME,
    values_from = C2_COUNT_MEN.,
    values_fill = 0
  ) %>%
  rename_with(~paste0("male_", make.names(.)), -DGUID)

# Join all three (total, female, male)
census_wide <- census_total %>%
  left_join(census_female, by = "DGUID") %>%
  left_join(census_male, by = "DGUID") 
  
# Join census info with shapefile 
census_wide <- bc_das %>%
  left_join(census_wide, by = "DGUID")

# Remove the "X...." prefix from column names
names(census_wide) <- gsub("^(total|male|female)_(X\\.)+", "\\1_", names(census_wide))
names(census_wide) <- gsub("\\.+", ".", names(census_wide)) 
names(census_wide) <- gsub("\\.$", "", names(census_wide)) 
names(census_wide) <- gsub("_\\.", "_", names(census_wide))


#  create age group 
census_wide <- census_wide %>%
  mutate (
    # TOTAL population age groups
    pop_total_65_74 = `total_65.to.69.years` + `total_70.to.74.years`,
    pop_total_75_84 = `total_75.to.79.years` + `total_80.to.84.years`,
    pop_total_85plus = `total_85.to.89.years` + `total_90.to.94.years` + 
      `total_95.to.99.years` +  `total_100.years.and.over`,
    pop_total_65plus = pop_total_65_74 + pop_total_75_84 + pop_total_85plus,
    pop_total_75plus = pop_total_75_84 + pop_total_85plus,
    pop_total_2021 = `total_0.to.4.years` + `total_5.to.9.years` + 
      `total_10.to.14.years` + `total_15.to.19.years` + 
      `total_20.to.24.years` + `total_25.to.29.years` + 
      `total_30.to.34.years` + `total_35.to.39.years` + 
      `total_40.to.44.years` + `total_45.to.49.years` + 
      `total_50.to.54.years` + `total_55.to.59.years` + 
      `total_60.to.64.years` + pop_total_65plus,
    
    # FEMALE population age groups
    pop_female_65_74 = `female_65.to.69.years` + `female_70.to.74.years`,
    pop_female_75_84 = `female_75.to.79.years` + `female_80.to.84.years`,
    pop_female_85plus = `female_85.to.89.years` + `female_90.to.94.years` + 
      `female_95.to.99.years` + `female_100.years.and.over`,
    pop_female_65plus = pop_female_65_74 + pop_female_75_84 + pop_female_85plus,
    pop_female_75plus = pop_female_75_84 + pop_female_85plus,

    # MALE population age groups
    pop_male_65_74 = `male_65.to.69.years` + `male_70.to.74.years`,
    pop_male_75_84 = `male_75.to.79.years` + `male_80.to.84.years`,
    pop_male_85plus = `male_85.to.89.years` + `male_90.to.94.years` + 
      `male_95.to.99.years` + `male_100.years.and.over`,
    pop_male_65plus = pop_male_65_74 + pop_male_75_84 + pop_male_85plus,
    pop_male_75plus = pop_male_75_84 + pop_male_85plus,
  )


summary(census_wide[c("pop_total_2021", "pop_total_65plus", "pop_total_75plus", 
                      "pop_female_65plus", "pop_male_65plus")])

test_da <- census_wide %>%
  filter(DAUID == "59010129") %>%
  select(DAUID, starts_with("pop_"))

print(test_da)

st_write(census_wide, "bc_das_population_for_raster.shp", delete_layer = TRUE)
saveRDS(census_wide, "bc_das_population_for_raster.rds")