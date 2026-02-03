library(sf)
library(dplyr)
library(tidyr)

bc_das <- st_read("lda_000b21a_e.shp") %>%
  filter(PRUID == "59") 

bc_das_data <- bc_das %>%
  st_drop_geometry()

census_original <-read.csv("age_dist_only.csv")

bc_das_pop_data <- census_original %>%
  filter(grepl("years|years", CHARACTERISTIC_NAME, ignore.case = TRUE)) %>%
  filter(grepl("to [0-9]|and over", CHARACTERISTIC_NAME)) %>%
  mutate(char_clean = trimws(CHARACTERISTIC_NAME))

# Find which groups have both female and male counts
age_groups_with_both <- bc_das_pop_data %>%
  group_by(char_clean) %>%
  summarise(
    has_female = any(!is.na(C3_COUNT_WOMEN.) & C3_COUNT_WOMEN.> 0),
    has_male = any(!is.na(C2_COUNT_MEN.) & C2_COUNT_MEN. > 0),
    .groups = "drop"
  ) %>%
  filter(has_male & has_female) %>%
  pull(char_clean)

bc_das_pop_filtered <- bc_das_pop_data %>%
  filter(char_clean %in% age_groups_with_both)

# Pivot female
census_female <- bc_das_pop_filtered %>%
  select(DGUID, char_clean, C3_COUNT_WOMEN.) %>%
  group_by(DGUID, char_clean) %>%
  summarise(count = sum(C3_COUNT_WOMEN., na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = char_clean,
    values_from = count,
    names_prefix = "female",
    values_fill = 0
  )

# Pivot male
census_male <- bc_das_pop_filtered %>%
  select(DGUID, char_clean, C2_COUNT_MEN.) %>%
  group_by(DGUID, char_clean) %>%
  summarise(count = sum(C2_COUNT_MEN., na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = char_clean,
    values_from = count,
    names_prefix = "male",
    values_fill = 0
  )

census_wide <- census_female %>%
  left_join(census_male, by = "DGUID")

male_cols <- names(census_wide)[grepl("^male[0-9]", names(census_wide))] 

for(col in male_cols) {
  age_group <- gsub("^male", "", col)  
  female_col <- paste0("female", age_group)
  total_col <- paste0("total", age_group)
  
  if(female_col %in% names(census_wide)) {
    census_wide[[total_col]] <- census_wide[[col]] + census_wide[[female_col]]
  }
}

# Join census info with shapefile 
bc_das_spatial <- bc_das %>%
  select(DGUID, DAUID, geometry) %>%
  distinct(DGUID, .keep_all = TRUE)

census_wide_spatial <- bc_das_spatial %>%
  left_join(census_wide, by = "DGUID")

#  create age group 
census_wide_spatial <- census_wide_spatial %>%
  mutate(
    # TOTAL population age groups
    pop_total_0_4 = `total0 to 4 years`,
    pop_total_5_9 = `total5 to 9 years`,
    pop_total_10_14 = `total10 to 14 years`,
    pop_total_15_19 = `total15 to 19 years`,
    pop_total_20_24 = `total20 to 24 years`,
    pop_total_25_29 = `total25 to 29 years`,
    pop_total_30_34 = `total30 to 34 years`,
    pop_total_35_39 = `total35 to 39 years`,
    pop_total_40_44 = `total40 to 44 years`,
    pop_total_45_49 = `total45 to 49 years`,
    pop_total_50_54 = `total50 to 54 years`,
    pop_total_55_59 = `total55 to 59 years`,
    pop_total_60_64 = `total60 to 64 years`,
    pop_total_65_69 = `total65 to 69 years`,
    pop_total_70_74 = `total70 to 74 years`,
    pop_total_75_79 = `total75 to 79 years`,
    pop_total_80_84 = `total80 to 84 years`,
    pop_total_85_89 = `total85 to 89 years`,
    pop_total_90_94 = `total90 to 94 years`,
    pop_total_95_99 = `total95 to 99 years`,
    pop_total_100plus = `total100 years and over`,
    
    pop_total_65_74 = pop_total_65_69 + pop_total_70_74,
    pop_total_75_84 = pop_total_75_79 + pop_total_80_84,
    pop_total_85plus = pop_total_85_89 + pop_total_90_94 + pop_total_95_99 + pop_total_100plus,
    pop_total_65plus = pop_total_65_74 + pop_total_75_84 + pop_total_85plus,
    pop_total_75plus = pop_total_75_84 + pop_total_85plus,
    
    # MALE population age groups
    pop_male_65_74 = `male65 to 69 years` + `male70 to 74 years`,
    pop_male_75_84 = `male75 to 79 years` + `male80 to 84 years`,
    pop_male_85plus = `male85 to 89 years` + `male90 to 94 years` + 
      `male95 to 99 years` + `male100 years and over`,
    
    pop_male_65plus = pop_male_65_74 + pop_male_75_84 + pop_male_85plus,
    pop_male_75plus = pop_male_75_84 + pop_male_85plus,
    
    # FEMALE population age groups
    pop_female_65_74 = `female65 to 69 years` + `female70 to 74 years`,
    pop_female_75_84 = `female75 to 79 years` + `female80 to 84 years`,
    
    pop_female_85plus = `female85 to 89 years` +  `female90 to 94 years` + 
      `female95 to 99 years` + `female100 years and over`,
    
    pop_female_65plus = pop_female_65_74 + pop_female_75_84 + pop_female_85plus,
    pop_female_75plus = pop_female_75_84 + pop_female_85plus,
  )

# Verify calculations
census_wide_spatial <- census_wide_spatial %>%
  mutate(
    pop_total_2021 = `total0 to 4 years` + `total5 to 9 years` + 
      `total10 to 14 years` + `total15 to 19 years` +
      `total20 to 24 years` + `total25 to 29 years` +
      `total30 to 34 years` + `total35 to 39 years` +
      `total40 to 44 years` + `total45 to 49 years` +
      `total50 to 54 years` + `total55 to 59 years` +
      `total60 to 64 years` + pop_total_65plus
  )

census_wide_spatial <- census_wide_spatial %>%
  mutate(
    check_total = pop_male_65plus + pop_female_65plus,
    difference = pop_total_65plus - check_total
  )

print("=== VERIFICATION ===")
print("Total should equal male + female (difference should be ~0):")
summary(census_wide_spatial$difference)

print("\nAge group summary:")
summary(census_wide_spatial[c("pop_total_2021", "pop_total_65plus", 
                              "pop_male_65plus", "pop_female_65plus")])

# Save
st_write(census_wide_spatial, "bc_das_population_for_raster.gpkg", delete_layer = TRUE)

print("\n✓ Saved to bc_das_population_for_raster.gpkg")

write.csv(census_wide, "rasterized_census.csv", row.names = FALSE)


st_write(census_wide, "bc_das_population_for_raster.gpkg", delete_layer = TRUE)
saveRDS(census_wide, "bc_das_population_for_raster.rds")