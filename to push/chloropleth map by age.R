library(sf)
library(ggplot2)
library(dplyr)
library(viridis)

Sys.setenv(SHAPE_RESTORE_SHX = "YES")

geo_Data <- st_read("lda_000b21a_e.shp")

age_data <- read.csv("age_dist_only.csv")

# remove non-bc data
library(stringr)
bc_only <- geo_Data %>%
  filter(str_detect(DGUID, "^2021S0512[5][9]"))


age_65plus <- c("65 to 69 years", "70 to 74 years", "75 to 79 years", 
                           "80 to 84 years", "85 to 89 years", "90 to 94 years", 
                           "95 to 99 years", "100 years and over")

age_45plus <- c("45 to 49 years", "50 to 54 years","55 to 59 years", "60 to 64 years",
                "65 to 69 years", "70 to 74 years", "75 to 79 years", 
                "80 to 84 years", "85 to 89 years", "90 to 94 years", 
                "95 to 99 years", "100 years and over")

age_85plus <- c("85 to 89 years", "90 to 94 years", 
                "95 to 99 years", "100 years and over")

age_data$CHARACTERISTIC_NAME <- trimws(age_data$CHARACTERISTIC_NAME)

#Calculate Total Pop, think about how we treat NA
total_pop <- age_data %>%
  group_by(GEO_NAME, DGUID) %>%
  summarise(
    Total_Pop = sum(C2_COUNT_MEN. + C3_COUNT_WOMEN., na.rm = TRUE),
    # Check if ALL values are NA
    All_NA = all(is.na(C2_COUNT_MEN.) & is.na(C3_COUNT_WOMEN.)),
    .groups = "drop"
  ) %>%
  mutate(
    # If all values are NA, set Total_Pop to NA instead of 0
    Total_Pop = ifelse(All_NA, NA, Total_Pop)
  )


# Calculate total population 65+ by location
pop_65plus <- age_data %>%
  filter(CHARACTERISTIC_NAME %in% age_65plus) %>%
  group_by(GEO_NAME, DGUID) %>%
  summarise(
    Pop_65plus = sum(C2_COUNT_MEN. + C3_COUNT_WOMEN., na.rm = TRUE),
    All_NA = all(is.na(C2_COUNT_MEN.) & is.na(C3_COUNT_WOMEN.)),
    .groups = "drop"
  ) %>%
  mutate(
    Pop_65plus = ifelse(All_NA, NA, Pop_65plus)
    )

# Calculate total population 45+ by location
pop_45plus <- age_data %>%
  filter(CHARACTERISTIC_NAME %in% age_45plus) %>%
  group_by(GEO_NAME, DGUID) %>%
  summarise(
    Pop_45plus = sum(C2_COUNT_MEN. + C3_COUNT_WOMEN., na.rm = TRUE),
    All_NA = all(is.na(C2_COUNT_MEN.) & is.na(C3_COUNT_WOMEN.)),
    .groups = "drop"
  ) %>%
  mutate(
    Pop_45plus = ifelse(All_NA, NA, Pop_45plus)
    )

# Calculate total population 85+ by location
pop_85plus <- age_data %>%
  filter(CHARACTERISTIC_NAME %in% age_85plus) %>%
  group_by(GEO_NAME, DGUID) %>%
  summarise(
    Pop_85plus = sum(C2_COUNT_MEN. + C3_COUNT_WOMEN., na.rm = TRUE),
    All_NA = all(is.na(C2_COUNT_MEN.) & is.na(C3_COUNT_WOMEN.)),
    .groups = "drop"
  ) %>%
  mutate(
    Pop_85plus = ifelse(All_NA, NA, Pop_85plus)
  )

#get the unique DGUID for each location
dguid <- age_data %>%
  select (GEO_NAME, DGUID) %>%
  distinct()

demo_data <- total_pop %>%
  left_join(pop_65plus %>% select(DGUID, Pop_65plus), by = "DGUID") %>%
  left_join(pop_85plus %>% select(DGUID, Pop_85plus), by = "DGUID") %>%
  left_join(pop_45plus %>% select(DGUID, Pop_45plus), by = "DGUID") %>%
  mutate(
    # Only calculate percentages if Total_Pop is not NA
    Percent_65plus = ifelse(!is.na(Total_Pop) & Total_Pop > 0, 
                            (Pop_65plus / Total_Pop) * 100, 
                            NA),
    Percent_85plus = ifelse(!is.na(Total_Pop) & Total_Pop > 0, 
                            (Pop_85plus / Total_Pop) * 100, 
                            NA),
    Percent_45plus = ifelse(!is.na(Total_Pop) & Total_Pop > 0, 
                            (Pop_45plus / Total_Pop) * 100, 
                            NA)
    
  ) %>%
  select(-All_NA) 



# View the results
head(demo_data)

# Create age midpoints for each group
age_midpoints <- data.frame(
  CHARACTERISTIC_NAME = c("0 to 4 years", "5 to 9 years", "10 to 14 years", 
                          "15 to 19 years", "20 to 24 years", "25 to 29 years",
                          "30 to 34 years", "35 to 39 years", "40 to 44 years",
                          "45 to 49 years", "50 to 54 years", "55 to 59 years",
                          "60 to 64 years", "65 to 69 years", "70 to 74 years",
                          "75 to 79 years", "80 to 84 years", "85 to 89 years",
                          "90 to 94 years", "95 to 99 years", "100 years and over"),
  age_midpoint = c(2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87, 92, 97, 100)
)

# Calculate weighted mean age (approximation of median)
median_age_data <- age_data %>%
  left_join(age_midpoints, by = "CHARACTERISTIC_NAME") %>%
  group_by(GEO_NAME) %>%
  summarise(
    Median_Age = weighted.mean(age_midpoint, 
                               w = C2_COUNT_MEN. + C3_COUNT_WOMEN., 
                               na.rm = TRUE)
  )

# View results
head(median_age_data)

#match names between two files
map_data <- bc_only %>%
  left_join(demo_data, by = "DGUID")

# Check the join
print(paste("Total areas in map:", nrow(map_data)))
print(paste("Areas with valid data:", sum(!is.na(map_data$Percent_65plus))))
print(paste("Areas with NA (suppressed data):", sum(is.na(map_data$Percent_65plus))))

#adjust projection
map_data_projected <- st_transform(map_data, crs = 3005) # BC Albers

#create Chloropleth map
p <- ggplot(map_data_projected) +  
  geom_sf(aes(fill = Percent_65plus), color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "% Population\n65+",
    na.value = "grey87",
    guide = guide_colorbar(
      barwidth = 1,
      barheight = 10,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(
    title = "Percentage of Population Aged 65+ by Dissemination Area",
    subtitle = "British Columbia, 2021",
    caption = "Source: Statistics Canada Census 2021 Gray areas indicate data suppression or unavailable enumeration"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 15),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position    = c(0.03, 0.72),
    legend.justification = c(0, 1),
    legend.key.size = unit(2, "cm"),
    legend.key.height = unit(8, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.caption = element_text(hjust = 0, size = 11, lineheight = 1.2)
  )

print(p)
ggsave("bc_da_pct_65plus.png", p, width = 8.2, height = 8, dpi = 300)