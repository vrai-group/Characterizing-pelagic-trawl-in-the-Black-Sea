# Load required libraries
library(dplyr)
library(geosphere)
library(tidyr)
library(purrr)
library(data.table)
library(mixtools)
library(stats)
library(sf)

library(dplyr)

# Filter Turkish vessels and speed threshold
turkey_data <- black %>%
  filter(grepl("^TC", CALLSIGN)) %>% 
  filter(SPEED <= 20)

# Function to filter out AIS points near ports
filter_ports <- function(data, ports, buffer_km = 1) {
  data %>%
    rowwise() %>%
    mutate(
      near_port = any(
        distHaversine(cbind(LONGITUDE, LATITUDE), cbind(ports$lon, ports$lat)) / 1000 < buffer_km
      )
    ) %>%
    filter(!near_port) %>%
    ungroup()
}

# Exclude port movements
turkey_data <- filter_ports(turkey_data, ports)

# 1) Assume LATITUDE and LONGITUDE columns exist in turkey_data
#    Convert to sf object (keep original columns with remove = FALSE)
turkey_data <- st_as_sf(
  turkey_data,
  coords = c("LONGITUDE", "LATITUDE"),
  crs = 4326,
  remove = FALSE
)

# 2) For many-to-many joins, "y" side should be a regular data.frame
#    Drop geometry to create non-spatial version
df_turkey_data <- turkey_data %>% 
  st_drop_geometry()

# 3) Construct dyads
dyadic1 <- turkey_data %>%
  group_by(DATE.TIME..UTC.) %>%
  filter(n() > 1) %>%
  summarise(
    dyad_combinations = list(combn(MMSI, 2, simplify = FALSE))
  ) %>%
  unnest(dyad_combinations) %>%
  transmute(
    DATE.TIME..UTC.,
    MMSI1 = map_chr(dyad_combinations, ~ as.character(.x[1])),
    MMSI2 = map_chr(dyad_combinations, ~ as.character(.x[2]))
  ) %>%
  filter(MMSI1 != MMSI2) %>% 
  # ---- 1st left_join (y = df_turkey_data, non-sf) ----
left_join(
  df_turkey_data %>%
    mutate(MMSI = as.character(MMSI)) %>% 
    rename(
      MMSI1 = MMSI,
      LATITUDE1 = LATITUDE,
      LONGITUDE1 = LONGITUDE,
      COURSE1 = COURSE,
      SPEED1 = SPEED
    ),
  by = c("DATE.TIME..UTC.", "MMSI1"),
  relationship = "many-to-many"
) %>%
  # ---- 2nd left_join (y = df_turkey_data, non-sf) ----
left_join(
  df_turkey_data %>%
    mutate(MMSI = as.character(MMSI)) %>%
    rename(
      MMSI2 = MMSI,
      LATITUDE2 = LATITUDE,
      LONGITUDE2 = LONGITUDE,
      COURSE2 = COURSE,
      SPEED2 = SPEED
    ),
  by = c("DATE.TIME..UTC.", "MMSI2"),
  relationship = "many-to-many"
) %>%
  group_by(MMSI1, MMSI2) %>%
  mutate(
    T = n(),  
    Distance = distHaversine(
      cbind(LONGITUDE1, LATITUDE1),
      cbind(LONGITUDE2, LATITUDE2)
    ),
    Proximity = sum(Distance < 1000, na.rm = TRUE) / T,
    DIh = sum(cos((COURSE1 - COURSE2) * pi / 180), na.rm = TRUE) / (T - 1),
    displacement1 = distHaversine(
      cbind(LONGITUDE1, LATITUDE1),
      cbind(lag(LONGITUDE1), lag(LATITUDE1))
    ),
    displacement2 = distHaversine(
      cbind(LONGITUDE2, LATITUDE2),
      cbind(lag(LONGITUDE2), lag(LATITUDE2))
    ),
    speed_diff = abs(SPEED1 - SPEED2) / (SPEED1 + SPEED2 + 1e-6),
    DId = sum(
      1 - (abs(displacement1 - displacement2) / (displacement1 + displacement2 + 1e-6)) * (1 - speed_diff),
      na.rm = TRUE
    ) / (T - 1),
    DItheta = sum(cos((COURSE1 - COURSE2) * pi / 180), na.rm = TRUE) / (T - 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(Distance) & !is.na(DId) & !is.na(DIh))

# 4) Convert to data.table and retain selected columns
library(data.table)

dt_dyadic <- as.data.table(dyadic1)
dyadic_results1 <- unique(
  dt_dyadic[, .(MMSI1, MMSI2, T, Distance, Proximity, DIh, DId, DItheta, speed_diff)]
)

# 5) Filter candidate dyads
candidate_dyadic1 <- dyadic_results1[
  Proximity > 0.5 &
    DIh > 0.5 &
    DId > 0.5 &
    DItheta > 0.5
]

cat("Candidate dyads:
")
print(candidate_dyadic1)

# 6) Vessel classification based on speed profiles
unique_vessels <- unique(c(candidate_dyadic1$MMSI1, candidate_dyadic1$MMSI2))
vessel_class <- list()

df_turkey_data_for_speed <- turkey_data %>% st_drop_geometry()

for (vid in unique_vessels) {
  vessel_speed <- df_turkey_data_for_speed %>%
    filter(MMSI == vid) %>%
    pull(SPEED)

  vessel_speed <- as.numeric(as.character(vessel_speed))
  vessel_speed <- vessel_speed[!is.na(vessel_speed)]

  if (length(vessel_speed) < 10) {
    vessel_class[[as.character(vid)]] <- "insufficient"
    next
  }

  set.seed(123)
  kmeans_result <- kmeans(vessel_speed, centers = 2)
  cluster_centers <- sort(kmeans_result$centers)

  fishing_speed_range <- c(2, 4)
  steaming_speed_range <- c(6, 10)

  if (cluster_centers[1] >= fishing_speed_range[1] &&
      cluster_centers[1] <= fishing_speed_range[2] &&
      cluster_centers[2] >= steaming_speed_range[1] &&
      cluster_centers[2] <= steaming_speed_range[2]) {
    vessel_class[[as.character(vid)]] <- "pelagic_trawl"
  } else {
    vessel_class[[as.character(vid)]] <- "other"
  }
}

vessel_class_dt <- data.table(
  MMSI = as.numeric(names(vessel_class)),
  vessel_type = unlist(vessel_class)
)

candidate_dyadic1[, MMSI1 := as.numeric(MMSI1)]
candidate_dyadic1[, MMSI2 := as.numeric(MMSI2)]

candidate_dyadic1 <- merge(
  candidate_dyadic1,
  vessel_class_dt,
  by.x = "MMSI1",
  by.y = "MMSI",
  all.x = TRUE
)
setnames(candidate_dyadic1, "vessel_type", "MMSI1_type")

candidate_dyadic1 <- merge(
  candidate_dyadic1,
  vessel_class_dt,
  by.x = "MMSI2",
  by.y = "MMSI",
  all.x = TRUE
)
setnames(candidate_dyadic1, "vessel_type", "MMSI2_type")

final_candidate_dyadic1 <- candidate_dyadic1[
  MMSI1_type == "pelagic_trawl" & MMSI2_type == "pelagic_trawl"
]

cat("Final candidate dyads (both vessels are pelagic_trawl):
")
print(final_candidate_dyadic1)

# Extract unique MMSI identifiers from final dyads
unique_mmsi <- unique(c(final_candidate_dyadic1$MMSI1, 
                        final_candidate_dyadic1$MMSI2))

# Report number of unique vessels
length(unique_mmsi)
