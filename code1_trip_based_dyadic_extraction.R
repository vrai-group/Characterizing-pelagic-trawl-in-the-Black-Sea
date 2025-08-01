# Load required libraries
library(dplyr)
library(data.table)
library(geosphere)
library(sf)
library(tidyr)
library(purrr)
library(mixtools)  # for normalmixEM
library(stats)     # for kmeans

# ---------------------------
# 1. Function: Port Segmentation (Trip Delimitation)
# ---------------------------
segment_trips_by_port <- function(vessel_df, port_csv, port_threshold = 500) {
  # vessel_df: Time series data for vessels; must contain at least "MMSI", "DATE.TIME..UTC.", "LONGITUDE", "LATITUDE" columns.
  # port_csv: Path to the port dataset CSV file; must contain at least "PORT_NAME", "LONGITUDE", "LATITUDE" columns.
  
  # Read port data and convert to sf object
  port_df <- fread(port_csv)
  if (!all(c("PORT_NAME", "LONGITUDE", "LATITUDE") %in% names(port_df))) {
    stop("Port CSV must contain 'PORT_NAME', 'LONGITUDE', and 'LATITUDE' columns.")
  }
  port_sf <- st_as_sf(port_df, coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE)
  port_sf <- st_transform(port_sf, crs = 32636)
  
  # Convert vessel data to sf object (original columns retained)
  vessel_sf <- st_as_sf(vessel_df, coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE)
  vessel_sf <- st_transform(vessel_sf, crs = 32636)
  
  # Calculate minimum distance to any port for each vessel record
  vessel_sf$min_port_distance <- apply(st_distance(vessel_sf, port_sf), 1, min)
  
  # Label as "in port" if below the defined threshold (e.g., 500 m)
  vessel_sf$in_port <- vessel_sf$min_port_distance < port_threshold
  
  # Segment trips by sorting per time and vessel ID
  vessel_sf <- vessel_sf %>%
    group_by(MMSI) %>%
    arrange(DATE.TIME..UTC.) %>%
    mutate(
      prev_in_port = lag(in_port, default = first(in_port)),
      departure = (!in_port & prev_in_port),  # Departure from port
      arrival   = (in_port & !prev_in_port)   # Arrival at port
    ) %>%
    # Create trip_id by incrementing counter at each departure;
    # Assign NA to in-port records
    mutate(
      trip_counter = cumsum(departure),
      trip_id = ifelse(in_port, NA, paste0(MMSI, "_", trip_counter))
    ) %>%
    ungroup()
  
  # Drop unnecessary columns and convert to data.table
  vessel_df_updated <- as.data.table(vessel_sf)
  vessel_df_updated[, c("min_port_distance", "prev_in_port", "departure", "arrival", "trip_counter") := NULL]
  
  return(vessel_df_updated)
}

# ---------------------------
# 2. Function: Coordinate Transformation
# ---------------------------
prepare_coordinates <- function(df, long_col = "LONGITUDE", lat_col = "LATITUDE", utm_epsg = 32636) {
  if (!all(c(long_col, lat_col) %in% names(df))) {
    stop("Required coordinate columns not found.")
  }
  sf_obj <- st_as_sf(df, coords = c(long_col, lat_col), crs = 4326, remove = FALSE)
  sf_obj <- st_transform(sf_obj, crs = utm_epsg)
  coords <- st_coordinates(sf_obj)
  df$x <- coords[,1]
  df$y <- coords[,2]
  return(df)
}

# ---------------------------
# 3. Function: Dyad Calculations (Trip-Based)
# ---------------------------
compute_dyads <- function(df, time_col = "DATE.TIME..UTC.", id_col = "MMSI") {
  # Ensure trip_id column exists
  necessary_cols <- c(time_col, id_col, "LATITUDE", "LONGITUDE", "COURSE", "SPEED", "trip_id")
  if (!all(necessary_cols %in% names(df))) {
    stop("Missing required columns: ", paste(setdiff(necessary_cols, names(df)), collapse = ", "))
  }
  
  # Use only non-port (in-trip) observations: trip_id must not be NA
  df_trip <- df %>% filter(!is.na(trip_id))
  
  # Generate dyad combinations for each time point within each trip
  dyads <- df_trip %>%
    group_by(trip_id, !!sym(time_col)) %>%
    filter(n() > 1) %>%
    summarise(dyad_combinations = list(combn(!!sym(id_col), 2, simplify = FALSE))) %>%
    unnest(dyad_combinations) %>%
    transmute(
      trip_id,
      !!time_col,
      MMSI1 = map_chr(dyad_combinations, ~ as.character(.x[1])),
      MMSI2 = map_chr(dyad_combinations, ~ as.character(.x[2]))
    ) %>%
    filter(MMSI1 != MMSI2) %>%
    # Join with vessel 1 data
    left_join(
      df_trip %>% 
        mutate(MMSI = as.character(!!sym(id_col))) %>% 
        rename(LATITUDE1 = LATITUDE, LONGITUDE1 = LONGITUDE, COURSE1 = COURSE, SPEED1 = SPEED) %>%
        rename(MMSI1 = !!sym(id_col)),
      by = c("trip_id", time_col, "MMSI1")
    ) %>%
    # Join with vessel 2 data
    left_join(
      df_trip %>% 
        mutate(MMSI = as.character(!!sym(id_col))) %>% 
        rename(LATITUDE2 = LATITUDE, LONGITUDE2 = LONGITUDE, COURSE2 = COURSE, SPEED2 = SPEED) %>%
        rename(MMSI2 = !!sym(id_col)),
      by = c("trip_id", time_col, "MMSI2")
    ) %>%
    group_by(trip_id, MMSI1, MMSI2) %>%
    mutate(
      T = n(),
      Distance = distHaversine(cbind(LONGITUDE1, LATITUDE1), cbind(LONGITUDE2, LATITUDE2)),
      Proximity = sum(Distance < 1000, na.rm = TRUE) / T,
      DIh = ifelse(T > 1, sum(cos((COURSE1 - COURSE2) * pi / 180), na.rm = TRUE) / (T - 1), NA_real_),
      displacement1 = distHaversine(cbind(LONGITUDE1, LATITUDE1), cbind(lag(LONGITUDE1), lag(LATITUDE1))),
      displacement2 = distHaversine(cbind(LONGITUDE2, LATITUDE2), cbind(lag(LONGITUDE2), lag(LATITUDE2))),
      speed_diff = abs(SPEED1 - SPEED2) / (SPEED1 + SPEED2 + 1e-6),
      DId = ifelse(T > 1, sum(1 - (abs(displacement1 - displacement2) / (displacement1 + displacement2 + 1e-6)) * (1 - speed_diff), na.rm = TRUE) / (T - 1), NA_real_),
      DItheta = ifelse(T > 1, sum(cos((COURSE1 - COURSE2) * pi / 180), na.rm = TRUE) / (T - 1), NA_real_)
    ) %>%
    ungroup() %>%
    filter(!is.na(Distance) & !is.na(DId) & !is.na(DIh)) %>%
    distinct(trip_id, MMSI1, MMSI2, T, Distance, Proximity, DIh, DId, DItheta, speed_diff)
  
  dyads_dt <- as.data.table(dyads)
  
  candidate_dyads <- dyads_dt[
    Proximity > 0.5 &
      DIh > 0.5 &
      DId > 0.5 &
      DItheta > 0.5
  ]
  return(candidate_dyads)
}

# ---------------------------
# 4. Function: Vessel Classification (Based on Speed Profile)
# ---------------------------
classify_vessels <- function(df, id_col = "MMSI", speed_col = "SPEED", min_obs = 10,
                             fishing_speed_range = c(2,4), steaming_speed_range = c(6,10), seed = 123) {
  vessel_ids <- unique(as.character(df[[id_col]]))
  vessel_class <- list()
  
  for (vid in vessel_ids) {
    vessel_speed <- df %>% filter((!!sym(id_col)) == vid) %>% pull(!!sym(speed_col))
    vessel_speed <- as.numeric(as.character(vessel_speed))
    vessel_speed <- vessel_speed[!is.na(vessel_speed)]
    
    if (length(vessel_speed) < min_obs) {
      vessel_class[[vid]] <- "insufficient"
      next
    }
    
    set.seed(seed)
    k_result <- tryCatch({
      kmeans(vessel_speed, centers = 2)
    }, error = function(e) {
      message("K-means error (vessel ", vid, "): ", e$message)
      return(NULL)
    })
    
    if (is.null(k_result)) {
      vessel_class[[vid]] <- "error"
      next
    }
    
    cluster_centers <- sort(as.numeric(k_result$centers))
    message("Vessel: ", vid, " - Cluster Centers: ", paste(cluster_centers, collapse = ", "))
    
    if (cluster_centers[1] >= fishing_speed_range[1] && cluster_centers[1] <= fishing_speed_range[2] &&
        cluster_centers[2] >= steaming_speed_range[1] && cluster_centers[2] <= steaming_speed_range[2]) {
      vessel_class[[vid]] <- "pelagic_trawl"
    } else {
      vessel_class[[vid]] <- "other"
    }
  }
  
  vessel_class_dt <- data.table(MMSI = as.numeric(names(vessel_class)),
                                vessel_type = unlist(vessel_class))
  return(vessel_class_dt)
}

# ---------------------------
# 5. Function: Add Vessel Classification to Candidate Dyads
# ---------------------------
add_vessel_class_to_dyads <- function(dyads, vessel_class_dt) {
  dyads[, MMSI1 := as.numeric(MMSI1)]
  dyads[, MMSI2 := as.numeric(MMSI2)]
  
  dyads <- merge(dyads, vessel_class_dt, by.x = "MMSI1", by.y = "MMSI", all.x = TRUE)
  setnames(dyads, "vessel_type", "MMSI1_type")
  
  dyads <- merge(dyads, vessel_class_dt, by.x = "MMSI2", by.y = "MMSI", all.x = TRUE)
  setnames(dyads, "vessel_type", "MMSI2_type")
  
  final_dyads <- dyads[MMSI1_type == "pelagic_trawl" & MMSI2_type == "pelagic_trawl"]
  return(final_dyads)
}

# ---------------------------
# 6. MAIN ANALYSIS WORKFLOW
# ---------------------------
# Assume 'turkey_data' is your vessel data and "port_data.csv" is the port dataset.
# Step 1: Coordinate transformation
turkey_data <- prepare_coordinates(turkey_data)

# Step 2: Trip segmentation by port logic (adds trip_id)
turkey_data <- segment_trips_by_port(turkey_data, port_csv = "port_data.csv", port_threshold = 500)

# Step 3: Compute dyads per trip
candidate_dyads <- compute_dyads(turkey_data)
cat("Trip-based candidate dyads:\n")
print(candidate_dyads)

# Step 4: Vessel classification (based on speed profile)
vessel_classification <- classify_vessels(turkey_data)
cat("Vessel classification results:\n")
print(vessel_classification)

# Step 5: Add vessel class info to candidate dyads
final_candidate_dyads <- add_vessel_class_to_dyads(candidate_dyads, vessel_classification)
cat("Final candidate dyads (both vessels are pelagic_trawl):\n")
print(final_candidate_dyads)
