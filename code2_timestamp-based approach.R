# =========================
# Timestamp-based approach
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(data.table)
  library(geosphere)
  library(sf)
  library(mclust)     # GMM + BIC (higher is better)
  library(mixtools)   # fallback EM (classic BIC lower is better)
  library(stats)
})

set.seed(42)

# -------------------------
# Parameters
# -------------------------
speed_cap      <- 20        # hard cap for SPEED (kn)
port_buffer_km <- 1         # distance to ports for exclusion
# Dyad thresholds (can be tuned / swept later)
th_prox  <- 0.50
th_dih   <- 0.50
th_did   <- 0.50
th_dith  <- 0.50

# -------------------------
# Helper: filter out points near ports
# -------------------------
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

# -------------------------
# 0) Load & clean AIS
# -------------------------
black$SPEED <- suppressWarnings(as.numeric(black$SPEED))

turkey_data <- black %>%
  filter(grepl("^TC", CALLSIGN)) %>%           # Turkish-flag callsign (TC- prefix)
  filter(is.finite(SPEED), SPEED <= speed_cap)

# Exclude port movements (assumes `ports` has columns: lon, lat)
turkey_data <- filter_ports(turkey_data, ports, buffer_km = port_buffer_km)

# -------------------------
# 1) Build dyads (all timestamps)
# -------------------------
turkey_sf <- st_as_sf(
  turkey_data,
  coords = c("LONGITUDE", "LATITUDE"),
  crs = 4326,
  remove = FALSE
)

df_turkey_data <- turkey_sf %>% st_drop_geometry()

dyadic1 <- turkey_sf %>%
  group_by(DATE.TIME..UTC.) %>%
  filter(n() > 1) %>%
  summarise(dyad_combinations = list(combn(MMSI, 2, simplify = FALSE)), .groups = "drop") %>%
  unnest(dyad_combinations) %>%
  transmute(
    DATE.TIME..UTC.,
    MMSI1 = map_chr(dyad_combinations, ~ as.character(.x[1])),
    MMSI2 = map_chr(dyad_combinations, ~ as.character(.x[2]))
  ) %>%
  filter(MMSI1 != MMSI2) %>%
  # join vessel 1
  left_join(
    df_turkey_data %>%
      mutate(MMSI = as.character(MMSI)) %>%
      rename(MMSI1 = MMSI,
             LATITUDE1 = LATITUDE, LONGITUDE1 = LONGITUDE,
             COURSE1 = COURSE, SPEED1 = SPEED),
    by = c("DATE.TIME..UTC.", "MMSI1")
  ) %>%
  # join vessel 2
  left_join(
    df_turkey_data %>%
      mutate(MMSI = as.character(MMSI)) %>%
      rename(MMSI2 = MMSI,
             LATITUDE2 = LATITUDE, LONGITUDE2 = LONGITUDE,
             COURSE2 = COURSE, SPEED2 = SPEED),
    by = c("DATE.TIME..UTC.", "MMSI2")
  ) %>%
  group_by(MMSI1, MMSI2) %>%
  mutate(
    T = n(),
    Distance = distHaversine(cbind(LONGITUDE1, LATITUDE1),
                             cbind(LONGITUDE2, LATITUDE2)),
    Proximity = sum(Distance < 1000, na.rm = TRUE) / T,
    DIh = ifelse(T > 1, sum(cos((COURSE1 - COURSE2) * pi / 180), na.rm = TRUE)/(T - 1), NA_real_),
    displacement1 = distHaversine(cbind(LONGITUDE1, LATITUDE1),
                                  cbind(dplyr::lag(LONGITUDE1), dplyr::lag(LATITUDE1))),
    displacement2 = distHaversine(cbind(LONGITUDE2, LATITUDE2),
                                  cbind(dplyr::lag(LONGITUDE2), dplyr::lag(LATITUDE2))),
    speed_diff = abs(SPEED1 - SPEED2) / (SPEED1 + SPEED2 + 1e-6),
    DId = ifelse(T > 1,
                 sum(1 - (abs(displacement1 - displacement2) / (displacement1 + displacement2 + 1e-6)) * (1 - speed_diff),
                     na.rm = TRUE) / (T - 1),
                 NA_real_),
    DItheta = ifelse(T > 1, sum(cos((COURSE1 - COURSE2) * pi / 180), na.rm = TRUE)/(T - 1), NA_real_)
  ) %>%
  ungroup() %>%
  filter(!is.na(Distance) & !is.na(DId) & !is.na(DIh))

dt_dyadic <- as.data.table(dyadic1)
dyadic_results1 <- unique(
  dt_dyadic[, .(MMSI1, MMSI2, T, Distance, Proximity, DIh, DId, DItheta, speed_diff)]
)

# -------------------------
# 2) Candidate dyads by thresholds
# -------------------------
cand_all_ts <- dyadic_results1[
  Proximity > th_prox &
    DIh > th_dih &
    DId > th_did &
    DItheta > th_dith
]

cat("Candidate dyads (all timestamps):", nrow(cand_all_ts), "\n")

# -------------------------
# 3) GMM classification for ALL candidate MMSIs
#    (i) BEST: k in {2,3,4}, BIC selection
#    (ii) K=4 forced
# -------------------------
df_speed <- df_turkey_data %>% as.data.table()
uni_mmsi <- sort(unique(c(dyadic_results1$MMSI1, dyadic_results1$MMSI2)))

bic_from_loglik <- function(loglik, k, n){
  p <- 3*k - 1
  -2*loglik + p*log(n)  # classic BIC: lower is better
}

get_best_mclust_fit <- function(x, Gset = 2:4){
  bic <- tryCatch(mclustBIC(x, G = Gset, modelNames = c("E","V")), error=function(e) NULL)
  if (is.null(bic) || !any(is.finite(bic))) return(NULL)
  
  cand <- expand.grid(model=c("E","V"), G=Gset, KEEP.OUT.ATTRS = FALSE)
  
  pull_bic <- function(m,g){
    if (is.null(dim(bic))) return(NA_real_)
    rn <- rownames(bic); cn <- colnames(bic)
    if (is.null(rn) || is.null(cn)) return(NA_real_)
    if (m %in% rn && as.character(g) %in% cn) as.numeric(bic[m, as.character(g)]) else NA_real_
  }
  
  cand$BIC <- mapply(pull_bic, cand$model, cand$G)
  cand <- cand[is.finite(cand$BIC),]
  if (!nrow(cand)) return(NULL)
  
  best <- cand[which.max(cand$BIC),]
  fit  <- tryCatch(Mclust(x, G = best$G, modelNames = best$model), error=function(e) NULL)
  if (is.null(fit) || is.null(fit$parameters$mean)) return(NULL)
  
  list(k = best$G,
       mu = as.numeric(fit$parameters$mean),
       sd = sqrt(as.numeric(fit$parameters$variance$sigmasq)),
       pi = as.numeric(fit$parameters$pro),
       rule = "mclust")
}

get_best_mixtools_fit <- function(x, Gset = 2:4){
  best <- NULL; best_bic <- Inf
  for (k in Gset){
    f <- tryCatch(normalmixEM(x, k=k, maxit=1000, verb=FALSE), error=function(e) NULL)
    if (!is.null(f)){
      bic <- bic_from_loglik(f$loglik, k, length(x))
      if (is.finite(bic) && bic < best_bic){
        best_bic <- bic
        best <- list(k=k, mu=as.numeric(f$mu), sd=as.numeric(f$sigma), pi=as.numeric(f$lambda), rule="mixtools")
      }
    }
  }
  best
}

classify_one_gmm_best <- function(speeds, min_n=50, tow_band=c(2,4), steam_band=c(6,10)){
  x <- as.numeric(speeds); x <- x[is.finite(x)]
  if (length(x) < min_n || sd(x) < 1e-6 || length(unique(x)) < 3)
    return(list(vessel_type=NA_character_, k=NA_integer_))
  
  fit <- get_best_mclust_fit(x, 2:4)
  if (is.null(fit)) fit <- get_best_mixtools_fit(x, 2:4)
  if (is.null(fit)) return(list(vessel_type=NA_character_, k=NA_integer_))
  
  mu <- fit$mu
  has_tow   <- any(mu >= tow_band[1]   & mu <= tow_band[2])
  has_steam <- any(mu >= steam_band[1] & mu <= steam_band[2])
  p24  <- mean(x >= tow_band[1]   & x <= tow_band[2])
  p610 <- mean(x >= steam_band[1] & x <= steam_band[2])
  
  vessel_type <- if (has_tow && has_steam && p24>=0.05 && p610>=0.05) "pelagic_trawl" else "other"
  list(vessel_type=vessel_type, k=fit$k)
}

classify_one_gmm_k4 <- function(speeds, min_n=50, tow_band=c(2,4), steam_band=c(6,10)){
  x <- as.numeric(speeds); x <- x[is.finite(x)]
  if (length(x) < min_n || sd(x) < 1e-6 || length(unique(x)) < 3)
    return(list(vessel_type=NA_character_))
  
  fit <- tryCatch(Mclust(x, G=4, modelNames=c("E","V")), error=function(e) NULL)
  if (is.null(fit) || is.null(fit$parameters$mean)){
    f <- tryCatch(normalmixEM(x, k=4, maxit=1000, verb=FALSE), error=function(e) NULL)
    if (is.null(f)) return(list(vessel_type=NA_character_))
    mu <- as.numeric(f$mu)
  } else {
    mu <- as.numeric(fit$parameters$mean)
  }
  
  has_tow   <- any(mu >= tow_band[1]   & mu <= tow_band[2])
  has_steam <- any(mu >= steam_band[1] & mu <= steam_band[2])
  p24  <- mean(x >= tow_band[1]   & x <= tow_band[2])
  p610 <- mean(x >= steam_band[1] & x <= steam_band[2])
  
  vessel_type <- if (has_tow && has_steam && p24>=0.05 && p610>=0.05) "pelagic_trawl" else "other"
  list(vessel_type=vessel_type)
}

speeds_by_mmsi <- df_speed[MMSI %in% uni_mmsi, .(SPEED = list(SPEED)), by = MMSI]

cls_best_all <- speeds_by_mmsi[, {
  c <- classify_one_gmm_best(unlist(SPEED))
  .(vessel_type = c$vessel_type, k_chosen = c$k)
}, by = MMSI]

cls_k4_all <- speeds_by_mmsi[, {
  c <- classify_one_gmm_k4(unlist(SPEED))
  .(vessel_type = c$vessel_type)
}, by = MMSI]

cat("BEST coverage:", nrow(cls_best_all), "vessels\n")
cat("K=4  coverage:", nrow(cls_k4_all),  "vessels\n")
cat("BEST labels:\n"); print(table(cls_best_all$vessel_type, useNA="ifany"))
cat("K=4  labels:\n"); print(table(cls_k4_all$vessel_type,  useNA="ifany"))

# -------------------------
# 4) Merge classifications → final candidates
# -------------------------
mk_final <- function(cand_dt, class_dt){
  cand <- copy(as.data.table(cand_dt))
  cand[, `:=`(MMSI1 = as.numeric(MMSI1), MMSI2 = as.numeric(MMSI2))]
  drop_cols <- grep("(^vessel_type$|_type$)", names(cand), value = TRUE)
  if (length(drop_cols)) cand[, (drop_cols) := NULL]
  
  cls <- copy(as.data.table(class_dt))[
    , .(MMSI = as.numeric(MMSI), vessel_type = as.character(vessel_type))
  ]
  cls <- unique(cls, by = "MMSI")
  
  c1 <- merge(cand, cls, by.x="MMSI1", by.y="MMSI", all.x=TRUE)
  data.table::setnames(c1, "vessel_type", "MMSI1_type")
  c2 <- merge(c1,   cls, by.x="MMSI2", by.y="MMSI", all.x=TRUE)
  data.table::setnames(c2, "vessel_type", "MMSI2_type")
  
  c2[MMSI1_type=="pelagic_trawl" & MMSI2_type=="pelagic_trawl"]
}

final_candidate_gmm_best <- mk_final(cand_all_ts, cls_best_all)
final_candidate_gmm_k4   <- mk_final(cand_all_ts, cls_k4_all)

cat("FINAL (BEST):", nrow(final_candidate_gmm_best),
    "| FINAL (K=4):", nrow(final_candidate_gmm_k4), "\n")

mmsi_best <- sort(unique(c(final_candidate_gmm_best$MMSI1, final_candidate_gmm_best$MMSI2)))
mmsi_k4   <- sort(unique(c(final_candidate_gmm_k4$MMSI1,   final_candidate_gmm_k4$MMSI2)))
cat("Unique MMSI (BEST):", length(mmsi_best), " | Unique MMSI (K=4):", length(mmsi_k4), "\n")
