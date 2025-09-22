<img align="left" width="500" src="vrai_logo_.jpg" style="margin-right:-230px"></br></br></br></br>
<h1> Characterizing pelagic trawl partnerships using coordination metrics and behavior classification in the Black Sea</h1></br>

## Overview

This repository contains supplementary code and documentation for the study:</br></br>
Title: **"Characterizing pelagic trawl partnerships using coordination metrics and behavior classification in the Black Sea"**</br></br>
The study investigates cooperative behavior among pelagic trawlers using AIS data. 
Two complementary R scripts are provided for identifying candidate dyads based on vessel proximity, course, displacement, and speed metrics.


<div style="text-align: center;">
    <img src="flowchart.png?raw=true" alt="Graphical Abstract" width="850" height="400" style="display: block; margin: 0 auto;"/>
</div>




## File Descriptions

- **`code1_trip_based_dyadic_extraction.R`**  
  **Purpose :** Identifies dyads based on segmented trips between ports.
  **Key steps :**
  - Spatial matching of vessel positions to port locations to identify trip segments.
  - Temporal dyadic reconstruction within trips.
  - Behavioral metric calculations (Proximity, DIh, DId, DItheta).
  - Candidate dyad selection and speed-based vessel classification.

- **`code2_timestamp-based_approach.R`**  
  **Purpose :** Dyadic extraction directly based on timestamp-level proximity, independent of trips.
  **Key steps :**
  - Filters by nationality and speed threshold.
  - Exclude port proximity using haversine buffer.
  - Pairwise dyad construction for each timestamp.
  - Calculation of spatial and directional coordination metrics.
  - Classification of vessels based on speed clusters.

## Required R Packages

Both scripts require the following R packages:

```R
install.packages(c("dplyr", "geosphere", "tidyr", "purrr", "data.table", "mixtools", "mclust", "sf"))
```

The script loads:
```R
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(geosphere)
library(sf)
library(mclust)     # BIC selection (higher is better in mclust)
library(mixtools)   # fallback EM, classic BIC (lower is better)
library(stats)
```

## Input Requirements

## Input Data

- **AIS Data :** data frame named `black` with columns:  
  `MMSI`, `DATE.TIME..UTC.`, `LONGITUDE`, `LATITUDE`, `SPEED`, `COURSE`, `CALLSIGN`

Note: `SPEED` is coerced to numeric; Turkish‑flagged vessels are selected via `CALLSIGN` starting with `"TC"`.

- **Ports data :** data frame named `ports` with columns:  
  `lon`, `lat` (WGS84). Used to remove points within a buffer of each port.

All geographic calculations use `geosphere::distHaversine` (meters). The script creates an `sf` layer (`crs = 4326`) only for convenience; all metric calculations remain in geographic coordinates (lon/lat).

## Running Instructions

**1.	Load AIS and port data into your R session.**
**2.	Run either script line-by-line** in an R console or RStudio.
**3.	Final outputs :**
  - A table of candidate dyads with coordination metrics.
  -	Vessel classifications based on speed.
  -	Optional: count of unique MMSIs involved in final dyads.

## Usage Instructions

### Trip-Based Approach (`code1`)

```R
# Load your AIS data
turkey_data <- fread("your_ais_data.csv")

# Coordinate conversion
turkey_data <- prepare_coordinates(turkey_data)

# Segment trips using port data
turkey_data <- segment_trips_by_port(turkey_data, port_csv = "port_data.csv")

# Compute candidate dyads
candidate_dyads <- compute_dyads(turkey_data)

# Classify vessels
vessel_classification <- classify_vessels(turkey_data)

# Filter only pelagic trawl dyads
final_dyads <- add_vessel_class_to_dyads(candidate_dyads, vessel_classification)
```

### Timestamp-Based Approach (`code2`)

```R
ais_data <- fread("your_ais_data.csv")
result_dyads <- identify_dyads_timestamp_based(ais_data)
```

> *Note: Code 2 is more robust to AIS gaps and was preferred in the study.*


#### Workflow
**1. Load & Clean AIS**
**2. Coerce `SPEED` to numeric and **cap** large values (`<= speed_cap`).**
**3. Keep Turkish‑flagged vessels (`CALLSIGN` starts with `"TC"`).**
**4. Remove points within `port_buffer_km` kilometers of any port:**
```R   
   filter_ports(data, ports, buffer_km = port_buffer_km)
``` 
**5. Build Dyads (All Timestamps)**
For each timestamp (`DATE.TIME..UTC.`) where `n > 1` vessels are present:
- Enumerate all unordered pairs with `combn(MMSI, 2)`.
- Join attributes for both vessels (`_1`, `_2` suffixes).
- Compute, per dyad, the sequence length `T`, per‑timestamp `Distance`, and aggregate metrics: `Proximity`, `DIh`, `DId`, `DItheta`.

The resulting table is de‑duplicated into `dyadic_results1` with columns:
`MMSI1, MMSI2, T, Distance, Proximity, DIh, DId, DItheta, speed_diff`.

**6. Candidate Dyads by Thresholds**
Filter dyads as:
```R 
cand_all_ts <- dyadic_results1[
  Proximity > th_prox & DIh > th_dih & DId > th_did & DItheta > th_dith
]
```
**7. Vessel Classification via Speed GMMs**
- Build the **speed vector** for each MMSI present in candidate dyads.
- **BEST** model selection: `mclustBIC` across `G ∈ {2,3,4}`; pick the model with **highest** BIC.  
  If `mclust` fails, fallback to `mixtools::normalmixEM` and select by **classic** BIC (lower is better).
- Label as pelagic_trawl if component means include both:
  - Towing mode within 2–4 kn, and with minimum sample shares `p24 ≥ 0.05` and `p610 ≥ 0.05`. Otherwise label “other”.
- Two classification outputs:
  - `cls_best_all` (BIC‑selected best G in {2,3,4})
  - `cls_k4_all` (forced `G = 4`)

**8. Merge Classifications → Final Dyads**
Join labels back to dyads and retain only pairs where both MMSIs are `pelagic_trawl`:
- `final_candidate_gmm_best`
- `final_candidate_gmm_k4`

Console summaries report counts of candidate dyads and unique MMSIs per strategy.

## Output

- Table of candidate dyads with coordination metrics
- Vessel classifications
- Optionally: count of unique MMSIs involved in final dyads


## How to Cite

If you use this repository or the scripts code in your research, please cite our work as follows:

Taner Yıldız, Nurdan Cömert, Adriano Mancini, Alessandro Galdelli, Anna Nora Tasetti</br> 
"Characterizing pelagic trawl partnerships using coordination metrics and behavior classification in the Black Sea",
2025




