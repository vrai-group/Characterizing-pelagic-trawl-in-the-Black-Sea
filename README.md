<img align="left" width="500" src="vrai_logo_.jpg" style="margin-right:-230px"></br></br></br></br>
<h1> Characterizing pelagic trawl partnerships using coordination metrics and behavior classification in the Black Sea</h1></br>

## Overview

This repository contains supplementary code and documentation for the study:
Title: **"Characterizing pelagic trawl partnerships using coordination metrics and behavior classification in the Black Sea"**
The study investigates cooperative behavior among pelagic trawlers using AIS data. 
Two complementary R scripts are provided for identifying candidate dyads based on vessel proximity, course, displacement, and speed metrics.


<div style="text-align: center;">
    <img src="GA1.jpeg?raw=true" alt="Graphical Abstract" width="850" height="400" style="display: block; margin: 0 auto;"/>
</div>




## File Descriptions

- **`code1_trip_based_dyadic_extraction.R`**  
  Identifies dyads based on segmented trips between ports.
  - Matches AIS points to ports
  - Segments fishing trips
  - Calculates behavior metrics (Proximity, DIh, DId, DItheta)
  - Selects dyad candidates and classifies vessels

- **`code2_timestamp-based_approach.R`**  
  Extracts dyads based directly on proximity at each timestamp.
  - Filters by nationality and speed
  - Removes port-adjacent points
  - Pairs vessels based on spatial/directional alignment
  - Clusters speed profiles for classification

## Required R Packages

Install the required packages:

```R
install.packages(c("dplyr", "geosphere", "tidyr", "purrr", "data.table", "mixtools", "sf"))
```

## Input Data

- **AIS Data:** CSV file with columns:  
  `MMSI`, `DATE.TIME..UTC.`, `LONGITUDE`, `LATITUDE`, `SPEED`, `COURSE`, `CALLSIGN`
- **Port Data (for Code 1 only):** CSV file with columns:  
  `PORT_NAME`, `LONGITUDE`, `LATITUDE`

## Usage Instructions

### Trip-Based Method (`code1`)

```R
turkey_data <- fread("your_ais_data.csv")
turkey_data <- prepare_coordinates(turkey_data)
turkey_data <- segment_trips_by_port(turkey_data, port_csv = "port_data.csv")
candidate_dyads <- compute_dyads(turkey_data)
vessel_classification <- classify_vessels(turkey_data)
final_dyads <- add_vessel_class_to_dyads(candidate_dyads, vessel_classification)
```

### Timestamp-Based Method (`code2`)

```R
ais_data <- fread("your_ais_data.csv")
result_dyads <- identify_dyads_timestamp_based(ais_data)
```

> *Note: Code 2 is more robust to AIS gaps and was preferred in the study.*

## Output

- Table of candidate dyads with coordination metrics
- Vessel classifications
- Optionally: count of unique MMSIs involved in final dyads


## How to Cite

If you use this repository or the scripts code in your research, please cite our work as follows:

Taner Yıldız1, Nurdan Cömert, Adriano Mancini, Alessandro Galdelli, Anna Nora Tasetti 
"Characterizing pelagic trawl partnerships using coordination metrics and behavior classification in the Black Sea",
2025




