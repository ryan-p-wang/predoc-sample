# ============================================================
# Title: Cleaning and Harmonizing S&P Global Data Center Data
# Author: Ryan Wang
# Project: Senior Essay – When AI comes to Town
# Last updated: 2025-10-29
#
# Description:
# This script processes the raw S&P Global Market Intelligence
# (S&P 451) data center dataset to produce a clean facility-level
# file for analysis.
#
# Inputs:
#   - senior-essay-sample/raw-data/sp_data_center.xlsx
#
# Outputs:
#   - senior-essay-sample/clean-data/sp_data_center_clean.csv
# ============================================================


# --- 0. Environment setup ----------------------------------------------------

rm(list = ls())

library(readxl)
library(dplyr)
library(stringr)
library(sf)
library(tigris)
library(here)

options(tigris_use_cache = TRUE)


# --- 1. Load raw S&P data ----------------------------------------------------

dc <- read_xlsx(
  here("senior-essay-sample", "raw-data", "sp_data_center.xlsx")
)

# Restrict to U.S. facilities
dc <- dc[dc$Country == "USA", ]


# --- 2. Construct lifecycle timing variables --------------------------------

dc <- dc %>%
  mutate(
    year_operational    = as.numeric(str_sub(`Operational Quarter`, 1, 4)),
    quarter_operational = as.numeric(str_extract(`Operational Quarter`, "(?<=Q)\\d")),
    quarter_operational = ifelse(
      is.na(quarter_operational) & !is.na(year_operational),
      NA,
      quarter_operational
    )
  )

dc <- dc %>%
  mutate(
    across(
      c(`Year Built`, `Year Retrofitted`, `Decommissioned Year`, `Decommissioned Quarter`),
      ~ as.numeric(.x)
    )
  ) %>%
  rename(
    year_built            = `Year Built`,
    year_retrofitted      = `Year Retrofitted`,
    year_decommissioned   = `Decommissioned Year`,
    quarter_decommissioned = `Decommissioned Quarter`
  )


dc <- dc %>%
  mutate(
    across(
      c(`Operational Status`, `Datacenter Type`, `Single/Multi Tenant`),
      as.factor
    )
  ) %>%
  rename(
    operational_status = `Operational Status`,
    center_type        = `Datacenter Type`,
    tenant_type        = `Single/Multi Tenant`
  )


# --- 3. Assign county & state using geospatial matching ----------------------

dc <- dc %>%
  mutate(
    across(c(Latitude, Longitude),
    ~ as.numeric(.x)
    )
  ) %>%
  rename(
    latitude  = Latitude,
    longitude = Longitude
  )

dc <- dc %>%
  mutate(
    row_id    = row_number(),
    has_coords = !is.na(latitude) & !is.na(longitude)
  )

# Load county shapefile (2023)
counties_sf <- tigris::counties(cb = TRUE, year = 2023, class = "sf") %>%
  st_transform(4326)

# Load state crosswalk for FIPS codes
state_xwalk <- tigris::fips_codes %>%
  distinct(state_code, state, state_name)

# Add spatial point object for facilities with coordinates
dc_pts <- dc %>%
  filter(has_coords) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Assign county polygons based on spatial point object
dc_joined <- st_join(dc_pts, counties_sf, join = st_within, left = TRUE)

# Keep only county, state abbreviations, and state names
dc_with_geo <- dc_joined %>%
  left_join(state_xwalk, by = c("STATEFP" = "state_code")) %>%
  transmute(
    across(everything()),
    county   = NAME,
    state_ab = state,
    state    = state_name
  ) %>%
  st_drop_geometry() %>%
  select(county, state_ab, state, row_id)

# Merge geocoded data back to main
dc <- dc %>%
  left_join(dc_with_geo, by = "row_id")


# --- 4. Manual county/state corrections --------------------------------------

# Manually fix rows without usable geospatial data
# Location correction based on database city/state information
dc[dc$row_id == 1817, ]$county <- "Cascade"
dc[dc$row_id == 1817, ]$state_ab <- "MT"
dc[dc$row_id == 1817, ]$state <- "Montana"

dc[dc$row_id == 2050, ]$county <- "Fairfax"
dc[dc$row_id == 2050, ]$state_ab <- "VA"
dc[dc$row_id == 2050, ]$state <- "Virginia"

dc[dc$row_id == 2073, ]$county <- "Fairfax"
dc[dc$row_id == 2073, ]$state_ab <- "VA"
dc[dc$row_id == 2073, ]$state <- "Virginia"

dc[dc$row_id == 2432, ]$county <- "Caldwell"
dc[dc$row_id == 2432, ]$state_ab <- "TX"
dc[dc$row_id == 2432, ]$state <- "Texas"

dc[dc$row_id == 2966, ]$county <- "Licking"
dc[dc$row_id == 2966, ]$state_ab <- "OH"
dc[dc$row_id == 2966, ]$state <- "Ohio"

dc[dc$row_id == 2967, ]$county <- "Licking"
dc[dc$row_id == 2967, ]$state_ab <- "OH"
dc[dc$row_id == 2967, ]$state <- "Ohio"

dc[dc$row_id == 3740, ]$county <- "Cass"
dc[dc$row_id == 3740, ]$state_ab <- "TX"
dc[dc$row_id == 3740, ]$state <- "Texas"

dc[dc$row_id == 3760, ]$county <- "Chaves"
dc[dc$row_id == 3760, ]$state_ab <- "NM"
dc[dc$row_id == 3760, ]$state <- "New Mexico"

dc[dc$row_id == 3830, ]$county <- "Hays"
dc[dc$row_id == 3830, ]$state_ab <- "TX"
dc[dc$row_id == 3830, ]$state <- "Texas"

dc[dc$row_id == 4277, ]$county <- "Galveston"
dc[dc$row_id == 4277, ]$state_ab <- "TX"
dc[dc$row_id == 4277, ]$state <- "Texas"

# Remove facilities still lacking county information
dc <- dc[!is.na(dc$county), ]

# Note: I found 6 disagreements between state derived from geospatial data
# and state reported by database. Coordinate check using Google Maps confirms
# validity of geospatial-derived state over reported state
dc$state_og <- dc$`State or Province`
dc_disagree <- dc %>%
  filter(!is.na(state), !is.na(state_og), state != state_og)


# --- 5. Select and rename final analysis variables ---------------------------
dc <- dc %>%
  select(
    Company, `Datacenter Name`, `Datacenter ID`, `MI Datacenter ID`,
    `MI Company ID`, `Facility Owner`, `Street Address`,
    `Total UPS Power (kW)`, `Net UPS Power (kW)`, `Net UPS Power Utilized`,
    `Total Utility Power (kW)`, PUE, `Kilowatts/Rack`, `Total Revenue ($M)`,
    year_built, year_operational, quarter_operational,
    year_retrofitted, year_decommissioned, quarter_decommissioned,
    operational_status, center_type, tenant_type,
    county, state_ab, state, `Postal Code`, latitude, longitude
  )

names(dc)[c(1:11, 13:14, 27)] <- c(
  "company", "datacenter_name", "datacenter_id",
  "mi_datacenter_id", "mi_company_id", "facility_owner",
  "street_address", "total_UPS_power", "net_UPS_power",
  "net_UPS_power_utilized", "total_utility_power",
  "kW_per_rack", "total_revenue", "zip_code"
)


# --- 6. Reconcile operational timing and status ------------------------------

# Backfill missing operational year for facilities marked Operational
dc <- dc %>%
  mutate(
    year_operational = if_else(
      is.na(year_operational) & operational_status == "Operational",
      year_built,
      year_operational
    )
  )

# Correct status for two facilities with year_built after 2025
# Status correction based on Google Map satellite imagery
dc[
  !is.na(dc$operational_status) &
    !is.na(dc$year_built) &
    dc$operational_status == "Operational" &
    dc$year_built == 2026,
]$operational_status <- "Under Construction"

dc[
  !is.na(dc$operational_status) &
    !is.na(dc$year_built) &
    dc$operational_status == "Operational" &
    dc$year_built == 2027,
]$operational_status <- "Planned"

dc$operational_status <- as.factor(dc$operational_status)


# --- 7. Export cleaned dataset -----------------------------------------------

write.csv(
  dc,
  here("senior-essay-sample", "clean-data", "sp_data_center_clean.csv"),
  row.names = FALSE
)