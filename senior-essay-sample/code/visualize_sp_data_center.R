# ============================================================
# Title: Visualizing S&P Global Data Center Patterns
# Author: Ryan Wang
# Project: Senior Essay – AI Data Center Crowding-Out
# Last updated: 2025-10-30
#
# Description:
# This script reads the cleaned S&P Global Market Intelligence
# (S&P 451) data center dataset and produces descriptive figures
# on:
#   1. The operational stock of U.S. data centers over time.
#   2. Total utility power demand of operational centers (GW).
#   3. County-level distribution of operational centers in 2025.
#   4. The distribution of rack power density (kW per rack) over
#      time, with emphasis on AI-relevant thresholds.
#
# Inputs:
#   - senior-essay-sample/clean-data/sp_data_center_clean.csv
#
# Outputs:
#   - senior-essay-sample/output/operational_stock.png
#   - senior-essay-sample/output/utility_demand.png
#   - senior-essay-sample/output/county_map.png
#   - senior-essay-sample/output/rack_power_density.png
# ============================================================


# --- 0. Environment setup ----------------------------------------------------

rm(list = ls())

library(dplyr)
library(ggplot2)
library(scales)
library(stringr)
library(sf)
library(tidyr)
library(patchwork)
library(tigris)
library(here)

options(tigris_use_cache = TRUE)


# --- 1. Load cleaned data ----------------------------------------------------

dc <- read.csv(
  here("senior-essay-sample", "clean-data", "sp_data_center_clean.csv"),
  stringsAsFactors = FALSE
)


# --- 2. Figure 1: Operational stock of U.S. data centers ---------------------

starts <- dc %>%
  filter(!is.na(year_operational)) %>%
  transmute(year = year_operational, delta = 1)

# Assume data centers are decommissioned at end of specified year
exits <- dc %>%
  filter(!is.na(year_decommissioned)) %>%
  transmute(year = year_decommissioned + 1, delta = -1)

net_add <- bind_rows(starts, exits) %>%
  group_by(year) %>%
  summarise(net_additions = sum(delta), .groups = "drop")

year_range <- seq(1980, 2025, by = 1)

timeline <- tibble(year = year_range) %>%
  left_join(net_add, by = "year") %>%
  mutate(
    net_additions      = replace_na(net_additions, 0),
    stock_operational  = cumsum(net_additions)
  )

(f1 <- ggplot(timeline, aes(x = year, y = stock_operational)) +
    geom_area(fill = "#2E86AB", alpha = 0.25) +
    geom_line(linewidth = 1.1, color = "#2E86AB") +
    labs(
      title    = "Operational Stock of U.S. Data Centers",
      subtitle = "Based on S&P 451 Database",
      x        = "Year",
      y        = "Number of Data Centers",
      caption  = "Stock = cumulative starts minus cumulative decommissions"
    ) +
    theme_minimal(base_size = 13) +
    scale_x_continuous(
      breaks = seq(
        min(timeline$year, na.rm = TRUE),
        max(timeline$year, na.rm = TRUE),
        by = 5
      )
    ) +
    theme(
      plot.background    = element_rect(fill = "white", color = NA),
      plot.title.position = "plot",
      plot.title         = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle      = element_text(hjust = 0.5),
      plot.caption       = element_text(hjust = 0.5, size = 10),
      axis.title         = element_text(face = "bold")
    ))

ggsave(
  filename = here("senior-essay-sample", "output", "fig1_operational_stock.png"),
  plot     = f1,
  width    = 8,
  height   = 5,
  dpi      = 300
)


# --- 3. Figure 2: Total utility power demand (GW) ----------------------------

starts <- dc %>%
  filter(!is.na(year_operational), !is.na(total_utility_power)) %>%
  transmute(
    year     = year_operational,
    delta_kw = +total_utility_power
  )

exits <- dc %>%
  filter(!is.na(year_decommissioned), !is.na(total_utility_power)) %>%
  transmute(
    year     = year_decommissioned + 1,
    delta_kw = -total_utility_power
  )

net_add <- bind_rows(starts, exits) %>%
  group_by(year) %>%
  summarize(net_additions_kw = sum(delta_kw, na.rm = TRUE), .groups = "drop") %>%
  arrange(year)

timeline <- tibble(year = year_range) %>%
  left_join(net_add, by = "year") %>%
  mutate(
    net_additions_kw      = replace_na(net_additions_kw, 0),
    total_utility_power_kw = cumsum(net_additions_kw)
  )

timeline$total_utility_power_gw <- timeline$total_utility_power_kw / 1e6

(f2 <- ggplot(timeline, aes(x = year, y = total_utility_power_gw)) +
    geom_area(fill = "#2E86AB", alpha = 0.25) +
    geom_line(linewidth = 1.1, color = "#2E86AB") +
    labs(
      title    = "Total Utility Power Demand of Operational U.S. Data Centers",
      subtitle = "Based on S&P 451 Database",
      x        = "Year",
      y        = "Power Demand (GW)"
    ) +
    theme_minimal(base_size = 13) +
    scale_x_continuous(
      breaks = seq(
        min(timeline$year, na.rm = TRUE),
        max(timeline$year, na.rm = TRUE),
        by = 5
      )
    ) +
    theme(
      plot.background    = element_rect(fill = "white", color = NA),
      plot.title.position = "plot",
      plot.title         = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle      = element_text(hjust = 0.5),
      axis.title         = element_text(face = "bold")
    ))

ggsave(
  filename = here("senior-essay-sample", "output", "fig2_utility_demand.png"),
  plot     = f2,
  width    = 8,
  height   = 5,
  dpi      = 300
)


# --- 4. Figure 3: County-level distribution of operational centers (2025) ----

dc_2025 <- dc %>%
  filter(
    !is.na(year_operational),
    year_operational <= 2025,
    is.na(year_decommissioned) | year_decommissioned >= 2025
  )

suffix_pat <- "\\s+(County|Parish|Borough|Census Area|City and Borough)$"

center_by_county <- dc_2025 %>%
  transmute(
    county_clean = str_to_lower(str_trim(str_remove(county, suffix_pat))),
    state_full   = str_to_lower(state)
  ) %>%
  filter(!is.na(county_clean), !is.na(state_full)) %>%
  count(state_full, county_clean, name = "n_sites") %>%
  mutate(
    site_category = case_when(
      n_sites >= 40 ~ "≥40",
      n_sites >= 20 ~ "20–39",
      n_sites >= 10 ~ "10–19",
      n_sites >= 5  ~ "5–9",
      n_sites >= 1  ~ "1–4",
      TRUE          ~ "0"
    ),
    site_category = factor(
      site_category,
      levels = c("0", "1–4", "5–9", "10–19", "20–39", "≥40")
    )
  )

counties_sf <- tigris::counties(cb = TRUE, year = 2023, class = "sf") %>%
  st_transform(4326) %>%
  mutate(
    county_clean = str_to_lower(str_trim(str_remove(NAME, suffix_pat))),
    state_full   = str_to_lower(STATE_NAME)
  )

map_sf <- counties_sf %>%
  left_join(center_by_county, by = c("state_full", "county_clean")) %>%
  mutate(
    n_sites = ifelse(is.na(n_sites), 0L, n_sites),
    site_category = ifelse(is.na(site_category), "0", as.character(site_category)),
    site_category = factor(
      site_category,
      levels = c("0", "1–4", "5–9", "10–19", "20–39", "≥40")
    )
  )

# Exclude non-continental U.S. FIPS
exc_fips <- c("02", "15", "72", "78", "60", "66", "69")

map_sf <- map_sf %>%
  filter(!STATEFP %in% exc_fips)

(f3 <- ggplot(map_sf) +
    geom_sf(aes(fill = site_category), color = NA) +
    geom_sf(fill = NA, color = "grey70", linewidth = 0.1) +
    scale_fill_manual(
      values = c(
        "0"     = "grey95",
        "1–4"   = "#c7dbf0",
        "5–9"   = "#8bbce8",
        "10–19" = "#559fd6",
        "20–39" = "#2f79b5",
        "≥40"   = "#0a4a7a"
      ),
      name = NULL,
      drop = FALSE
    ) +
    labs(
      title    = "Distribution of Operational U.S. Data Centers, by County",
      subtitle = "Operational in 2025",
      caption  = "Operational in 2025 = Operational Year ≤ 2025 and Decommission Year ≥ 2025."
    ) +
    coord_sf(xlim = c(-125, -66.5), ylim = c(24.5, 49.5), expand = FALSE) +
    theme_void(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position  = "bottom",
      legend.key.width = unit(1.1, "cm"),
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle    = element_text(hjust = 0.5),
      plot.caption     = element_text(hjust = 0.5)
    ))

ggsave(
  filename = here("senior-essay-sample", "output", "fig3_county_map.png"),
  plot     = f3,
  width    = 8,
  height   = 5,
  dpi      = 300
)


# --- 5. Figure 4: Rack power density distribution ----------------------------

dc_1980 <- dc %>%
  filter(
    !is.na(kW_per_rack),
    kW_per_rack > 0,
    !is.na(year_operational),
    year_operational >= 1980,
    year_operational <= 2025
  )

dc_2020 <- dc %>%
  filter(
    !is.na(kW_per_rack),
    kW_per_rack > 0,
    !is.na(year_operational),
    year_operational >= 2020,
    year_operational <= 2025
  )

f3_1980 <- ggplot(dc_1980, aes(x = kW_per_rack)) +
  geom_histogram(
    bins  = 30,
    fill  = "#2E86AB",
    color = "white",
    alpha = 0.85
  ) +
  geom_vline(xintercept = 20, color = "red3", linetype = "dashed", linewidth = 0.5) +
  annotate(
    "text",
    x = 20.3, y = 250, label = "AI threshold", color = "red3",
    angle = 90, vjust = 1.2, hjust = 0, size = 3.5, fontface = "bold"
  ) +
  xlim(0, 30) +
  labs(
    subtitle = "Operational between 1980–2025",
    x        = "kW per Rack",
    y        = "Count of Data Centers"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.minor   = element_blank(),
    plot.title.position = "plot",
    plot.subtitle      = element_text(hjust = 0.5),
    axis.title         = element_text(face = "bold")
  )

f3_2020 <- ggplot(dc_2020, aes(x = kW_per_rack)) +
  geom_histogram(
    bins  = 30,
    fill  = "#2E86AB",
    color = "white",
    alpha = 0.85
  ) +
  geom_vline(xintercept = 20, color = "red3", linetype = "dashed", linewidth = 0.5) +
  annotate(
    "text",
    x = 20.3, y = 75, label = "AI threshold", color = "red3",
    angle = 90, vjust = 1.2, hjust = 0, size = 3.5, fontface = "bold"
  ) +
  xlim(0, 30) +
  labs(
    subtitle = "Operational between 2020–2025",
    x        = "kW per Rack",
    y        = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.minor   = element_blank(),
    plot.title.position = "plot",
    plot.subtitle      = element_text(hjust = 0.5),
    axis.title         = element_text(face = "bold")
  )

(f4 <- (f3_1980 | f3_2020) +
    plot_annotation(
      title = "Distribution of Data Center Rack Power Density",
      theme = theme(
        plot.title     = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.caption   = element_text(hjust = 0.5, size = 10),
        plot.background = element_rect(fill = "white", color = NA)
      )
    ))

ggsave(
  filename = here("senior-essay-sample", "output", "fig4_rack_power_density.png"),
  plot     = f4,
  width    = 8,
  height   = 5,
  dpi      = 300
)