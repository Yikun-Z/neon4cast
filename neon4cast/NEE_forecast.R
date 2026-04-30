# Packages ----------------------------------------------------
library(tidyverse)
library(lubridate)
library(neon4cast)
library(scoringRules)
library(arrow)
library(workflows)
library(recipes)
library(parsnip)
library(ranger)

# Setup -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No site_id was provided.")

site_to_run <- args[1]
message("Running forecast for site: ", site_to_run)

ref_date <- as_date(now(tzone = "UTC")) - days(2)
forecast_end_date <- ref_date + days(35)

my_var <- "nee"
n_boot <- 5
my_model_id <- "rf_hybrid_model"

drivers <- c(
  "air_temperature",
  "surface_downwelling_shortwave_flux_in_air",
  "relative_humidity"
)

# Helper function ---------------------------------------------

make_daily_weather <- function(x) {
  x |>
    mutate(date = as_date(datetime)) |>
    group_by(site_id, date, parameter, variable) |>
    summarize(daily_mean = mean(prediction, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = variable, values_from = daily_mean) |>
    mutate(
      temp_c = air_temperature - 273.15,
      es = 0.6108 * exp((17.27 * temp_c) / (temp_c + 237.3)),
      vpd = es * (1 - relative_humidity / 100),
      doy = yday(date),
      sin_doy = sin(2 * pi * doy / 365),
      cos_doy = cos(2 * pi * doy / 365)
    )
}

# Targets -----------------------------------------------------

target_url <- "https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/terrestrial_daily-targets.csv.gz"

targets <- read_csv(target_url, show_col_types = FALSE) |>
  filter(variable == my_var, site_id == site_to_run) |>
  mutate(date = as_date(datetime)) |>
  arrange(site_id, date) |>
  group_by(site_id) |>
  mutate(nee_lag1 = lag(observation, 1)) |>
  ungroup() |>
  drop_na(observation)

if (nrow(targets) == 0) stop("No target data for site: ", site_to_run)

# Training drivers --------------------------------------------

noaa_train_raw <- neon4cast::noaa_stage3() |>
  filter(
    site_id == site_to_run,
    variable %in% drivers,
    datetime >= as_datetime("2023-01-01"),
    datetime <= as_datetime(forecast_end_date)
  ) |>
  collect()

if (nrow(noaa_train_raw) == 0) {
  stop("No NOAA Stage 3 training drivers for site: ", site_to_run)
}

train_weather <- make_daily_weather(noaa_train_raw) |>
  group_by(site_id, date) |>
  summarize(
    temp_c = mean(temp_c, na.rm = TRUE),
    surface_downwelling_shortwave_flux_in_air =
      mean(surface_downwelling_shortwave_flux_in_air, na.rm = TRUE),
    vpd = mean(vpd, na.rm = TRUE),
    sin_doy = first(sin_doy),
    cos_doy = first(cos_doy),
    .groups = "drop"
  )

train_data <- targets |>
  filter(date < as_date(ref_date)) |>
  left_join(train_weather, by = c("site_id", "date")) |>
  drop_na(
    observation,
    nee_lag1,
    temp_c,
    surface_downwelling_shortwave_flux_in_air,
    vpd,
    sin_doy,
    cos_doy
  )

if (nrow(train_data) == 0) {
  stop("No complete training data for site: ", site_to_run)
}

last_nee <- targets |>
  filter(date < as_date(ref_date)) |>
  slice_max(date, n = 1, with_ties = FALSE) |>
  select(site_id, nee_lag1 = observation)

if (nrow(last_nee) == 0) {
  stop("No latest NEE observation for site: ", site_to_run)
}

# Forecast drivers --------------------------------------------

noaa_forecast_raw <- neon4cast::noaa_stage2(
  start_date = as.character(ref_date)
) |>
  filter(site_id == site_to_run, variable %in% drivers) |>
  collect()

if (nrow(noaa_forecast_raw) == 0) {
  stop("No NOAA Stage 2 forecast drivers for site: ", site_to_run)
}

forecast_data <- make_daily_weather(noaa_forecast_raw) |>
  filter(date >= as_date(ref_date), date < as_date(forecast_end_date)) |>
  left_join(last_nee, by = "site_id") |>
  drop_na(
    temp_c,
    surface_downwelling_shortwave_flux_in_air,
    vpd,
    sin_doy,
    cos_doy,
    nee_lag1
  )

if (nrow(forecast_data) == 0) {
  stop("No complete forecast data for site: ", site_to_run)
}

# RF bootstrap ensemble ---------------------------------------

forecast_ens <- vector("list", n_boot)

for (b in seq_len(n_boot)) {
  message("Bootstrap member: ", b, " / ", n_boot)
  
  boot_train <- train_data |>
    slice_sample(prop = 1, replace = TRUE)
  
  rf_fit <- workflow() |>
    add_recipe(
      recipe(
        observation ~ temp_c +
          surface_downwelling_shortwave_flux_in_air +
          vpd +
          sin_doy +
          cos_doy +
          nee_lag1,
        data = boot_train
      )
    ) |>
    add_model(
      rand_forest(trees = 200) |>
        set_engine("ranger") |>
        set_mode("regression")
    ) |>
    fit(data = boot_train)
  
  train_preds <- predict(rf_fit, new_data = boot_train)
  resid_sd <- sd(boot_train$observation - train_preds$.pred, na.rm = TRUE)
  if (is.na(resid_sd) || resid_sd == 0) resid_sd <- 1e-6
  
  forecast_ens[[b]] <- forecast_data |>
    bind_cols(predict(rf_fit, new_data = forecast_data)) |>
    mutate(
      prediction = .pred + rnorm(n(), 0, resid_sd),
      boot_id = b
    ) |>
    select(site_id, date, parameter, boot_id, prediction)
}

# Format output -----------------------------------------------

final_forecast <- bind_rows(forecast_ens) |>
  mutate(
    project_id = "neon4cast",
    model_id = my_model_id,
    datetime = as_datetime(date),
    reference_datetime = as_datetime(ref_date),
    duration = "P1D",
    family = "ensemble",
    parameter = paste0("m", parameter, "_b", boot_id),
    variable = my_var
  ) |>
  select(
    project_id,
    model_id,
    datetime,
    reference_datetime,
    duration,
    site_id,
    family,
    parameter,
    variable,
    prediction
  )

if (nrow(final_forecast) == 0) {
  stop("Final forecast is empty for site: ", site_to_run)
}

out_file <- paste0(
  "terrestrial_daily-",
  as_date(ref_date),
  "-",
  my_model_id,
  "-",
  site_to_run,
  ".csv"
)

write_csv(final_forecast, out_file)
message("Forecast written to: ", out_file)
