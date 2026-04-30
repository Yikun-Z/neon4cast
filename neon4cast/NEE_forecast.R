#install.packages("ranger")
#install.packages("tidyverse")
#install.packages("tidymodels")
#install.packages("lubridate")
#install.packages("neon4cast")
#install.packages("scoringRules")
#install.packages("arrow")

library(tidyverse)
library(tidymodels)
library(lubridate)
library(neon4cast)
library(scoringRules)
library(arrow)

# ============================================================
# 0. Site passed from GitHub Actions matrix
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No site_id was provided. Please pass one site_id from GitHub Actions matrix.")
}

site_to_run <- args[1]
message("Running forecast for site: ", site_to_run)

# ============================================================
# 1. Set up
# ============================================================

# Use the latest available Stage 2 forecast date
ref_date <- as_date(now(tzone = "UTC")) - days(2)
forecast_end_date <- ref_date + days(35)

my_var <- "nee"
n_boot <- 5
my_model_id <- "rf_hybrid_model"

# ============================================================
# 2. Target variable
# ============================================================

target_url <- "https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/terrestrial_daily-targets.csv.gz"

targets <- read_csv(target_url, show_col_types = FALSE) |> 
  filter(
    variable == my_var,
    site_id == site_to_run
  ) |> 
  mutate(date = as_date(datetime)) |> 
  arrange(site_id, date) |>
  group_by(site_id) |>
  mutate(nee_lag1 = lag(observation, 1)) |>
  ungroup() |>
  drop_na(observation)

if (nrow(targets) == 0) {
  stop("No target data available for site: ", site_to_run)
}

# ============================================================
# 3. NOAA Stage 3 drivers for training
# ============================================================

noaa_train_raw <- neon4cast::noaa_stage3() |>
  filter(site_id == site_to_run) |>
  filter(variable %in% c(
    "air_temperature",
    "surface_downwelling_shortwave_flux_in_air",
    "relative_humidity"
  )) |>
  filter(
    datetime >= as_datetime("2023-01-01"),
    datetime <= as_datetime(forecast_end_date)
  ) |>
  collect()

if (nrow(noaa_train_raw) == 0) {
  stop("No NOAA Stage 3 training drivers available for site: ", site_to_run)
}

noaa_train_daily <- noaa_train_raw |>
  mutate(date = as_date(datetime)) |>
  group_by(site_id, date, parameter, variable) |>
  summarize(
    daily_mean = mean(prediction, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = variable,
    values_from = daily_mean
  ) |>
  mutate(
    temp_c = air_temperature - 273.15,
    es = 0.6108 * exp((17.27 * temp_c) / (temp_c + 237.3)),
    vpd = es * (1 - relative_humidity / 100),
    doy = yday(date),
    sin_doy = sin(2 * pi * doy / 365),
    cos_doy = cos(2 * pi * doy / 365)
  )

# Training weather: average NOAA ensemble for historical training
train_weather <- noaa_train_daily |>
  group_by(site_id, date) |>
  summarize(
    temp_c = mean(temp_c, na.rm = TRUE),
    surface_downwelling_shortwave_flux_in_air = mean(surface_downwelling_shortwave_flux_in_air, na.rm = TRUE),
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
  stop("No complete training data available for site: ", site_to_run)
}

# ============================================================
# 4. Latest observed NEE before forecast date
# ============================================================

last_nee <- targets |>
  filter(date < as_date(ref_date)) |>
  group_by(site_id) |>
  slice_max(date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(site_id, nee_lag1 = observation)

if (nrow(last_nee) == 0) {
  stop("No latest NEE observation available for site: ", site_to_run)
}

# ============================================================
# 5. NOAA Stage 2 drivers for forecast
# ============================================================

noaa_forecast_raw <- neon4cast::noaa_stage2(
  start_date = as.character(ref_date)
) |>
  filter(site_id == site_to_run) |>
  filter(variable %in% c(
    "air_temperature",
    "surface_downwelling_shortwave_flux_in_air",
    "relative_humidity"
  )) |>
  collect()

if (nrow(noaa_forecast_raw) == 0) {
  stop("No NOAA Stage 2 forecast drivers available for site: ", site_to_run)
}

noaa_forecast_daily <- noaa_forecast_raw |>
  mutate(date = as_date(datetime)) |>
  group_by(site_id, date, parameter, variable) |>
  summarize(
    daily_mean = mean(prediction, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = variable,
    values_from = daily_mean
  ) |>
  mutate(
    temp_c = air_temperature - 273.15,
    es = 0.6108 * exp((17.27 * temp_c) / (temp_c + 237.3)),
    vpd = es * (1 - relative_humidity / 100),
    doy = yday(date),
    sin_doy = sin(2 * pi * doy / 365),
    cos_doy = cos(2 * pi * doy / 365)
  )

forecast_data <- noaa_forecast_daily |>
  filter(
    date >= as_date(ref_date),
    date < as_date(forecast_end_date)
  ) |>
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
  stop("No complete forecast data available for site: ", site_to_run)
}

# ============================================================
# 6. RF ensemble forecast
# ============================================================

forecast_ens <- list()

for (b in 1:n_boot) {
  
  message("Bootstrap member: ", b, " / ", n_boot)
  
  boot_train <- train_data |> 
    slice_sample(prop = 1, replace = TRUE)
  
  rf_fit_b <- workflow() |>
    add_recipe(recipe(
      observation ~ temp_c +
        surface_downwelling_shortwave_flux_in_air +
        vpd +
        sin_doy +
        cos_doy +
        nee_lag1,
      data = boot_train
    )) |>
    add_model(
      rand_forest(trees = 200) |>
        set_engine("ranger") |>
        set_mode("regression")
    ) |>
    fit(data = boot_train)
  
  train_preds_b <- predict(rf_fit_b, new_data = boot_train)
  
  resid_sd <- sd(
    boot_train$observation - train_preds_b$.pred,
    na.rm = TRUE
  )
  
  if (is.na(resid_sd) || resid_sd == 0) {
    resid_sd <- 1e-6
  }
  
  preds_b <- predict(rf_fit_b, new_data = forecast_data)
  
  forecast_ens[[b]] <- forecast_data |>
    bind_cols(preds_b) |>
    mutate(
      prediction = .pred + rnorm(n(), 0, resid_sd),
      boot_id = b
    ) |>
    select(site_id, date, parameter, boot_id, prediction)
}

forecast_ensemble <- bind_rows(forecast_ens)

# ============================================================
# 7. EFI-format final forecast
# ============================================================

final_forecast <- forecast_ensemble |>
  mutate(
    parameter = paste0("m", parameter, "_b", boot_id),
    datetime = as_datetime(date),
    reference_datetime = as_datetime(ref_date),
    family = "ensemble",
    duration = "P1D",
    variable = my_var,
    model_id = my_model_id,
    project_id = "neon4cast"
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

# ============================================================
# 8. Output file
# ============================================================

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
