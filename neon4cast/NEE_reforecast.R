# install.packages("ranger")

library(tidyverse)
library(tidymodels)
library(lubridate)
library(duckdbfs)
library(neon4cast)
library(dplyr)
library(scoringRules)
library(arrow)

# set up
ref_dates <- seq.Date(
  from = as.Date("2024-01-01"),
  to   = as.Date("2024-12-31"),
  by   = "7 days"
)

my_var <- "nee"
n_boot <- 5
my_model_id <- "rf_hybrid_model"

# target variable
target_url <- "https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/terrestrial_daily-targets.csv.gz"

targets <- read_csv(target_url, show_col_types = FALSE) |> 
  filter(variable == my_var) |> 
  mutate(date = as_date(datetime)) |> 
  arrange(site_id, date) |>
  group_by(site_id) |>
  mutate(nee_lag1 = lag(observation, 1)) |>
  ungroup() |>
  drop_na(observation)

# baseline models
clim_ds <- arrow::open_dataset(
  "s3://anonymous@bio230014-bucket01/challenges/forecasts/bundled-parquet/project_id=neon4cast/duration=P1D/variable=nee/model_id=climatology?endpoint_override=sdsc.osn.xsede.org"
)

pers_ds <- arrow::open_dataset(
  "s3://anonymous@bio230014-bucket01/challenges/forecasts/bundled-parquet/project_id=neon4cast/duration=P1D/variable=nee/model_id=persistenceRW?endpoint_override=sdsc.osn.xsede.org"
)

# CRPS function
crps_tidy_ensemble <- function(prediction, observation) {
  pred <- prediction[!is.na(prediction)]
  obs <- observation[1]
  if (length(pred) == 0 || is.na(obs)) return(NA_real_)
  crps_sample(y = obs, dat = pred)
}

# store all outputs
all_forecasts <- list()

# Target 1: reforecast loop
noaa_raw_all <- neon4cast::noaa_stage3() |>
  filter(variable %in% c("air_temperature", 
                         "surface_downwelling_shortwave_flux_in_air",
                         "relative_humidity")) |> 
  filter(datetime >= as_datetime("2023-01-01"),
         datetime <= as_datetime("2025-02-05")) |> 
  collect()

for (i in seq_along(ref_dates)) {
  
  ref_date <- as_datetime(ref_dates[i])
  forecast_end_date <- ref_date + days(35)
  
  message("Running reference date: ", as.Date(ref_date))
  
  # predicting variable
  noaa_raw <- noaa_raw_all |>
    filter(datetime <= forecast_end_date)
  
  noaa_daily <- noaa_raw |>
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
  
  # data division
  train_weather <- noaa_daily |>
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
    drop_na(observation, temp_c, surface_downwelling_shortwave_flux_in_air, vpd, sin_doy, cos_doy, nee_lag1)
  
  last_nee <- targets |>
    filter(date < as_date(ref_date)) |>
    group_by(site_id) |>
    slice_max(date, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(site_id, nee_lag1 = observation)
  
  forecast_data <- noaa_daily |>
    filter(date >= as_date(ref_date),
           date < as_date(forecast_end_date)) |>
    left_join(last_nee, by = "site_id") |>
    drop_na(temp_c, surface_downwelling_shortwave_flux_in_air, vpd, sin_doy, cos_doy, nee_lag1)
  
  # RF model
  forecast_ens <- list()
  
  for (b in 1:n_boot) {
    boot_train <- train_data |> 
      slice_sample(prop = 1, replace = TRUE)
    
    rf_fit_b <- workflow() |>
      add_recipe(recipe(
        observation ~ temp_c + surface_downwelling_shortwave_flux_in_air + vpd + sin_doy + cos_doy + nee_lag1,
        data = boot_train
      )) |>
      add_model(
        rand_forest(trees = 200) |>
          set_engine("ranger") |>
          set_mode("regression")
      ) |>
      fit(data = boot_train)
    
    train_preds_b <- predict(rf_fit_b, new_data = boot_train)
    resid_sd <- sd(boot_train$observation - train_preds_b$.pred, na.rm = TRUE)
    
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
  
  final_forecast_i <- forecast_ensemble |>
    mutate(
      parameter = paste0("m", parameter, "_b", boot_id),
      datetime = as_datetime(date),
      reference_datetime = ref_date,
      family = "ensemble",
      duration = "P1D",
      variable = my_var,
      model_id = my_model_id,
      project_id = "neon4cast"
    ) |>
    select(
      project_id, model_id, datetime, reference_datetime,
      duration, site_id, family, parameter, variable, prediction
    )
  
  all_forecasts[[i]] <- final_forecast_i
}

final_forecast_all <- bind_rows(all_forecasts)

write_csv(final_forecast_all, "terrestrial_daily-reforecast_all_2024.csv")

# Target 2: evaluation

# prepare target data for evaluation
targets_eval <- targets |>
  mutate(datetime = as.Date(datetime)) |>
  select(datetime, site_id, variable, observation)

# my model CRPS
my_forecasts <- final_forecast_all |>
  mutate(
    
    
    datetime = as.Date(datetime),
    reference_datetime = as.Date(reference_datetime)
  )

# Figure 1: Prediction time series for each site across the reforecast period
site_prediction_summary <- my_forecasts |>
  group_by(site_id, datetime, variable) |>
  summarize(
    pred_q05 = quantile(prediction, 0.05, na.rm = TRUE),
    pred_q25 = quantile(prediction, 0.25, na.rm = TRUE),
    pred_q50 = quantile(prediction, 0.50, na.rm = TRUE),
    pred_q75 = quantile(prediction, 0.75, na.rm = TRUE),
    pred_q95 = quantile(prediction, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(targets_eval, by = c("datetime", "site_id", "variable"))

p_site_timeseries <- ggplot(site_prediction_summary,
                            aes(x = datetime)) +
  geom_ribbon(aes(ymin = pred_q05, ymax = pred_q95),
              alpha = 0.15) +
  geom_ribbon(aes(ymin = pred_q25, ymax = pred_q75),
              alpha = 0.25) +
  geom_line(aes(y = pred_q50), linewidth = 0.4) +
  geom_point(aes(y = observation), size = 0.4, alpha = 0.5, na.rm = TRUE) +
  facet_wrap(~ site_id, scales = "free_y") +
  labs(
    x = "Date",
    y = "NEE",
    title = "Site-level reforecast prediction time series",
    subtitle = "Line = median prediction; dark band = 50% interval; light band = 90% interval; points = observations"
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 7))

print(p_site_timeseries)
ggsave("site_prediction_timeseries_2024.png",
       p_site_timeseries, width = 12, height = 9)

write_csv(site_prediction_summary, "site_prediction_summary_2024.csv")

my_crps <- my_forecasts |>
  left_join(targets_eval, by = c("datetime", "site_id", "variable")) |>
  drop_na(observation) |>
  mutate(lead_time = as.numeric(datetime - reference_datetime)) |>
  group_by(model_id, reference_datetime, site_id, datetime, lead_time) |>
  summarize(crps = crps_tidy_ensemble(prediction, observation), .groups = "drop")

# baseline forecasts
score_ds <- arrow::open_dataset(
  "s3://anonymous@bio230014-bucket01/challenges/scores/bundled-parquet/project_id=neon4cast/duration=P1D/variable=nee?endpoint_override=sdsc.osn.xsede.org"
)

baseline_crps <- score_ds |>
  filter(
    model_id %in% c("climatology", "persistenceRW"),
    reference_datetime %in% as.Date(ref_dates))|>
  select(model_id, reference_datetime, site_id, datetime, crps) |>
  collect() |>
  mutate(
    datetime = as.Date(datetime),
    reference_datetime = as.Date(reference_datetime),
    lead_time = as.numeric(datetime - reference_datetime)
  )

all_crps <- bind_rows(
  my_crps,
  baseline_crps
)

crps_summary <- all_crps |>
  group_by(model_id, lead_time) |>
  summarize(mean_crps = mean(crps, na.rm = TRUE), .groups = "drop")

# Figure 2: CRPS vs forecast horizon
p <- ggplot(crps_summary,
            aes(x = lead_time, y = mean_crps, color = model_id)) +
  geom_line(linewidth = 1) +
  labs(x = "Forecast horizon (days)",
       y = "Mean CRPS",
       color = "Model",
       title = "Reforecast performance comparison") +
  theme_bw()

print(p)

ggsave("crps_comparison_2024.png", p, width = 7, height = 5)

# Figure 3: Seasonal performance analysis
crps_by_month <- all_crps |>
  mutate(
    month_num = month(datetime),
    month = month(datetime, label = TRUE, abbr = TRUE)
  ) |>
  group_by(model_id, month_num, month) |>
  summarize(mean_crps = mean(crps, na.rm = TRUE), .groups = "drop")

p_crps_month <- ggplot(crps_by_month,
                       aes(x = month_num, y = mean_crps, color = model_id)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  labs(
    x = "Target month",
    y = "Mean CRPS",
    color = "Model",
    title = "Seasonal pattern in reforecast performance"
  ) +
  theme_bw()

print(p_crps_month)
ggsave("crps_vs_month_2024.png", p_crps_month, width = 8, height = 5)

write_csv(crps_by_month, "crps_by_month_2024.csv")

# Figure 4: Example forecast trajectories with uncertainty ribbons
# This is a case-study figure, not a global performance metric. It uses one representative reference date and several sites to show how uncertainty unfolds across the 35-day forecast horizon.
case_ref_date <- as.Date("2024-07-01")
if (!case_ref_date %in% as.Date(ref_dates)) {
  case_ref_date <- as.Date(ref_dates)[which.min(abs(as.Date(ref_dates) - case_ref_date))]
}

case_sites <- my_forecasts |>
  filter(reference_datetime == case_ref_date) |>
  distinct(site_id, datetime) |>
  count(site_id, name = "n_dates") |>
  arrange(desc(n_dates)) |>
  slice_head(n = 3) |>
  pull(site_id)

ribbon_case <- my_forecasts |>
  filter(reference_datetime == case_ref_date,
         site_id %in% case_sites) |>
  group_by(site_id, reference_datetime, datetime) |>
  summarize(
    pred_q05 = quantile(prediction, 0.05, na.rm = TRUE),
    pred_q25 = quantile(prediction, 0.25, na.rm = TRUE),
    pred_q50 = quantile(prediction, 0.50, na.rm = TRUE),
    pred_q75 = quantile(prediction, 0.75, na.rm = TRUE),
    pred_q95 = quantile(prediction, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(targets_eval, by = c("datetime", "site_id")) |>
  mutate(lead_time = as.numeric(datetime - reference_datetime))

p_ribbon_case <- ggplot(ribbon_case,
                        aes(x = lead_time)) +
  geom_ribbon(aes(ymin = pred_q05, ymax = pred_q95),
              alpha = 0.15) +
  geom_ribbon(aes(ymin = pred_q25, ymax = pred_q75),
              alpha = 0.25) +
  geom_line(aes(y = pred_q50), linewidth = 0.7) +
  geom_point(aes(y = observation), size = 1.2, alpha = 0.7, na.rm = TRUE) +
  facet_wrap(~ site_id, scales = "free_y") +
  labs(
    x = "Forecast horizon (days)",
    y = "NEE",
    title = paste0("Example forecast uncertainty trajectories: reference date ", case_ref_date),
    subtitle = "Line = median prediction; dark band = 50% interval; light band = 90% interval; points = observations"
  ) +
  theme_bw()

print(p_ribbon_case)
ggsave("forecast_uncertainty_ribbon_typical_sites_2024.png",
       p_ribbon_case, width = 9, height = 5)

write_csv(ribbon_case, "forecast_uncertainty_ribbon_typical_sites_2024.csv")
