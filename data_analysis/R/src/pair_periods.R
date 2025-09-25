library(tidyverse)

pair_periods <- function(df,
                         hist = "1981_2005",
                         fut  = "2076_2100",
                         vars,
                         by   = c("model","label"),
                         keep = c("color","fill","border","shape","linetype"),
                         agg  = \(x) mean(x, na.rm = TRUE)) {
  
  hist_tbl <- df |>
    filter(PERIOD == hist) |>
    group_by(across(all_of(by))) |>
    summarise(
      across(all_of(vars), agg),
      across(any_of(keep), ~ dplyr::first(na.omit(.x))),
      .groups = "drop"
    ) |>
    rename_with(~ paste0(.x, "_hist"), all_of(vars))    # only metrics get _hist
  
  fut_tbl <- df |>
    filter(PERIOD == fut) |>
    group_by(across(all_of(by))) |>
    summarise(
      across(all_of(vars), agg),
      .groups = "drop"
    ) |>
    rename_with(~ paste0(.x, "_fut"), all_of(vars)) |>  # only metrics get _fut
    select(all_of(by), ends_with("_fut"))               # DROP styling cols here
  
  inner_join(hist_tbl, fut_tbl, by = by)
}

build_scatter_data <- function(paired_df, var,
                               by = c("model","label"),
                               keep = c("color","fill","border","shape","linetype")) {
  paired_df |>
    transmute(
      model = interaction(.data[[by[1]]], drop = TRUE),
      label = .data[[by[2]]],
      across(any_of(keep)),
      x = .data[[paste0(var, "_hist")]],
      y = .data[[paste0(var, "_fut")]]
    )
}
