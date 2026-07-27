pacman::p_load(dplyr, tidyr, purrr, stringr, sf, leaflet, 
               car, here, ggplot2, rlang, readr, mapview,
               readxl, raster, caret, gt, ranger, tibble)



################################################################################
# Basic descriptive statistics helpers
################################################################################

#' Geometric standard deviation
#'
#' @param col Numeric vector to compute the geometric standard deviation for.
#' @param na_rm Logical, whether `NA` values should be stripped. 
#'
#' @return A single numeric value: the geometric standard deviation of `col`.
#'
#' @examples
#' \dontrun{
#' GSD1(c(1, 2, 3, NA, 5))
#' }
GSD1 <- function(col, na_rm = TRUE) {
  GSD_col <- exp(sd(log(col), na.rm = na_rm))
  GSD_col
}

#' Geometric mean
#'
#' @param col Numeric vector to compute the geometric mean for.
#' @param na_rm Logical, whether `NA` values should be stripped. 
#'
#' @return A single numeric value: the geometric mean of `col`.
GM1 <- function(col, na_rm = TRUE) {
  GM_col <- exp(mean(log(col), na.rm = na_rm)) 
  GM_col
}

#' Coefficient of variation
#'
#' @param col Numeric vector.
#' @param na_rm Logical, whether `NA` values should be stripped. 
#'
#' @return A single numeric value: `sd(col) / mean(col)`.
CV1 <- function(col, na_rm = TRUE) {
  cv <- sd(col, na.rm = na_rm) / mean(col, na.rm = na_rm)
  cv
}

#' Standard error of the mean
#'
#' @param col Numeric vector.
#' @param na.rm Logical, whether to drop `NA` values before computing the
#'   standard error.
#'
#' @return A single numeric value: the standard error of `col`.
stderr <- function(col, na.rm = FALSE) {
  if (na.rm) col <- na.omit(col)
  sqrt(var(col) / length(col))
}

#' Significance stars for a p-value
#'
#' @param x A single numeric p-value (or `NA`/`NULL`/`" "`).
#'
#' @return A character string: `"***"` (p < 0.001), `"**"` (p < 0.01),
#'   `"*"` (p < 0.05), `"."` (p < 0.1), or `" "` otherwise (including for
#'   `NA`, `NULL`, or an already-blank input).
sig_star <- function(x) {
  if (x < 0.001) {
    return ("***")
  } else if (x < 0.01) {
    return ("**")
  } else if (x < 0.05) {
    return ("*")
  } else if (x < 0.1) {
    return (".")
  } else if(is.na(x)) {
    return (" ")
  } else if(is.null(x)) {
    return (" ")
  } else if(x == " ") {
    return (" ")
  } else {
    return (" ")
  }
}

#' Multiple assignment operator
#'
#' Assigns multiple values returned in a list (`rhs`) to multiple variable
#' names given on the left-hand side, e.g.
#' `c(a, b, c) := list(1, 2, 3)`. Adapted from
#' <https://stackoverflow.com/a/1829651>.
#'
#' @param lhs An unquoted call of the form `c(var1, var2, ...)` naming the
#'   variables to assign to.
#' @param rhs A list (or function/formula) whose elements are assigned, in
#'   order, to the names in `lhs`.
#'
#' @return Invisibly `NULL`; used for its side effect of assigning variables
#'   into the calling frame.
':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir = frame)
    return(invisible(NULL)) 
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in seq_len(length(lhs)))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir = frame)
  return(invisible(NULL)) 
}



################################################################################
# Functions for model building
################################################################################

#' Extract coefficient-level statistics and direction-of-effect check
#'
#' Extracts slope, standard error, t-value, p-value and significance stars for
#' every predictor in a fitted `lm` model (excluding the intercept), and joins
#' on the desired direction-of-effect table to flag whether each predictor's
#' sign matches expectation.
#'
#' @param my_model A fitted `lm` object.
#' @param sig_star Function used to convert a p-value to a significance-star
#'   string; normally `sig_star` itself.
#' @param direction_of_effect_table Data frame with columns `param`, `sign`,
#'   and `val` describing the expected direction of effect for each family of
#'   predictors (see `add_sign_check()`).
#'
#' @return A data frame with one row per predictor (intercept excluded) and
#'   columns `value`, `slope`, `stde`, `tval`, `prob`, `sig`, `param`, `sign`,
#'   `val`, and `obj` (logical: is the direction of effect satisfied?).
#'
extract_model_data <- function(my_model, 
                               sig_star, 
                               direction_of_effect_table) {
  if (!inherits(my_model, "lm")) {
    stop("`my_model` must be a fitted lm object.", call. = FALSE)
  }
  if (!is.function(sig_star)) {
    stop("`sig_star` must be a function.", call. = FALSE)
  }
  required_cols <- c("param", "sign", "val")
  if (!all(required_cols %in% names(direction_of_effect_table))) {
    stop(
      "`direction_of_effect_table` must contain columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }
  coef_table <- summary(my_model)$coefficients
  data_extracted <- tibble::tibble(
    value = rownames(coef_table),
    slope = format(coef_table[, 1], scientific = TRUE),
    stde  = round(coef_table[, 2], 4),
    tval  = round(coef_table[, 3], 4),
    prob  = round(coef_table[, 4], 4),
    sig   = vapply(coef_table[, 4], sig_star, character(1))
    ) %>% 
    filter(value != "(Intercept)")
  
  data_extracted <- add_sign_check(data_extracted, direction_of_effect_table)
}

#' Check whether a slope's sign matches the expected direction of effect
#'
#' @param sign Character string naming a comparison function, e.g. `"<"` or
#'   `">"`, taken from the `sign` column of the direction-of-effect table.
#' @param slope Numeric coefficient value to test.
#' @param val Numeric value to compare `slope` against (typically `0`), or
#'   `NA` if no direction is enforced for this parameter.
#'
#' @return `TRUE` if `val` is `NA` (no constraint) or if
#'   `do.call(sign, list(slope, val))` is `TRUE`; otherwise `FALSE`.
#'
check_sign_from_table <- function(sign, slope, val) {
  valid_signs <- c("<", "<=", ">", ">=", "is.na")
  if (is.na(val)) {
    return(TRUE)
  }
  if (!sign %in% valid_signs) {
    stop(
      "`sign` must be one of: ",
      paste(valid_signs, collapse = ", "),
      call. = FALSE
    )
  }
  do.call(sign, list(slope, val))
}

#' Join direction-of-effect expectations onto extracted model coefficients
#'
#' Derives a `param` key from each coefficient's variable name (by stripping
#' the `_buffer...` suffix), joins the expected sign/value from
#' `direction_of_effect_table`, and evaluates whether each coefficient's sign
#' is as expected.
#'
#' @param data_extracted Data frame of extracted coefficients (as produced by
#'   `extract_model_data()`), must contain a `value` column with the raw
#'   predictor name (e.g. `"tree_cover_buffer_5000"`).
#' @param direction_of_effect_table Data frame with columns `param`, `sign`,
#'   `val` giving the expected direction of effect per base parameter.
#'
#' @return `data_extracted` with additional columns `param`, `sign`, `val`,
#'   and `obj` (logical, whether the direction of effect is satisfied).
#'
add_sign_check <- function(data_extracted, 
                           direction_of_effect_table) {
  required_cols <- c("value", "slope")
  missing_cols <- setdiff(required_cols, names(data_extracted))
  if (length(missing_cols) > 0) {
    stop(
      "`data_extracted` is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  required_doe_cols <- c("param", "sign", "val")
  missing_doe_cols <- setdiff(required_doe_cols, names(direction_of_effect_table))
  if (length(missing_doe_cols) > 0) {
    stop(
      "`direction_of_effect_table` is missing required columns: ",
      paste(missing_doe_cols, collapse = ", "),
      call. = FALSE
    )
  }
  data_extracted$param <- sub("_buffer.*", "\\1", data_extracted$value)
  data_extracted <- left_join(data_extracted, direction_of_effect_table, by = "param")
  for(i in seq_len(nrow(data_extracted))) {
    data_extracted$obj[i] <- check_sign_from_table(data_extracted$sign[i], 
                                                   as.numeric(as.character(data_extracted$slope[i])), 
                                                   as.numeric(as.character(data_extracted$val[i])))
  }
  data_extracted
}

#' Build a linear model from a set of selected columns
#'
#' Constructs a linear model using all predictor columns in `model_data`
#' except `CODE`, the response variable, and columns whose names contain
#' `"predicted"`.
#'
#' @param model_data Data frame containing `CODE`, the response variable, and
#'   the predictor columns to include.
#' @param response_variable Character string naming the response column
#'   (e.g. `"PM2.5"`).
#'
#' @return A fitted `lm` object.
create_model <- function(model_data, response_variable) {
  if (!is.data.frame(model_data)) {
    stop(
      "`model_data` must be a data.frame.",
      call. = FALSE
    )
  }
  if (!response_variable %in% names(model_data)) {
    stop(
      "`response_variable` is not a column in `model_data`.",
      call. = FALSE
    )
  }
  new_data <- model_data %>% 
    dplyr::select(everything(), 
                  -all_of(response_variable),
                  -"CODE", 
                  -contains(c("predicted")))
  others <- names(new_data)
  if (length(others) == 0) {
    stop(
      "No predictor variables remain after excluding CODE, the response variable, and predicted columns.",
      call. = FALSE
    )
  }
  if(ncol(new_data) == 1){  
    eqtn <- as.formula(paste(response_variable, others, sep = " ~ "))
  } else { 
    eqtn <- as.formula(paste(response_variable, paste(others, collapse = "+"), sep = " ~ "))
  }
  my_model <- lm(eqtn, data = model_data)
  my_model
}

#' Drop predictors with p-value above a threshold
#'
#' @param my_model A fitted `lm` object.
#' @param data_with_selected_variables Data frame containing `CODE`, the
#'   response variable, and every predictor currently in `my_model`.
#' @param p_threshold Numeric p-value threshold; predictors with p-value strictly
#'   greater than `p_threshold` are dropped.
#'   
#' @details
#' Predictor removal is based on the full-precision p-values returned by
#' `summary(my_model)$coefficients`. P-values are not rounded before
#' comparison with the specified threshold.
#'
#' @return `data_with_selected_variables` with the offending predictor
#'   columns removed (columns only - rows/observations are untouched).
#'
remove_p_value <- function(my_model, 
                           data_with_selected_variables, 
                           p_threshold) {
  if (!inherits(my_model, "lm")) {
    stop("`my_model` must be a fitted lm object.", call. = FALSE)
  }
  if (!is.data.frame(data_with_selected_variables)) {
    stop(
      "`data_with_selected_variables` must be a data frame.",
      call. = FALSE
    )
  }
  if (!is.numeric(p_threshold) || length(p_threshold) != 1 || p_threshold < 0 || p_threshold > 1) {
    stop(
      "`p_threshold ` must be a single numeric value between 0 and 1.",
      call. = FALSE
    )
  }
  p_values <- summary(my_model)$coefficients[, 4]
  p_remove <- names(p_values)[
    names(p_values) != "(Intercept)" &
      p_values > p_threshold
  ]
  if (length(p_remove) == 0) {
    message("No predictors exceeded the p-value threshold.")
    return(data_with_selected_variables)
  }
  message("Removing predictors with p-value > ", p_threshold, ": ", 
          paste(p_remove, collapse = ", "))
  if(length(p_remove) > 0) {
    data_cleaned <- data_with_selected_variables[ , 
                                                  - which(names(data_with_selected_variables) %in% 
                                                            p_remove)]
  } else {
    data_cleaned <- data_with_selected_variables
  }
  data_cleaned
}

#' Drop predictors with excessive variance inflation factor (VIF)
#'
#' @param my_model A fitted `lm` object with at least two predictors (VIF is
#'   undefined for single-predictor models).
#' @param data_with_selected_variables Data frame containing `CODE`, the
#'   response variable, and every predictor currently in `my_model`.
#' @param vif_threshold Numeric VIF threshold. Predictors with VIF values
#'   strictly greater than `vif_threshold` are removed.
#'
#' @details
#' Variance inflation factors are calculated using `car::vif()`. The
#' function requires a fitted linear model with at least two predictors.
#'
#' @return A data frame with predictors exceeding the specified VIF
#' threshold removed.
vif_function <- function(my_model, data_with_selected_variables, vif_threshold) {
  if (!inherits(my_model, "lm")) {
    stop("`my_model` must be a fitted lm object.", call. = FALSE)
  }
  if (!is.data.frame(data_with_selected_variables)) {
    stop(
      "`data_with_selected_variables` must be a data frame.",
      call. = FALSE
    )
  }
  if (!is.numeric(vif_threshold) ||
      length(vif_threshold) != 1 ||
      vif_threshold <= 0) {
    stop(
      "`vif_threshold` must be a single positive numeric value.",
      call. = FALSE
    )
  }
  n_predictors <- length(coef(my_model)) - 1
  if (n_predictors < 2) {
    stop(
      "VIF requires a model with at least two predictors.",
      call. = FALSE
    )
  }
  vif_variable <- vif(my_model)
  message("VIF values:\n",
          paste(names(vif_variable), round(vif_variable, 2), collapse = "\n")
  )  
  influential <- names(vif_variable)[vif_variable > vif_threshold]
  if (length(influential) == 0) {
    message("No predictors exceeded the VIF threshold.")
    return(data_with_selected_variables)
  }
  message("Removing predictors with VIF > ",
          vif_threshold, ": ", paste(influential, collapse = ", "))
  if(!identical(influential, character(0))) {
    data_cleaned <- data_with_selected_variables[ , 
                                                  -which(names(data_with_selected_variables) %in% 
                                                           influential)]
  } else {
    data_cleaned <- data_with_selected_variables
  }
  return(data_cleaned)
}

#' Leave-one-out cross-validation (LOOCV)
#'
#' Refits the model once per observation, leaving that observation out of
#' training and predicting it, then records that fold's in-sample R^2 /
#' adjusted R^2 alongside the held-out prediction. Assumes every row is a
#' unique site/`CODE`.
#' 
#' @details
#' Leave-one-out cross-validation is deterministic and does not involve any
#' random sampling. The returned `predicted_loocv_r2` and
#' `predicted_loocv_r2_adj` correspond to the training model fitted after
#' excluding the held-out observation; they are not out-of-sample R² values.
#'
#' @param data_final_selected_variables Data frame with `CODE`, the response
#'   variable, and the final set of selected predictor columns.
#' @param response_variable Character string naming the response column.
#'
#' @return `data_final_selected_variables` with three additional columns:
#'   `predicted_loocv` (the held-out prediction for that row), and
#'   `predicted_loocv_r2` / `predicted_loocv_r2_adj` (the *training* model's
#'   own R^2 / adjusted R^2 for the fold that excluded this row -- **not** an
#'   out-of-sample R^2 for the held-out point itself).
loop_loocv <- function(data_final_selected_variables, response_variable) {
  if (!is.data.frame(data_final_selected_variables)) {
    stop(
      "`data_final_selected_variables` must be a data frame.",
      call. = FALSE
    )
  }
  if (!response_variable %in% names(data_final_selected_variables)) {
    stop(
      "`response_variable` must be a column in `data_final_selected_variables`.",
      call. = FALSE
    )
  }
  if (nrow(data_final_selected_variables) < 2) {
    stop(
      "LOOCV requires at least two observations.",
      call. = FALSE
    )
  }
  data_final_selected_variables <- data_final_selected_variables %>% 
    mutate(predicted_loocv = NA_real_, 
           predicted_loocv_r2 = NA_real_,
           predicted_loocv_r2_adj = NA_real_)
  for(i in seq_len(nrow(data_final_selected_variables))) {
    data_train <- data_final_selected_variables[-i, ]
    my_model <- create_model(data_train, response_variable)
    data_final_selected_variables$predicted_loocv[i] <- predict(my_model, data.frame(data_final_selected_variables[i, ]))
    data_final_selected_variables$predicted_loocv_r2[i] <- summary(my_model)$r.squared
    data_final_selected_variables$predicted_loocv_r2_adj[i] <- summary(my_model)$adj.r.squared
  }
  data_final_selected_variables
}

#' K-fold cross-validation
#'
#' @param data_final_selected_variables Data frame with `CODE`, the response
#'   variable, and the final set of selected predictor columns.
#' @param response_variable Character string naming the response column.
#' @param k Integer, number of folds (default `10`).
#' @param seed Optional integer random seed used when assigning
#' observations to folds. Supplying the same seed produces identical
#' fold assignments and cross-validation results across runs.
#'
#' @details
#' Observations with missing values are removed before cross-validation.
#' When `seed` is supplied, fold allocation is reproducible.
#'
#' @return A data frame (row-bound across folds) with a `predicted_10fold`
#'   column (out-of-fold prediction), `predicted_10fold_r2` /
#'   `predicted_10fold_r2_adj` (that fold's training R^2/adjusted R^2), and a
#'   `fold` index column.
#'
loop_kfold <- function(data_final_selected_variables, response_variable,
                       k = 10, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.numeric(k) || length(k) != 1 || k < 2) {
    stop(
      "`k` must be a single integer greater than or equal to 2.",
      call. = FALSE
    )
  }
  if (k > nrow(data_final_selected_variables)) {
    stop(
      "`k` cannot be greater than the number of observations.",
      call. = FALSE
    )
  }
  data_final_selected_variables <- na.omit(data_final_selected_variables)
  data_final_selected_variables <- data_final_selected_variables[sample(seq_len(nrow(data_final_selected_variables))), ]
  folds <- cut(seq(1, nrow(data_final_selected_variables)), breaks = k, labels = FALSE)
  new_data <- data.frame()
  for(i in seq_len(k)) {
    training_indexes <- which(folds == i, arr.ind = TRUE)
    data_train <- data_final_selected_variables[-training_indexes, ]
    my_model <- create_model(data_train, response_variable)
    cal_new_data <- data.frame(data_final_selected_variables[training_indexes, ])
    cal_new_data$predicted_10fold <- predict(my_model, cal_new_data)
    cal_new_data$predicted_10fold_r2 <- summary(my_model)$r.squared
    cal_new_data$predicted_10fold_r2_adj <- summary(my_model)$adj.r.squared
    cal_new_data$fold <- i
    new_data <- bind_rows(new_data, cal_new_data)
  }
  new_data
}



################################################################################
# Functions for spatial data manipulation
################################################################################

#' Convert a lat/long data frame to an `sf` object (no reprojection)
#'
#' @param file_df Data frame containing longitude/latitude columns plus any
#'   other attributes.
#' @param wgs Character string giving the coordinate reference system of the
#'   lat/long data, e.g. `"+proj=longlat +datum=WGS84 +no_defs"`.
#' @param long Character string, name of the longitude column.
#' @param lat Character string, name of the latitude column.
#'
#' @return An `sf` object in the `wgs` coordinate system.
convert_sf <- function(file_df, wgs, long, lat) {
  missing_cols <- setdiff(c(long, lat), names(file_df))
  if (length(missing_cols) > 0) {
    stop(
      "convert_sf(): column(s) not found in `file_df`: ",
      paste(missing_cols, collapse = ", "), ".\n",
      "Available columns: ", paste(names(file_df), collapse = ", "),
      call. = FALSE
    )
  }
  n_na <- sum(is.na(file_df[[long]]) | is.na(file_df[[lat]]))
  if (n_na > 0) {
    stop(
      "convert_sf(): ", n_na, " row(s) have NA in `", long, "` or `", lat,
      "`. st_as_sf() cannot build geometry for these - filter them out or ",
      "resolve the missing coordinates before converting.",
      call. = FALSE
    )
  }
  file_sf <- st_as_sf(file_df, coords = c(long, lat), crs = st_crs(wgs))
  file_sf
}

#' Convert a lat/long data frame to a projected `sf` object
#'
#' Convert a file with latitude and longitude columns to an `sf` object and
#' transform it into a projected coordinate system (usually the one used for
#' distance/area calculations, e.g. a UTM zone).
#'
#' @param file_df Data frame containing (at minimum) longitude and latitude
#'   columns and any other attribute columns.
#' @param wgs Character string giving the source (geographic) coordinate
#'   reference system, e.g. `"+proj=longlat +datum=WGS84 +no_defs"`.
#' @param UTM_proj Character string giving the target projected coordinate
#'   reference system, e.g. `"+proj=utm +zone=43 +datum=WGS84"`.
#' @param long Character string, name of the longitude column in `file_df`.
#' @param lat Character string, name of the latitude column in `file_df`.
#'
#' @return An `sf` object with geometry in the `UTM_proj` coordinate system.
#'
convert_sf_proj <- function(file_df, wgs, UTM_proj, long, lat) {
  file_sf_proj <- convert_sf(file_df, wgs, long, lat)
  file_sf_proj <- st_transform(file_sf_proj, crs = UTM_proj)
  file_sf_proj
}

#' Read and reproject a vector file
#'
#' Reads vector data (shapefile, GeoPackage, GeoJSON, etc.) from disk and
#' transforms it to a specified coordinate reference system.
#'
#' @param file_sf Path to a vector file readable by `sf::st_read()` (e.g. a
#'   `.shp` path), assumed to already have a geographic coordinate system.
#' @param proj Character string (or CRS object) giving the target coordinate
#'   reference system.
#'
#' @return An `sf` object reprojected to `proj`.
#'
#' @section Reproducibility:
#' Depends on the file at `file_sf` existing at that exact relative/absolute
#' path at run time; no existence check or informative error is raised if it's
#' missing, so failures downstream may be harder to trace back to a bad path.
proj_sf <- function(file_sf, proj) {
  if (!file.exists(file_sf)) {
    stop(
      "proj_sf(): file not found at '", file_sf, "'.\n",
      "Check the path (e.g. via here::here()) and, for shapefiles, that the ",
      "matching .shx/.dbf/.prj sidecar files are present in the same folder.",
      call. = FALSE
    )
  }
  file_sf_proj <- st_read(file_sf, quiet = TRUE)
  file_sf_proj <- st_transform(file_sf_proj, crs = proj)
  file_sf_proj
}

#' Build a named vector of buffer distances
#'
#' Generates a named numeric vector suitable for use with `buffer_points()`,
#' where names encode the buffer radius (e.g. `"rail_buffer_500m"`).
#'
#' @param name Character string, prefix to use for each buffer's name (e.g.
#'   `"rail_buffer_"`).
#' @param buffer_len Numeric vector of buffer radii in metres, e.g.
#'   `c(500, 1000, 5000)`.
#'
#' @return A named numeric vector `buffering`, where `names(buffering)` are of
#'   the form `"<name><radius>m"` and the values are the radii themselves.
lur_buffer_maker <- function(name, buffer_len) {
  # Create a list of buffer lengths
  buff_labs <- paste0(name, buffer_len, "m")
  buffering <- buffer_len
  names(buffering) <- buff_labs
  buffering
}

#' Create buffers of multiple radii around a set of points
#'
#' The buffer width is the radius measured outward from each point's centre.
#'
#' @param buffer Named numeric vector of buffer radii (metres), typically
#'   produced by `lur_buffer_maker()`.
#' @param file_sf_proj A projected `sf` (multi)point object to buffer around.
#'   Must already be in a projected (metric) CRS for the radii to be
#'   meaningful in metres.
#' @param use_names Logical, passed through to `mapply(USE.NAMES = ...)`.
#'   Default `TRUE`, so the returned list is named using `names(buffer)`
#'   (e.g. `"rail_buffer_500m"`). Set `FALSE` to get an unnamed list instead.
#' @param simplify Logical, passed through to `mapply(SIMPLIFY = ...)`.
#'   Default `FALSE`, which keeps the return type a `list` of `sf` objects
#'   regardless of how many buffer radii are supplied -- see Reproducibility
#'   section below before changing this.
#'
#' @return By default (`simplify = FALSE`), a named list of `sf` polygon
#'   objects, one per entry in `buffer`, with names inherited from
#'   `names(buffer)` (unless `use_names = FALSE`).
buffer_points <- function(buffer, file_sf_proj,
                          use_names = TRUE,
                          simplify = FALSE) {
  if (!is.numeric(buffer)) {
    stop(
      "`buffer` must be a numeric vector of radii in metres.",
      call. = FALSE
    )
  }
  if (any(buffer < 0, na.rm = TRUE)) {
    stop(
      "`buffer` cannot contain negative distances.",
      call. = FALSE
    )
  }
  if (is.null(names(buffer)) || any(names(buffer) == "")) {
    stop(
      "buffer_points(): `buffer` must be a fully named numeric vector ",
      "(e.g. the output of lur_buffer_maker()), so results can be looked up ",
      "by name downstream (e.g. buffers$rail_buffer_1000m).",
      call. = FALSE
    )
  }
  if (sf::st_is_longlat(file_sf_proj)) {
    stop(
      "buffer_points(): `file_sf_proj` is in a geographic (lat/long) CRS. ",
      "Buffer radii are assumed to be in metres, which only makes sense for ",
      "a projected CRS - reproject first with convert_sf_proj() or ",
      "sf::st_transform().",
      call. = FALSE
    )
  }
  buffers <- mapply(FUN = sf::st_buffer,
                    dist = buffer,
                    MoreArgs = list(x = file_sf_proj),
                    SIMPLIFY = simplify, 
                    USE.NAMES = use_names)
  buffers
}

#' Extract inverse-distance variables (airport, industries)
#'
#' Calculates the Euclidean distance from each point in `file_sf_proj` to the
#' nearest airport runway point and the nearest industry point, and stores
#' each as an inverse distance (`1 / distance`).
#' @details
#' If a monitoring site coincides exactly with an airport or industry
#' location (distance = 0 m), the distance is replaced with 1 m before
#' calculating the inverse distance. This avoids infinite values while
#' preserving the interpretation of very close proximity.
#'
#' @param file_sf_proj A projected `sf` (multi)point object, one row per
#'   monitoring site / prediction location.
#' @param airport A projected `sf` point object representing runway
#'   location(s). Passed to `st_distance()` against every row of
#'   `file_sf_proj`.
#' @param industries A projected `sf` point object representing industry
#'   locations.
#'   
#' @return `file_sf_proj` with two additional numeric columns,
#'   `inverse_distance_industries` and `inverse_distance_airport`.
extract_dist_variable <- function(file_sf_proj, 
                                  airport = NULL, 
                                  industries = NULL) {
  if (!inherits(file_sf_proj, "sf")) {
    stop("`file_sf_proj` must be an sf object.", call. = FALSE)
  }
  if (!inherits(airport, "sf")) {
    stop("`airport` must be an sf object.", call. = FALSE)
  }
  if (!inherits(industries, "sf")) {
    stop("`industries` must be an sf object.", call. = FALSE)
  }
  if (sf::st_is_longlat(file_sf_proj)) {
    stop(
      "`file_sf_proj` must use a projected CRS because distances are measured in metres.",
      call. = FALSE
    )
  }
  if (sf::st_crs(file_sf_proj) != sf::st_crs(industries)) {
    stop(
      "`industries` must use the same crs as `file_sf_proj`.",
      call. = FALSE
    )
  }
  if (sf::st_crs(file_sf_proj) != sf::st_crs(airport)) {
    stop(
      "`airport` must use the same crs as `file_sf_proj`.",
      call. = FALSE
    )
  }
  airport_dist <- as.numeric(st_distance(file_sf_proj, airport))
  airport_dist <- pmax(airport_dist, 1)
  file_sf_proj <- file_sf_proj %>%
    mutate(
      inverse_distance_industries = NA_real_,
      inverse_distance_airport = 1 / airport_dist
    )
  for(i in seq_len(nrow(file_sf_proj))) {
    min_dist <- which.min(st_distance(file_sf_proj[i, ], industries))
    file_sf_proj[i, "inverse_distance_industries"] <- 
      st_distance(file_sf_proj[i, ], industries[min_dist, ])
  }
  file_sf_proj$inverse_distance_industries <-
    pmax(file_sf_proj$inverse_distance_industries, 1)
  file_sf_proj <- file_sf_proj %>% 
    mutate(inverse_distance_industries =  
             1 / inverse_distance_industries)
  file_sf_proj
}

#' Extract elevation and AOD raster values at each point
#' 
#' @details
#' Elevation values are square-root transformed after extraction.
#' Negative elevation values will produce `NaN` according to the
#' behaviour of `sqrt()`.
#'
#' @param file_sf_proj A projected `sf` (multi)point object.
#' @param dem A `raster` object: digital elevation model.
#' @param aod A `raster` object: aerosol optical depth.
#' @param wgs Character string, geographic CRS matching the rasters
#'   (`dem`/`aod` are assumed to be in this CRS for extraction).
#' @param UTM_proj Character string, the projected CRS to return the object in.
#'
#' @return `file_sf_proj`, reprojected to `wgs` for extraction and back to
#'   `UTM_proj`, with two additional columns: `elevation` (square-rooted DEM
#'   value) and `aod`.
extract_raster <- function(file_sf_proj, dem, aod, wgs, UTM_proj) {
  dem_crs <- st_crs(dem)
  if (!identical(st_crs(wgs), dem_crs)) {
    warning(
      "`dem` CRS does not match the supplied `wgs` CRS.",
      call. = FALSE
    )
  }
  if (!inherits(file_sf_proj, "sf")) {
    stop("`file_sf_proj` must be an sf object.", call. = FALSE)
  }
  if (!inherits(dem, "RasterLayer")) {
    stop("`dem` must be a RasterLayer.", call. = FALSE)
  }
  if (!inherits(aod, "RasterLayer")) {
    stop("`aod` must be a RasterLayer.", call. = FALSE)
  }
  file_sf <- st_transform(file_sf_proj, crs = wgs)
  elev <- raster::extract(dem, file_sf)
  aod_vals <- raster::extract(aod, file_sf)
  if (any(elev < 0, na.rm = TRUE)) {
    warning(
      "Negative elevation values detected. ",
      "Square root will produce NaN values.",
      call. = FALSE
    )
  }
  file_sf <- file_sf %>% 
    mutate(elevation = sqrt(elev),
           aod = aod_vals)
  file_sf_proj <- st_transform(file_sf, crs = UTM_proj)
  file_sf_proj
}

#' Extract railway length within buffers of each site
#'
#' @param buffering_railway Named numeric vector of buffer radii (metres),
#'   typically from `lur_buffer_maker()`.
#' @param file_sf_proj A projected `sf` (multi)point object (monitoring
#'   sites).
#' @param railway A projected `sf` line object of railway geometry.
#' 
#' @details
#' Railway length is calculated as the total length of railway geometry
#' intersecting each circular buffer. Buffers with no intersecting railway
#' receive a value of 0.
#'
#' @return A data frame with one row per `CODE` and one column per buffer
#'   radius, named `rail_buffer_<radius>`, giving the total length (metres) of
#'   railway intersecting that buffer.
extract_railway <- function(buffering_railway, file_sf_proj, railway) {
  if (!inherits(file_sf_proj, "sf")) {
    stop("`file_sf_proj` must be an sf object.", call. = FALSE)
  }
  if (!inherits(railway, "sf")) {
    stop("`railway` must be an sf object.", call. = FALSE)
  }
  if (sf::st_is_longlat(file_sf_proj)) {
    stop(
      "`file_sf_proj` must use a projected CRS because railway lengths are measured in metres.",
      call. = FALSE
    )
  }
  if (sf::st_crs(file_sf_proj) != sf::st_crs(railway)) {
    stop(
      "`railway` and `file_sf_proj` must have the same CRS.",
      call. = FALSE
    )
  }
  buffers <- buffer_points(buffering_railway, file_sf_proj)
  buffers_railway <- mapply(FUN = sf::st_intersection,
                            x = buffers,
                            MoreArgs = list(y = railway),
                            SIMPLIFY = FALSE, 
                            USE.NAMES = TRUE)
  intersection_df <- dplyr::bind_rows(
    lapply(seq_along(buffers_railway), function(i) {
      x <- buffers_railway[[i]]
      x$buffer_m <- names(buffers_railway)[i]
      x
    })
  )
  ## Extract buffer radius from names such as "rail_buffer_500m"
  df_railway <- intersection_df %>%
    mutate(
      buffer_m_1 = sub("rail_buffer_(\\d+)m", "\\1", buffer_m),
      len = as.numeric(sf::st_length(geometry))
    ) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(buffer_m_1, CODE, len) %>%
    group_by(buffer_m_1, CODE) %>%
    summarise(len = sum(len), .groups = "drop") %>%
    pivot_wider(
      names_from = buffer_m_1,
      values_from = len,
      values_fill = 0
    )
  names(df_railway)[-1] <- paste0("rail_buffer_", names(df_railway)[-1])
  df_railway
}

#' Stepwise, direction-of-effect-supervised linear regression
#'
#' Iteratively adds one variable at a time to a running LUR model, at each
#' step keeping only candidate models where *every* included predictor
#' preserves its expected direction of effect, and stops once the best
#' available improvement in adjusted R^2 falls below `change_val`.
#'
#' @param col_interest Character vector of all candidate predictor column
#'   names to consider adding.
#' @param original_para Character string, the name of the single predictor
#'   with the highest univariate adjusted R^2 (and correct direction of
#'   effect) - the model's starting point.
#' @param original_r2 Numeric, the univariate adjusted R^2 of
#'   `original_para`.
#' @param response_variable Character string naming the response column.
#' @param data_with_variables Data frame with `CODE`, the response variable,
#'   and every column in `col_interest`.
#' @param direction_of_effect_table Data frame of expected sign of effect per
#'   base parameter (see `add_sign_check()`).
#' @param change_val Numeric threshold: stop adding variables once the best
#'   achievable increase in adjusted R^2 is smaller than this.
#'
#' @return A list of five elements (unnamed, unpacked via the custom `:=`
#'   operator elsewhere in this file):
#'   \enumerate{
#'     \item `var_list` - character vector of variables in the final model.
#'     \item `change_r` - the adjusted-R^2 improvement at the step where the
#'       loop stopped.
#'     \item `highest_r2` - the best adjusted R^2 achieved.
#'     \item `data_wo_slop_a` - data frame of every candidate model tried at
#'       every step whose predictors all preserved direction of effect.
#'     \item `data_r2` - data frame of only the single best model kept at
#'       each step.
#'   }
run_lur_model <- function(col_interest, original_para, 
                          original_r2, response_variable, 
                          data_with_variables, 
                          direction_of_effect_table, 
                          change_val) {
  ## Make an empty dataframe and try to extract results 
  data_with <- data.frame()
  ## keep track of all the variables with the right direction of effect
  data_wo_slop_a <- data.frame()
  name_col <- as.vector(col_interest)
  ## track the r2 change 
  data_r2 <- data.frame()
  ## Start a variable list 
  var_list <- c(original_para)
  change_r <- 100 ## positive number or number greater than 0.01 for the loop to run at least once
  thres <- change_val ## Define the threshold to come out of the loop
  ## Keep track of the original R2 value
  highest_r2 <- original_r2
  ## Use a while loop to break when change in R2 is less than 0.01
  while(change_r >= thres) {
    ## Loop for the model building and change data_with_variables to selected columns only
    for(i in name_col) {
      ## Check if the parameter is already in the var list, if there then ignore it use others
      ## If available ignore it use others
      if(i %in% var_list) {
        data_with <- data_with
      } else {
        others <- paste(paste(var_list, collapse = " + "), i, sep = " + ")
        ## Generate a formula each time
        equ <- as.formula(paste(response_variable, others, sep = " ~ "))
        ## Apply multiple linear regression
        lm_step <- lm(equ, data = data_with_variables)
        vr <- round((summary(lm_step))$coefficients[, 4], 4)
        if(all(!is.nan(vr)) & all(!is.na(vr))) {
          ## Apply extract_model_data function to extract slope, std error, t value, significance etc and also check the slope sign
          data_corr <- extract_model_data(lm_step, sig_star, direction_of_effect_table)
          ## Extract R2, AIC, RMSE and also track the equation
          data_corr <- data_corr %>% 
            mutate(r2 = round(summary(lm_step)$adj.r.squared, 4), aic = round(AIC(lm_step), 4), 
                   rmse = summary(lm_step)$sigma, eqtn = others)
          data_with <- rbind(data_with, data_corr)
        } else {
          data_with <- data_with
        }
      }
    }
    ## Now grouped by equation remove where slope change and also select the highest R2 change?
    data_wo_slope <- data_with %>%
      group_by(eqtn) %>%
      filter(!any(obj == FALSE)) %>% 
      arrange(desc(r2))
    ## Check if the removed table has any variables or not
    if (dim(data_wo_slope)[1] == 0) {
      data_with <- data_with
    } else {
      data_with <- data.frame()
    }
    ## Check the highest gained R2 
    highest_r2_l <- (data_wo_slope %>% 
                       arrange(desc(r2)))$r2[1]
    ## Check the parameters which gave this high R2 
    highest_para <- data_wo_slope %>% 
      filter(r2 == highest_r2_l)
    change_r <- highest_r2_l - highest_r2 
    if(change_r < 0.01) {
      ## Adding a break statement to move out of loop right away if the condition is satisfied
      break
    } else {
      ## Or else keep tracking the data
      data_r2 <- rbind(data_r2, highest_para)
      data_wo_slop_a <- rbind(data_wo_slop_a, data_wo_slope)
      ## Check if the R2 changed then replace the highest r2 with the new one otherwise keep the same
      if(highest_r2_l > highest_r2) {
        highest_r2 <- highest_r2_l
      } else {
        highest_r2 <- highest_r2
      }
      ## Now in the list of variables check to add in the variable list or to remove it
      for(i in highest_para$value) {
        if(i %in% var_list) {
          ## If added dont add to the list again
          var_list <- var_list
        } else {
          ## Append variable list each time if not added 
          var_list <- append(var_list, i) 
        }
      }
    }
  }
  return(list(var_list, change_r, highest_r2, 
              data_wo_slop_a, data_r2))
}

#' Reshape and group land-use/land-cover (LULC) classes
#'
#' Reshapes a wide LULC-by-buffer data frame into long form, extracts the
#' ESA WorldCover class code and buffer radius from each column name, groups
#' related classes together (20/30/90 -> shrubland), renames classes to
#' human-readable labels, and pivots back to wide form.
#' 
#' @details
#' ESA WorldCover classes 20 (Shrubland), 30 (Grassland) and
#' 90 (Herbaceous wetland) are combined into a single "shrubland"
#' category before aggregation.
#'
#' @param df_lulc Data frame with a `CODE` column and one column per
#'   class/buffer combination (e.g. `"10_buffer_500m"`), as derived from
#'   Earth Engine LULC extraction.
#'
#' @return A wide data frame with `CODE` and one column per
#'   `<land_use_label>_buffer_<radius>` combination, values summed within each
#'   group.
derive_lulc_as_df <- function(df_lulc) {
  lulc_lookup <- tibble::tribble(
    ~name, ~buff,
    "10", "tree_cover",
    "20", "shrubland",
    "40", "cropland",
    "50", "builtup",
    "60", "bare_land",
    "70", "snow_ice",
    "80", "per_water_bod",
    "95", "mang",
    "100", "moss_lichen"
  )
  if (!is.data.frame(df_lulc)) {
    stop("`df_lulc` must be a data frame.", call. = FALSE)
  }
  if (!"CODE" %in% names(df_lulc)) {
    stop("`df_lulc` must contain a `CODE` column.", call. = FALSE)
  }
  df <- df_lulc %>%
    pivot_longer(-CODE, names_to = "Parameter", values_to = "Value") %>% 
    mutate(name = str_extract(Parameter, "(\\d)+(?=_buffer)"),
           name_1 = str_match(Parameter, "buffer_\\s*(.*?)\\s*m")[, 2]) %>% 
    mutate(
      name = case_when(
        name %in% c("20", "30", "90") ~ "20",
        TRUE ~ name)
    ) %>%
    dplyr::select(everything(), - Parameter) %>%
    group_by(CODE, name, name_1) %>%
    summarize(Value = sum(Value, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(lulc_lookup, by = "name") %>% 
    mutate(Parameter = paste0(buff, "_buffer_", name_1)) %>% 
    ungroup() %>% 
    dplyr::select(CODE, Parameter, Value) %>% 
    pivot_wider(names_from = Parameter, 
                values_from = Value,
                values_fill = 0)
}

#' Extract LULC buffer areas for one site (ESA WorldCover, 10 m)
#'
#' For a single site (row `i` of `file_sf`), masks and crops the LULC raster
#' to each buffer radius in `buffer`, then computes the area of each land-use
#' class present within that buffer.
#'
#' @param i Integer row index of the site (within `file_sf`) to process.
#' @param buffer Named numeric vector of buffer radii, from
#'   `lur_buffer_maker()`.
#' @param file_sf A projected `sf` (multi)point object with a `CODE` column.
#' @param lulc A `raster` object: ESA WorldCover (or similar) land-use/land
#'   cover layer, 10 m resolution.
#' @param data_lulc Data frame to append this site's results onto (typically
#'   an empty `data.frame()` on the first call, then the accumulating result
#'   passed back in on subsequent calls in a calling loop).
#'
#' @return `data_lulc` with additional rows for site `i`: columns `CODE`,
#'   `Area`, `variable` (class code + buffer name), `Group.1` (raw class
#'   code).
extract_lulc_parameters <- function(i, buffer, file_sf, lulc, data_lulc) {
  if (!inherits(file_sf, "sf")) {
    stop("`file_sf` must be an sf object.", call. = FALSE)
  }
  if (!inherits(lulc, "Raster")) {
    stop("`lulc` must be a raster object.", call. = FALSE)
  }
  if (!"CODE" %in% names(file_sf)) {
    stop("`file_sf` must contain a `CODE` column.", call. = FALSE)
  }
  ## Create a temp folder to store intermediate data
  ## Code to reduce burden of c drive, make a temp folder in D drive and unlink the folder each time, saves space
  OWNER <- "del"
  dir.create(file.path("D:/", OWNER), showWarnings = FALSE)
  raster::rasterOptions(tmpdir = file.path("D:/", OWNER))
  ## Generate the buffers for each point
  ## width is the radius from the center
  ## Generate buffers of the above specified length for each of the multipoint observation
  buffers <- buffer_points(buffer, file_sf[i, ])
  lulc_clip_buffer <- mapply(FUN = raster::mask,
                             mask = buffers,
                             MoreArgs = list(x = lulc),
                             SIMPLIFY = FALSE, 
                             USE.NAMES = TRUE)
  for(buffer_name in names(lulc_clip_buffer)) {
    ## Crop the masked raster using each of the buffers and the buffer masked lulc one by one
    cropped_lulc <- raster::crop(lulc_clip_buffer[[buffer_name]], buffers[[buffer_name]])
    ## Calculate area / mean / sum for each land use type and track the buffer and type of land use
    lulc_summary <- as.data.frame(raster::aggregate(raster::getValues(raster::area(cropped_lulc, weights = FALSE)), 
                                               by = list(raster::getValues(cropped_lulc)), FUN = sum)) %>% 
      mutate(CODE = file_sf[i, "CODE"][[1]], variable = paste0(Group.1, "_", buffer_name)) %>% 
      dplyr::select(CODE, "Area" = x, variable, Group.1) %>% 
      filter(Group.1 != 0)
    ## bind for all the buffers and all other the points / CODE / observations
    data_lulc <- bind_rows(data_lulc, lulc_summary)
  }
  data_lulc
}

#' Pivot extracted LULC data to wide form and write to CSV
#'
#' @param data_lulc Data frame produced by
#'   `extract_lulc_parameters()`.
#' @param name File path where the CSV will be written.
#'
#' @return Invisibly returns the data frame written to `name`.
write_lulc_vars <- function(data_lulc, name) {
  if (!is.data.frame(data_lulc)) {
    stop(
      "`data_lulc` must be a data frame.",
      call. = FALSE
    )
  }
  required_cols <- c("Area", "variable", "CODE")
  missing_cols <- setdiff(required_cols, names(data_lulc))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.character(name) || length(name) != 1) {
    stop(
      "`name` must be a single file path.",
      call. = FALSE
    )
  }
  df_lulc <- data_lulc %>% 
    dplyr::select(Area, variable, CODE) %>% 
    tidyr::pivot_wider(names_from = variable, 
                values_from = Area,
                values_fill = 0)
  write.csv(df_lulc, name,
            row.names = FALSE)
  invisible(df_lulc)
}

#' Extract mean NDVI within buffers, for every site, and write to CSV
#' 
#' @details
#' Temporary raster files are written to `tmp_dir`, which defaults to
#' `tempdir()` for portability across operating systems.
#'
#' @param buffer Named numeric vector of buffer radii, from
#'   `lur_buffer_maker()`.
#' @param file_sf A projected `sf` (multi)point object with a `CODE` column.
#' @param ndvi A `raster` object: NDVI layer (e.g. derived from Sentinel-2,
#'   10 m).
#' @param name File path/name to write the resulting CSV to.
#' @param tmp_dir Directory used to store temporary raster files.
#'   Defaults to `tempdir()`.
#'
#' @return Invisibly returns the data frame written to `name`.
extract_ndvi_parameters <- function(buffer, file_sf, ndvi, name,
                                    tmp_dir = tempdir()) {
  if (!inherits(file_sf, "sf")) {
    stop("`file_sf` must be an sf object.", call. = FALSE)
  }
  if (!inherits(ndvi, "Raster")) {
    stop("`ndvi` must be a raster object.", call. = FALSE)
  }
  if (!"CODE" %in% names(file_sf)) {
    stop("`file_sf` must contain a `CODE` column.", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1) {
    stop("`name` must be a single file path.", call. = FALSE)
  }
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  old_options <- raster::rasterOptions(tmpdir = tmp_dir)
  on.exit({
    raster::rasterOptions(tmpdir = old_options$tmpdir)
    unlink(tmp_dir, recursive = TRUE)
  }, add = TRUE)  
  data_ndvi <- data.frame()
  for(i in seq_len(nrow(file_sf))) {
    buffers <- buffer_points(buffer, file_sf[i, ])
    ndvi_buffer <- mapply(FUN = raster::mask,
                          mask = buffers,
                          MoreArgs = list(x = ndvi),
                          SIMPLIFY = FALSE, USE.NAMES = TRUE)
    for(buffer_name in names(ndvi_buffer)) {
      ndvi_values <- raster::getValues(ndvi_buffer[[buffer_name]])
      mean_ndvi <- mean(ndvi_values, na.rm = TRUE)
      table_generated <- data.frame(CODE = file_sf[i, "CODE"][[1]], 
                                    variable = buffer_name, mean_ndvi = mean_ndvi) 
      data_ndvi <- dplyr::bind_rows(data_ndvi, table_generated)
    } 
  }
  df_ndvi <- data_ndvi %>% 
    dplyr::select(mean_ndvi, variable, CODE) %>% 
    pivot_wider(names_from = variable, 
                values_from = mean_ndvi,
                values_fill = NA_real_)
  write.csv(df_ndvi, file = name, row.names = FALSE)
  invisible(df_ndvi)
}

#' Extract summed population within buffers, for every site, and write to CSV
#'
#' @details
#' Temporary raster files are written to `tmp_dir`, which defaults to
#' `tempdir()` for portability across operating systems.
#' 
#' @param buffer Named numeric vector of buffer radii, from
#'   `lur_buffer_maker()`. **Not currently used** - see Reproducibility
#'   note below.
#' @param file_sf A projected `sf` (multi)point object with a `CODE` column.
#' @param pop A `raster` object: population layer.
#' @param name File path/name to write the resulting CSV to.
#' @param tmp_dir Directory used to store temporary raster files.
#'   Defaults to `tempdir()`.
#'
#' @return Invisibly, nothing meaningful; called for the side effect of
#'   writing `name` to disk (summed population per site per buffer radius,
#'   wide form).
extract_pop_vars <- function(buffer, file_sf, pop, name,
                             tmp_dir = tempdir()) {
  if (!inherits(file_sf, "sf")) {
    stop("`file_sf` must be an sf object.", call. = FALSE)
  }
  if (!inherits(pop, "Raster")) {
    stop("`pop` must be a raster object.", call. = FALSE)
  }
  if (!"CODE" %in% names(file_sf)) {
    stop("`file_sf` must contain a `CODE` column.", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1) {
    stop("`name` must be a single file path.", call. = FALSE)
  }
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  old_options <- raster::rasterOptions(tmpdir = tmp_dir)
  on.exit({
    raster::rasterOptions(tmpdir = old_options$tmpdir)
    unlink(tmp_dir, recursive = TRUE)
  }, add = TRUE)  
  data_population <- data.frame()
  for(i in seq_len(nrow(file_sf))) {
    buffers <- buffer_points(buffering, file_sf[i, ])
    pop_buffer <- mapply(FUN = raster::mask,
                         mask = buffers,
                         MoreArgs = list(x = pop),
                         SIMPLIFY = FALSE, USE.NAMES = TRUE)
    for(buffer_name in names(pop_buffer)) {
      population_values <- raster::getValues(pop_buffer[[buffer_name]])
      summed_population <- sum(
        population_values,
        na.rm = TRUE
      )
      table_generated <- data.frame(
        CODE = file_sf[i, "CODE"][[1]],
        variable = buffer_name,
        sum_pop = summed_population
      )
      data_population <- dplyr::bind_rows(
        data_population,
        table_generated
      )
    } 
  }
  df_population <- data_population %>% 
    dplyr::select(sum_pop, variable, CODE) %>% 
    tidyr::pivot_wider(names_from = variable, 
                       values_from = sum_pop)
  write.csv(df_population, name,
            row.names = FALSE)
  invisible(df_population)
}

################################################################################
# Functions for road data manipulation
################################################################################
# road <- road %>% 
#   filter(fclass %in% c("service", "tertiary", "secondary", "residential", "unclassified", "motorway",
#                        "primary", "trunk", "secondary_link", "living_street", "trunk_link", "track",
#                        "primary_link", "motorway_link", "tertiary_link")) 

#' Intersect roads (OSM) with site buffers
#'
#' @param buffer Named numeric vector of buffer radii, from
#'   `lur_buffer_maker()`. Note the buffer names are expected to look like
#'   `"roads_buffer_<radius>m"` for the regex in this function to correctly
#'   extract the radius afterwards.
#' @param file_sf A projected `sf` (multi)point object with a `CODE` column.
#' @param road A projected `sf` line object of OSM road geometry, expected to
#'   have an `fclass` column (OSM road-class label).
#'
#' @return An `sf` data frame of road segments clipped to each buffer, with
#'   `buffer_m` (raw buffer name), `buffer_m_1` (extracted numeric radius as
#'   character), and `len` (segment length in metres).
extract_road_vars <- function(buffer, file_sf, road) {
  if (!inherits(file_sf, "sf")) {
    stop("`file_sf` must be an sf object.", call. = FALSE)
  }
  if (!inherits(road, "sf")) {
    stop("`road` must be an sf object.", call. = FALSE)
  }
  if (!"CODE" %in% names(file_sf)) {
    stop("`file_sf` must contain a `CODE` column.", call. = FALSE)
  }
  if (!"fclass" %in% names(road)) {
    stop("`road` must contain an `fclass` column.", call. = FALSE)
  }
  if (sf::st_crs(file_sf) != sf::st_crs(road)) {
    stop(
      "`road` and `file_sf` must use the same CRS.",
      call. = FALSE
    )
  }
  buffers <- buffer_points(buffer, file_sf)
  ## Intersect roads with each buffer
  road_buffer <- mapply(
    FUN = sf::st_intersection,
    x = buffers,
    MoreArgs = list(y = road),
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE
  )
  ## Combine while preserving the buffer name
  intersection_df <- dplyr::bind_rows(
    lapply(seq_along(road_buffer), function(i) {
      x <- road_buffer[[i]]
      x$buffer_m <- names(road_buffer)[i]
      x
    })
  )
  ## Extract buffer radius and calculate road lengths
  data_road <- intersection_df %>%
    mutate(
      buffer_m_1 = sub("roads_buffer_(\\d+)m", "\\1", buffer_m),
      len = as.numeric(sf::st_length(geometry))
    )
  data_road
}

#' Sum road length per buffer radius per site
#'
#' @param data_road An `sf`/data frame produced by `extract_road_vars()` (or
#'   a subset of it filtered by `fclass`), with columns `buffer_m_1`, `CODE`,
#'   `len`.
#'
#' @return A wide data frame with `CODE` and one column per buffer radius,
#'   giving total road length (metres) summed within that buffer.
extract_road_filtered <- function(data_road) {
  data_road_filter <- data_road %>% 
    dplyr::select(buffer_m_1, CODE, len) %>% 
    group_by(buffer_m_1, CODE) %>% 
    summarise_all(list(~ sum(., na.rm = TRUE))) %>% 
    pivot_wider(names_from = buffer_m_1, 
                values_from = len,
                values_fill = 0)
  data_road_filter
}

################################################################################
# names of different road types
################################################################################

road_buff <- "roads_buffer_"
res_road_buff <- "roads_res_buffer_"
high_road_buff <- "roads_high_buffer_"
ter_road_buff <- "roads_ter_buffer_"

# Residential
# res <- "residential"
# liv <- "living_street"

# Major roads, primary / highway
# pri <- "primary"
# motor <- "motorway"
# trunk <- "trunk" 

# Arterial, tertiary
# sec <- "secondary"
# ter <- "tertiary"
# mo_li <- "motorway_link"
# pr_li <- "primary_link"
# ter_li <- "tertiary_link"
# tr_li <- "trunk_link"
# ser <- "service"

#' Default OpenStreetMap road classifications
#'
#' A named list defining the OpenStreetMap (`fclass`) values used to classify
#' roads into residential, highway and arterial categories when constructing
#' road-length predictor variables.
#'
#' @format A named list with three elements:
#' \describe{
#'   \item{residential}{Character vector of residential road classes.}
#'   \item{highway}{Character vector of major highway classes.}
#'   \item{arterial}{Character vector of arterial road classes.}
#' }
road_classes <- list(
  residential = c(
    "residential",
    "living_street"
  ),
  highway = c(
    "primary",
    "motorway",
    "trunk"
  ),
  arterial = c(
    "secondary",
    "tertiary",
    "service",
    "primary_link",
    "motorway_link",
    "trunk_link",
    "tertiary_link"
  )
)

#' Summarise road length within buffers
#'
#' Calculates the total road length within each buffer radius and monitoring
#' site. Optionally filters to a subset of OpenStreetMap road classes before
#' summarising.
#'
#' @param data A data frame produced by `extract_road_vars()`. Must contain
#'   columns `CODE`, `buffer_m_1`, `len`, and `fclass`.
#' @param classes Optional character vector of OpenStreetMap road classes
#'   (values of `fclass`) to retain before summarising. If `NULL`
#'   (default), all road classes are included.
#' @param prefix Character string used as the prefix for the output column
#'   names.
#'
#' @return A wide data frame with one row per monitoring site (`CODE`) and
#'   one column per buffer radius containing total road length (metres).
summarise_roads <- function(data,
                            classes = NULL,
                            prefix) {
  required_cols <- c("CODE", "buffer_m_1", "len", "fclass")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.null(classes)) {
    data <- dplyr::filter(data, fclass %in% classes)
  }
  out <- extract_road_filtered(data)
  names(out)[-1] <- paste0(prefix, names(out)[-1])
  out
}

#' Create road-length predictor variables
#'
#' Summarises road length within each buffer radius for four groups:
#'
#' * all roads
#' * residential roads
#' * major highways
#' * arterial roads
#'
#' The resulting data frames are merged by monitoring site (`CODE`) to create
#' the final road predictor table used in the land-use regression model.
#'
#' @param data_road An `sf` object returned by `extract_road_vars()`.
#' @param road_classes Named list defining the OpenStreetMap road classes used
#'   for each road category. The default groups roads into
#'   `"residential"`, `"highway"` and `"arterial"` classes.
#'
#' @return A data frame with one row per monitoring site (`CODE`) containing
#'   road-length predictor variables for each buffer radius.
#'
#' @examples
#' road_classes
#'
#' road_predictors <- extract_road_parameters(
#'   data_road,
#'   road_classes = road_classes
#' )
extract_road_parameters <- function(
    data_road,
    road_classes = road_classes
) {
  data_road <- sf::st_drop_geometry(data_road)
  all_roads <-
    summarise_roads(
      data_road,
      prefix = "roads_buffer_"
    )
  residential <-
    summarise_roads(
      data_road,
      road_classes$residential,
      "roads_res_buffer_"
    )
  highway <-
    summarise_roads(
      data_road,
      road_classes$highway,
      "roads_high_buffer_"
    )
  arterial <-
    summarise_roads(
      data_road,
      road_classes$arterial,
      "roads_ter_buffer_"
    )
  reduce(
    list(all_roads, residential, highway, arterial),
    full_join,
    by = "CODE"
  )
}


#' Tune a random forest model over a grid of hyperparameters
#'
#' Fits a sequence of random forest models using combinations of
#' hyperparameters supplied in `hyper_grid`, and evaluates each model using
#' both out-of-bag (OOB) and training-set performance metrics.
#'
#' @param hyper_grid A data frame containing one row per hyperparameter
#'   combination. Must contain the columns:
#'   `num.trees`, `splitrule`, `mtry`, `max.depth`,
#'   `min.node.size`, `replace`, and `sample.fraction`.
#' @param formula A model formula supplied either as a formula object or a
#'   character string (e.g. `"PM2.5 ~ road + ndvi"`).
#' @param train A data frame containing the training data.
#' @param response Character string giving the name of the response variable
#'   in `train`.
#' @param seed Integer random seed passed to `ranger()`.
#'   Defaults to `108`.
#'
#' @details
#' For each row of `hyper_grid`, a random forest model is fitted using
#' `ranger::ranger()`. The following performance statistics are calculated:
#'
#' * Out-of-bag RMSE
#' * Out-of-bag R-squared
#' * Training RMSE
#' * Training R-squared
#'
#' These metrics are appended to `hyper_grid`.
#'
#' @return
#' The input `hyper_grid` with four additional columns:
#' * `OOB_RMSE`
#' * `OOB_R2`
#' * `valid_RMSE`
#' * `valid_R2`
tune_rf <- function(hyper_grid,
                    formula,
                    train,
                    response,
                    seed = 108) {
  required_cols <- c(
    "num.trees",
    "splitrule",
    "mtry",
    "max.depth",
    "min.node.size",
    "replace",
    "sample.fraction"
  )
  missing_cols <- setdiff(required_cols, names(hyper_grid))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in `hyper_grid`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!response %in% names(train)) {
    stop(
      "`response` is not a column in `train`.",
      call. = FALSE
    )
  }
  if (is.character(formula)) {
    formula <- as.formula(formula)
  }
  hyper_grid$OOB_RMSE <- NA_real_
  hyper_grid$OOB_R2 <- NA_real_
  hyper_grid$valid_RMSE <- NA_real_
  hyper_grid$valid_R2 <- NA_real_
  obs <- train[[response]]
  for (i in seq_len(nrow(hyper_grid))) {
    fit <- ranger::ranger(
      formula = formula,
      data = train,
      num.trees = hyper_grid$num.trees[i],
      splitrule = hyper_grid$splitrule[i],
      mtry = hyper_grid$mtry[i],
      max.depth = hyper_grid$max.depth[i],
      min.node.size = hyper_grid$min.node.size[i],
      replace = hyper_grid$replace[i],
      sample.fraction = hyper_grid$sample.fraction[i],
      respect.unordered.factors = "order",
      seed = seed,
      verbose = FALSE
    )
    pred <- predict(fit, train)$predictions
    hyper_grid$OOB_RMSE[i] <- sqrt(fit$prediction.error)
    hyper_grid$OOB_R2[i] <- fit$r.squared
    hyper_grid$valid_RMSE[i] <-
      sqrt(mean((obs - pred)^2))
    hyper_grid$valid_R2[i] <-
      1 - mean((obs - pred)^2) / stats::var(obs)
  }
  hyper_grid
}