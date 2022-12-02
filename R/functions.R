# Important functions file, run this in the beginning before starting the analysis
library(tidyverse)
library(sf)
library(raster)
library(readxl)
library(car)
library(data.table)
library(here)
library(caret)

# Define the projection system 1
wgs <- "+proj=longlat +datum=WGS84 +no_defs"
# Define the projection system 2 which is usually a projected coordinate system to calculate buffer and distances in meters
UTM_proj <- "+proj=utm +zone=43 +datum=WGS84" 

# functions to calculate basic statistics, parameter for all function is a col: column to calculate;
# return the calculated columns
# GSD1 : Goemetric Standard Deviation
GSD1 <- function(col, na.rm = TRUE) {
  GSD_col <- exp(sd(log(col), na.rm = TRUE))
  return(GSD_col)
}
# GM1 : Goemetric Mean
GM1 <- function(col, na.rm = TRUE) {
  GM_col <- exp(mean(log(col), na.rm = TRUE)) 
  return(GM_col)
}
# CV1 : Coefficient of variation
CV1 <- function(col, na.rm = TRUE) {
  cv <- sd(col, na.rm = TRUE) / mean(col, na.rm = TRUE)
  return(cv)
}
# stderr : Standard Error
stderr <- function(col, na.rm = FALSE) {
  if (na.rm) col <- na.omit(col)
  sqrt(var(col) / length(col))
}

# Functions for spatial data manipulation
# Convert a file with latitude and longitude to a `sf` object and transform the projection of that file (usualy the projected coordinate system)
# Parameters in this function are file : file with lat, long column and other attributes; wgs : coordinate system of the lat long format;
# UTM proj: resultant projection; lat: column name of latitude column; long: column name of longitude column
# returns a sf object with the specified UTM_proj projection
convert_sf_proj <- function(file, wgs, UTM_proj, long, lat) {  
  # Convert file to sf
  file <- st_as_sf(file, coords = c(long, lat), crs = st_crs(wgs))
  # Transform the data's projection
  file <- st_transform(file, crs = UTM_proj)
}

# Convert a file with latitude and longitude to a `sf` object 
# Parameters in this function are file : file with lat, long column and other attributes; 
# wgs : coordinate system of the lat long format;
# lat: column name of latitude column; long: column name of longitude column
# returns a sf object with wgs projection
convert_sf <- function(file, wgs, long, lat) {  
  # Convert file to sf
  file <- st_as_sf(file, coords = c(long, lat), crs = st_crs(wgs))
}

# Convert shapefile to a particular coordinate system proj
# Parameters in this function are file : a shapefile. geopackage, geojson etc, a vector data; 
# proj : the coordinate system required;
# returns a sf object with the specified proj projection
sf_proj <- function(file, proj) {
  proj_file <- st_read(file)
  proj_file <- st_transform(proj_file, crs = proj)
}
# Generate a named vector for lur buffers
# parameters are name : name of the variable to create buffer on;
# len_buffer : a vector of buffers in meters to be generated
# returns a named numeric vector with the names as "name_buffer_500m", etc.
lur_buffer_maker <- function(name, len_buffer) {
  # Create a list of buffer engths
  buff_labs <- paste0(name, len_buffer, "m")
  buffering <- len_buffer
  names(buffering) <- buff_labs
  buffering
}


# Generate a named vector for lur buffers
# parameters are file : projected file with latitude longitude or projected sptial points to create buffers on;
# buffer : a vector of buffers in meters to be generated
# returns a named numeric vector with the names as "name_buffer_500m", etc.
# width of the buffer is the radius from the center
# Generate buffers of the above specified length for each of the points 
buffer_points <- function(buffer, file) {
  buffers <- mapply(FUN = sf::st_buffer,
                    dist = buffer,
                    MoreArgs = list(x = file),
                    SIMPLIFY = FALSE, 
                    USE.NAMES = TRUE)
}
buffers <- buffer_points(buffering, file_sf)


# Assign at once multiple returned values from a function to multiple variables
# https://stackoverflow.com/a/1829651
':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) 
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir = frame)
  return(invisible(NULL)) 
}

# Get significance code for regression coefficients
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

# Extract the results (slope, tval, prob etc) from a linear regression model and also check slope sign for each variable
# all_data : is the linear regression model;
# sig_star : the function above to generate the significance star;
# param_t is the parameter list of the desired sign of effect of each of the variable
# returns a data frame which has the all 
vars <- function(all_data, sig_star, param_t) {
  data_ext <- data.frame(
    slope = format((summary(all_data))$coefficients[, 1], scientific = TRUE),
    stde = round((summary(all_data))$coefficients[, 2], 4),
    tval = round((summary(all_data))$coefficients[, 3], 4),
    prob = round((summary(all_data))$coefficients[, 4], 4),
    sig = unlist(lapply((
      summary(all_data)
    )$coefficients[, 4], FUN = sig_star))
  )
  data_ext <- cbind(value = rownames(data_ext), data_ext) %>% 
    filter(value != "(Intercept)")
  data_ext <- add_sign_check(data_ext, param_t)
}

# Function to check sign of the slope for direction of effect
check_s <- function(sign, slope, val) {
  if(!is.na(val)){
    do.call(sign, list(slope, val)) # Change this to 0 
  } else if(is.na(val)){
    TRUE # Change this to NA 
  }
} 

# st_distance calculates Euclidean distance to each point
# Extract values of the roads, airport, industries corresponding for each point
dist_variable_extraction <- function(file_sf, airport_sf = NULL, industries = NULL) {
  file_sf <- file_sf %>%
    mutate(distance_industries = NA,
           distance_airport = as.numeric((1 / st_distance(., airport_sf, by_element = TRUE))))
  for(i in 1:nrow(file_sf)) {
    min_dist <- which.min(st_distance(file_sf[i, ], industries))
    file_sf[i, "distance_industries"] <- st_distance(file_sf[i, ], industries[min_dist, ])
  }
  file_sf <- file_sf %>% 
    mutate(distance_industries =  1 / distance_industries)
  return(file_sf)
}


extract_raster <- function(g_proj_file, dem, aod, wgs, UTM_proj) {
  g_proj_file <- st_transform(g_proj_file, crs = wgs)
  file <- g_proj_file %>% 
    mutate(elevation = sqrt(raster::extract(dem, .)),
           aod = raster::extract(aod, .))
  proj_file <- st_transform(file, crs = UTM_proj)
}

# Function to add sign from the existing LUR sheet and then check for sign and add an obj variable which shows sign preserved or not
add_sign_check <- function(slope_table, param_t) {
  slope_table$param <- sub("_buffer.*", "\\1", slope_table$value)
  slope_table <- left_join(slope_table, param_t, by = "param")
  for(i in 1:nrow(slope_table)) {
    slope_table$obj[i] <- check_s(slope_table$sign[i], 
                                  as.numeric(as.character(slope_table$slope[i])), 
                                  as.numeric(as.character(slope_table$val[i])))
  }
  return(slope_table)
}

create_model <- function(data, resp) {
  new_data <- data %>% 
    dplyr::select(everything(), -!!resp, -"CODE", -contains(c("predicted")))
  others <- names(new_data)
  if(ncol(new_data) == 1){  
    equ <- as.formula(paste(resp, others, sep = " ~ "))
  } else { 
    equ <- as.formula(paste(resp, paste(others, collapse = "+"), sep = " ~ "))
  }
  all_data <- lm(as.formula(equ), data = data)
  return(all_data)
}

# Derive railway parameters
railway_fun <- function(buffering, file_sf, railway) {
  buffers <- buffer_points(buffering, file_sf)
  buffers_r <- mapply(FUN = sf::st_intersection,
                      x = buffers,
                      MoreArgs = list(y = railway),
                      SIMPLIFY = FALSE, 
                      USE.NAMES = TRUE)
  df_1 <- do.call(rbind, lapply(buffers_r, as.data.frame))
  df <- cbind(buffer_m = rownames(df_1), df_1)
  df$buffer_m_1 <- sub(".*buffer_*(.*?) *m.*", "\\1", df$buffer_m)
  df <- st_as_sf(df)
  df$len <- as.numeric(as.character(st_length(df)))
  df_railway <- as.data.frame(df) %>% 
    dplyr::select(buffer_m_1, CODE, len) %>% 
    group_by(buffer_m_1, CODE) %>% 
    summarise_all(list(~ sum(., na.rm = TRUE))) %>% 
    pivot_wider(names_from = buffer_m_1, values_from = len)
  colnames(df_railway)[-1] <- paste0("rail_buffer_", colnames(df_railway)[-1])
  return(df_railway)
}

run_lur_model <- function(col_interest, original_para, original_r2, resp, 
                          df_sum_selected, param_t, change_val) {
  # Make an empty data frame and try to extract results 
  data_with <- data.frame()
  # keep track of all the variables with the right direction of effect
  data_wo_slop_a <- data.frame()
  name_col <- as.vector(col_interest)
  # remember the r2 change 
  data_r2 <- data.frame()
  # Start a variable list 
  var_list <- c(original_para)
  data_exclude <- data.frame()
  change_r <- 100 # positive number or number greater than 0.01 for the loop to run at least once
  thres <- change_val # Define your threshold to come out of the loop
  # Keep track of the original R2 value
  highest_r2 <- original_r2
  # Use a while loop to break when change in R2 is less than 0.01
  while(change_r >= thres) {
    # Loop for the model building and change df_sum_selected to selected columns only
    for(i in name_col) {
      # Check if the parameter is already in the var list, if there then ignore it use others
      # If available ignore it use others
      if(i %in% var_list) {
        data_with <- data_with
        data_exclude <- data_exclude
      } else {
        others <- paste(paste(var_list, collapse = " + "), i, sep = " + ")
        # Generate a formula each time
        equ <- as.formula(paste(resp, others, sep = " ~ "))
        # Apply multiple linear regression
        lm_step <- lm(equ, data = df_sum_selected)
        vr <- round((summary(lm_step))$coefficients[, 4], 4)
        if(all(!is.nan(vr)) & all(!is.na(vr))) {
          # Apply vars function to extract slope, std error, t value, significance etc and also check the slope sign
          data_corr <- vars(lm_step, sig_star, param_t)
          # Extract R2, AIC, RMSE and also track the equation
          data_corr <- data_corr %>% 
            mutate(r2 = round(summary(lm_step)$adj.r.squared, 4), aic = round(AIC(lm_step), 4), 
                   rmse = summary(lm_step)$sigma, eqtn = others)
          data_with <- rbind(data_with, data_corr)
        } else {
          data_ex <- as.data.frame(others)
          data_exclude <- rbind(data_exclude, data_ex)
          data_with <- data_with
        }
      }
    }
    # Now grouped by equation remove where slope change and also select the highest R2 change?
    data_wo_slope <- data_with %>%
      group_by(eqtn) %>%
      filter(!any(obj == FALSE)) %>% 
      arrange(desc(r2))
    # Check if the removed table has any variables or not
    if (dim(data_wo_slope)[1] == 0) {
      data_with <- data_with
    } else {
      data_with <- data.frame()
    }
    # Check the highest gained R2 
    highest_r2_l <- (data_wo_slope %>% 
                       arrange(desc(r2)))$r2[1]
    # Check the parameters which gave this high R2 
    highest_para <- data_wo_slope %>% 
      filter(r2 == highest_r2_l)
    change_r <- highest_r2_l - highest_r2 
    if(change_r < 0.01) {
      # Add a break statement here to come out of loop right away if the condition is satisfied
      break
    } else {
      # Or else keep noting the data
      data_r2 <- rbind(data_r2, highest_para)
      data_wo_slop_a <- rbind(data_wo_slop_a, data_wo_slope)
      # Check if the R2 changed then replace the highest r2 with the new one otherwise keep the same
      if(highest_r2_l > highest_r2) {
        highest_r2 <- highest_r2_l
      } else {
        highest_r2 <- highest_r2
      }
      # Now in the list of variables check to add in the variable list or to remove it
      for(i in highest_para$value) {
        if(i %in% var_list) {
          # If added dont add to the list again
          var_list <- var_list
        } else {
          # Append variable list each time if not added 
          var_list <- append(var_list, i) 
        }
      }
    }
  }
  return(list(var_list, change_r, highest_r2, data_wo_slop_a, data_r2, data_exclude))
}

# Remove variables with insignificant p value
remove_p_value <- function(model, df_sum_selected, no) {
  p_values <- round((summary(model))$coefficients[, 4], 4)
  influential_p_remove <- names(p_values)[(p_values > no)]
  print(c("The columns removed due to p-value greater than 0.1 are ", influential_p_remove))
  if(!identical(influential_p_remove, character(0))) {
    df_sum_selected_f <- df_sum_selected[ , - which(names(df_sum_selected) %in% influential_p_remove)]
  } else {
    df_sum_selected_f <- df_sum_selected
  }
  return(df_sum_selected_f)
}

## vif_of_reg <- 1 / (1 - summary(all_data)$adj.r.squared)
vif_function <- function(model, df_sum_selected, no) {
  vif_variable <- vif(model)
  print(c("The values of vif are ", vif_variable))
  # Find the columns with high VIF and remove them from the database
  influential <- names(vif_variable)[(vif_variable > no)]
  if(!identical(influential, character(0))) {
    df_sum_selected_f <- df_sum_selected[ , -which(names(df_sum_selected) %in% influential)]
  } else {
    df_sum_selected_f <- df_sum_selected
  }
  return(df_sum_selected_f)
}

# Change the LULC extracted data 
derive_lulc_as_df <- function(df_lulc) {
  df <- df_lulc %>%
    pivot_longer(-CODE, names_to = "Parameter", values_to = "Value") %>% 
    mutate(name = str_extract(Parameter, "(\\d)+(?=_buffer)"),
           name_1 = str_match(Parameter, "buffer_\\s*(.*?)\\s*m")[, 2]) %>% 
    mutate(name = ifelse((name == 20 | name == 30 | name == 90), 20, name)) %>%
    dplyr::select(everything(), - Parameter) %>%
    group_by(CODE, name, name_1) %>%
    summarize(Value = sum(Value, na.rm = TRUE)) %>%
    mutate(buff = case_when(
      name == 10 ~ 'tree_cover', name == 20 ~ 'shrubland', name == 40 ~ 'cropland', 
      name == 50 ~ 'builtup', name == 60 ~ 'bare_land', name == 70 ~ 'snow_ice', 
      name == 80 ~ 'per_water_bod', name == 95 ~ 'mang', name == 100 ~ 'moss_lichen', 
      TRUE ~ NA_character_)) %>% 
    mutate(Parameter = paste0(buff, "_buffer_", name_1)) %>% 
    ungroup() %>% 
    dplyr::select(CODE, Parameter, Value) %>% 
    pivot_wider(names_from = Parameter, values_from = Value)
}

# Loop LOOCV
loocv_loop <- function(data, resp) {
  data <- data %>% 
    mutate(predicted_loocv = NA, 
           predicted_loocv_r2 = NA,
           predicted_loocv_r2_adj = NA)
  for(i in 1:nrow(data)) {
    data_train <- data[-i, ]
    mdl <- create_model(data_train, resp)
    data$predicted_loocv[i] <- predict(mdl, data.frame(data[i, ]))
    data$predicted_loocv_r2[i] <- summary(mdl)$r.squared
    data$predicted_loocv_r2_adj[i] <- summary(mdl)$adj.r.squared
  }
  return(data)
}

# Loop 10 fold
foldk_loop <- function(data, resp, k = 10) {
  data <- data[sample(nrow(data)), ]
  folds <- cut(seq(1, nrow(data)), breaks = k, labels = FALSE)
  new_data <- data.frame()
  for(i in 1:k) {
    training_indexes <- which(folds == i, arr.ind = TRUE)
    data_train <- data[-training_indexes, ]
    mdl <- create_model(data_train, resp)
    cal_new_data <- data.frame(data[training_indexes, ])
    cal_new_data$predicted_10fold <- predict(mdl, cal_new_data)
    cal_new_data$predicted_10fold_r2 <- summary(mdl)$r.squared
    cal_new_data$predicted_10fold_r2_adj <- summary(mdl)$adj.r.squared
    new_data <- rbind(new_data, cal_new_data)
  }
  return(new_data)
}