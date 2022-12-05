# functions file, run this in the beginning before starting the analysis
library(tidyverse)
library(sf)
library(raster)
library(readxl)
library(car)
library(data.table)
library(here)
library(caret)
library(tmap)
library(gt)





################################################################################
# functions to calculate basic statistics, parameter for all function is a col: column to calculate;
# return the calculated columns
################################################################################

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




################################################################################
# Functions for spatial data manipulation
# Convert a file with latitude and longitude to a `sf` object and transform the projection of that file (usualy the projected coordinate system)
# Parameters in this function are file_df : file with lat, long column and other attributes; wgs : coordinate system of the lat long format;
# UTM proj: resultant projection; lat: column name of latitude column; long: column name of longitude column
# returns a sf object with the specified UTM_proj projection
################################################################################
convert_sf_proj <- function(file_df, wgs, UTM_proj, long, lat) {  
  # Convert file to sf
  file_sf_proj <- st_as_sf(file_df, coords = c(long, lat), crs = st_crs(wgs))
  # Transform the data's projection
  file_sf_proj <- st_transform(file_sf_proj, crs = UTM_proj)
}




################################################################################
# Convert a file with latitude and longitude to a `sf` object 
# Parameters in this function are file_df : file with lat, long column and other attributes; 
# wgs : coordinate system of the lat long format;
# lat: column name of latitude column; long: column name of longitude column
# returns a sf object with wgs projection
################################################################################
convert_sf <- function(file_df, wgs, long, lat) {  
  # Convert file_df to sf
  file_sf <- st_as_sf(file_df, coords = c(long, lat), crs = st_crs(wgs))
}




################################################################################
# Convert shapefile to a particular coordinate system proj
# Parameters in this function are file_sf : a shapefile. geopackage, geojson etc, a vector data; 
# proj : the coordinate system required;
# returns a sf object with the specified proj projection
################################################################################
sf_proj <- function(file_sf, proj) {
  file_sf_proj <- st_read(file_sf, quiet = TRUE)
  file_sf_proj <- st_transform(file_sf_proj, crs = proj)
}




################################################################################
# Generate a named vector for lur buffers
# parameters are name : name of the variable to create buffer on;
# buffer_len : a vector of buffers in meters to be generated
# returns a named numeric vector with the names as "name_buffer_500m", etc.
################################################################################
lur_buffer_maker <- function(name, buffer_len) {
  # Create a list of buffer engths
  buff_labs <- paste0(name, buffer_len, "m")
  buffering <- buffer_len
  names(buffering) <- buff_labs
  buffering
}




################################################################################
# Generate a named vector for lur buffers
# parameters are file_sf_proj : projected file with latitude longitude or projected spatial points to create buffers on;
# buffer : a vector of buffers in meters to be generated
# returns a named numeric vector with the names as "name_buffer_500m", etc.
# width of the buffer is the radius from the center
# Generate buffers of the above specified length for each of the points 
################################################################################
buffer_points <- function(buffer, file_sf_proj) {
  buffers <- mapply(FUN = sf::st_buffer,
                    dist = buffer,
                    MoreArgs = list(x = file_sf_proj),
                    SIMPLIFY = FALSE, 
                    USE.NAMES = TRUE)
}




################################################################################
# Assign at once multiple returned values from a function to multiple variables
# https://stackoverflow.com/a/1829651
################################################################################
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




################################################################################
# st_distance calculates Euclidean distance to each point
# Extract distance from the airport, industries corresponding for each point
# file_sf_proj : a projected sf object usually a multipoint object
# airport : a projected sf object (point of runway)
# industries : a projected sf object (point of multiple industries)
# returns a dataframe with `inverse_distance_industries`, `inverse_distance_airport` and other attributes of the original `file_sf`
################################################################################
dist_variable_extraction <- function(file_sf_proj, airport_sf = NULL, industries = NULL) {
  file_sf_proj <- file_sf_proj %>%
    mutate(inverse_distance_industries = NA,
           inverse_distance_airport = as.numeric((1 / st_distance(., airport_sf, by_element = TRUE))))
  for(i in 1:nrow(file_sf_proj)) {
    min_dist <- which.min(st_distance(file_sf_proj[i, ], industries))
    file_sf_proj[i, "inverse_distance_industries"] <- st_distance(file_sf_proj[i, ], industries[min_dist, ])
  }
  file_sf_proj <- file_sf_proj %>% 
    mutate(inverse_distance_industries =  1 / inverse_distance_industries)
  return(file_sf_proj)
}




################################################################################
# Extract distance from the airport, industries corresponding for each point
# file_sf : a projected sf object usually a multipoint object;
# aod : Aerosol optical depth in raster format;
# dem : Digital Elevation Model in raster format;
# wgs : coordinate system of the raster files;
# UTM_proj : projected coordinate system of the file_sf or the multipoint object
# returns a dataframe with `elevation`, `aod` and other attributes of the original `file_sf`
################################################################################
extract_raster <- function(file_sf, dem, aod, wgs, UTM_proj) {
  file_sf <- st_transform(file_sf, crs = wgs)
  file_sf <- file_sf %>% 
    mutate(elevation = sqrt(raster::extract(dem, .)),
           aod = raster::extract(aod, .))
  file_sf_proj <- st_transform(file_sf, crs = UTM_proj)
}




################################################################################
# Get significance code for regression coefficients
################################################################################
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




################################################################################
# Extract the results (slope, tval, prob etc) from a linear regression model and also check slope sign for each variable
# my_model : is the linear regression model;
# sig_star : the function above to generate the significance star;
# direction_of_effect_table : it is the parameter list of the desired sign of effect 
# returns a data frame which has the model's data / parameters extracted along with whether each predictor's direction of effect is conserved or not
################################################################################
model_data_extraction <- function(my_model, sig_star, direction_of_effect_table) {
  data_extracted <- data.frame(
    slope = format((summary(my_model))$coefficients[, 1], scientific = TRUE),
    stde = round((summary(my_model))$coefficients[, 2], 4),
    tval = round((summary(my_model))$coefficients[, 3], 4),
    prob = round((summary(my_model))$coefficients[, 4], 4),
    sig = unlist(lapply((
      summary(my_model)
    )$coefficients[, 4], FUN = sig_star))
  )
  data_extracted <- cbind(value = rownames(data_extracted), data_extracted) %>% 
    filter(value != "(Intercept)")
  data_extracted <- add_sign_check(data_extracted, direction_of_effect_table)
}




################################################################################
# Function to check sign of the slope for direction of effect of each predictor from the direction_of_effect_table
# sign : it is the sign column or the condition column in the direction_of_effect_table
# slope : it is the slope column or the condition column in the direction_of_effect_table
# value : it is the value column or the condition column in the direction_of_effect_table
################################################################################
check_sign_from_table <- function(sign, slope, val) {
  if(!is.na(val)){
    do.call(sign, list(slope, val)) # Change this to 0 
  } else if(is.na(val)){
    TRUE # Change this to NA 
  }
} 




################################################################################
# Function to add sign from the existing LUR sheet and then check for sign and add an obj variable which shows sign preserved or not
# data_extracted : data extracted from model with the parameters and slope value, use this to check for the direction of effect for each parameter
# param : the column which is used to join the `data_extracted` and `direction_of_effect_table`
# `param` extracted from the parameter name string by removing the buffer length
# returns the dataframe along with the direction of effect and whether it is satisfied or not
################################################################################
add_sign_check <- function(data_extracted, direction_of_effect_table) {
  data_extracted$param <- sub("_buffer.*", "\\1", data_extracted$value)
  data_extracted <- left_join(data_extracted, direction_of_effect_table, by = "param")
  for(i in 1:nrow(data_extracted)) {
    data_extracted$obj[i] <- check_sign_from_table(data_extracted$sign[i], 
                                  as.numeric(as.character(data_extracted$slope[i])), 
                                  as.numeric(as.character(data_extracted$val[i])))
  }
  return(data_extracted)
}




################################################################################
# Function to create the linear regression model 
# data_extracted : data extracted from model with the parameters and slope value, use this to check for the direction of effect for each parameter
# model_data : dataframe with the relevant columns for the model building (includes : CODE, response variable, other parameters for the particular model)
# response_variable : response variable like PM2.5 etc
# returns the `lm` model built for the parameters; a check for only one variable is also given 
################################################################################
create_model <- function(model_data, response_variable) {
  new_data <- model_data %>% 
    dplyr::select(everything(), -!!response_variable, -"CODE", -contains(c("predicted")))
  others <- names(new_data)
  if(ncol(new_data) == 1){  
    eqtn <- as.formula(paste(response_variable, others, sep = " ~ "))
  } else { 
    eqtn <- as.formula(paste(response_variable, paste(others, collapse = "+"), sep = " ~ "))
  }
  my_model <- lm(as.formula(eqtn), data = model_data)
  return(my_model)
}




################################################################################
# function to extract buffer variables for railway 
# buffering_railway : different buffer lengths named vector 
# file_sf_proj : a projected sf object usually a multipoint object
# railway : vector data of railways usually line data
# returns a dataframe with lengths of different buffers of railways for each point of interest / unique `CODE`
################################################################################
railway_fun <- function(buffering_railway, file_sf_proj, railway) {
  buffers <- buffer_points(buffering_railway, file_sf_proj)
  buffers_railway <- mapply(FUN = sf::st_intersection,
                      x = buffers,
                      MoreArgs = list(y = railway),
                      SIMPLIFY = FALSE, 
                      USE.NAMES = TRUE)
  df_1 <- do.call(rbind, lapply(buffers_railway, as.data.frame))
  df <- cbind(buffer_m = rownames(df_1), df_1)
  df_railway <- df %>% 
    mutate(buffer_m_1 = sub(".*buffer_*(.*?) *m.*", "\\1", buffer_m)) %>% 
    st_as_sf(.) %>% 
    mutate(len = as.numeric(as.character(st_length(.)))) %>% 
    as.data.frame(.) %>% 
    dplyr::select(buffer_m_1, CODE, len) %>% 
    group_by(buffer_m_1, CODE) %>% 
    summarise_all(list(~ sum(., na.rm = TRUE))) %>% 
    pivot_wider(names_from = buffer_m_1, values_from = len)
  colnames(df_railway)[-1] <- paste0("rail_buffer_", colnames(df_railway)[-1])
  return(df_railway)
}




################################################################################
# function to run the step wise supervised linear regression 
# col_interest : all columns to go through the linear regression
# original_para : the parameter which gave the highest R2 and the desired direction of effect
# original_r2 : original_para's R2 value
# response_variable : response variable like PM2.5 etc
# direction_of_effect_table : it is the parameter list of the desired sign of effect 
# data_with_variables : data with the `CODE` column, response variable column and all the predictors / variables
# change_val : is the change in the R2 value to stop the model from further building
# returns var_list, change_r, highest_r2, data_wo_slop_a, data_r2
# var_list : variable list with the highest R2 and desired sign of effect for all parameters
# change_r : the values of change in R2 when the loop stopped
# highest_r2 : highest R2 obtained from the var_list model
# data_wo_slop_a : data of all the models built at each step
# data_r2 : data of best models at each iteration of adding variables 
################################################################################
run_lur_model <- function(col_interest, original_para, original_r2, response_variable, 
                          data_with_variables, direction_of_effect_table, change_val) {
  # Make an empty data frame and try to extract results 
  data_with <- data.frame()
  # keep track of all the variables with the right direction of effect
  data_wo_slop_a <- data.frame()
  name_col <- as.vector(col_interest)
  # remember the r2 change 
  data_r2 <- data.frame()
  # Start a variable list 
  var_list <- c(original_para)
  change_r <- 100 # positive number or number greater than 0.01 for the loop to run at least once
  thres <- change_val # Define your threshold to come out of the loop
  # Keep track of the original R2 value
  highest_r2 <- original_r2
  # Use a while loop to break when change in R2 is less than 0.01
  while(change_r >= thres) {
    # Loop for the model building and change data_with_variables to selected columns only
    for(i in name_col) {
      # Check if the parameter is already in the var list, if there then ignore it use others
      # If available ignore it use others
      if(i %in% var_list) {
        data_with <- data_with
      } else {
        others <- paste(paste(var_list, collapse = " + "), i, sep = " + ")
        # Generate a formula each time
        equ <- as.formula(paste(response_variable, others, sep = " ~ "))
        # Apply multiple linear regression
        lm_step <- lm(equ, data = data_with_variables)
        vr <- round((summary(lm_step))$coefficients[, 4], 4)
        if(all(!is.nan(vr)) & all(!is.na(vr))) {
          # Apply model_data_extraction function to extract slope, std error, t value, significance etc and also check the slope sign
          data_corr <- model_data_extraction(lm_step, sig_star, direction_of_effect_table)
          # Extract R2, AIC, RMSE and also track the equation
          data_corr <- data_corr %>% 
            mutate(r2 = round(summary(lm_step)$adj.r.squared, 4), aic = round(AIC(lm_step), 4), 
                   rmse = summary(lm_step)$sigma, eqtn = others)
          data_with <- rbind(data_with, data_corr)
        } else {
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
  return(list(var_list, change_r, highest_r2, data_wo_slop_a, data_r2))
}




################################################################################
# function to remove variables with p value greater than 0.1
# my_model : is the linear regression model
# data_with_selected_variables : data frame with the selected parameters in my_model, CODE, and the response variable
# no is the value of p for removing from the data
# returns a data frame similar to data_with_selected_variables after removing the variables which have p value greater than 0.1
################################################################################
remove_p_value <- function(my_model, data_with_selected_variables, no) {
  p_values <- round((summary(my_model))$coefficients[, 4], 4)
  p_remove <- names(p_values)[(p_values > no)]
  print(c("The columns removed due to p-value greater than ", no, " are ", p_remove))
  if(!identical(p_remove, character(0))) {
    data_cleaned <- data_with_selected_variables[ , - which(names(data_with_selected_variables) %in% p_remove)]
  } else {
    data_cleaned <- data_with_selected_variables
  }
  return(data_cleaned)
}




################################################################################
# function to remove variables with influential vif
# my_model : is the linear regression model
# data_with_selected_variables : data frame with the selected parameters in my_model, CODE, and the response variable
# no is the value of vif for removing from the data
# returns a data frame similar to data_with_selected_variables after removing the influential variables or vif > 3
################################################################################
vif_function <- function(my_model, data_with_selected_variables, no) {
  vif_variable <- vif(my_model)
  print(c("The values of vif are ", vif_variable))
  # Find the columns with high VIF and remove them from the dataframe
  influential <- names(vif_variable)[(vif_variable > no)]
  if(!identical(influential, character(0))) {
    data_cleaned <- data_with_selected_variables[ , -which(names(data_with_selected_variables) %in% influential)]
  } else {
    data_cleaned <- data_with_selected_variables
  }
  return(data_cleaned)
}




################################################################################
# function to add column names and group some groups of the Land use classification 
# df_lulc : dataframe with different land uses - 10, 20, 30... , and `CODE` column and the corresponding buffers
# returns a dataframe which is grouped and has `CODE` and other land use variables with the buffers
################################################################################
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




################################################################################
# function to carry out leave one out cross validation (loocv)
# data_final_selected_variables : data frame with `CODE`, response variable, and the final model selected parameters only
# response_variable : response variable like the PM2.5 etc
# returns a data frame with the `CODE` column and loocv r2, adjusted r2, and loocv prediction
################################################################################
loocv_loop <- function(data_final_selected_variables, response_variable) {
  data_final_selected_variables <- data_final_selected_variables %>% 
    mutate(predicted_loocv = NA, 
           predicted_loocv_r2 = NA,
           predicted_loocv_r2_adj = NA)
  for(i in 1:nrow(data_final_selected_variables)) {
    data_train <- data_final_selected_variables[-i, ]
    mdl <- create_model(data_train, response_variable)
    data_final_selected_variables$predicted_loocv[i] <- predict(mdl, data.frame(data_final_selected_variables[i, ]))
    data_final_selected_variables$predicted_loocv_r2[i] <- summary(mdl)$r.squared
    data_final_selected_variables$predicted_loocv_r2_adj[i] <- summary(mdl)$adj.r.squared
  }
  return(data_final_selected_variables)
}

