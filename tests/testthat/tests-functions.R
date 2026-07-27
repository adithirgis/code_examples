pacman::p_load(testthat, sf, raster, dplyr)


test_that("GM1 computes the geometric mean correctly", {
  expect_equal(GM1(c(1, 2, 4, 8)), exp(mean(log(c(1, 2, 4, 8)))))
  expect_equal(GM1(c(2, 2, 2, 2)), 2)
  expect_equal(GM1(c(1, 2, 4)), 2, tolerance = 1e-8)
})

test_that("GM1 handles NA according to na_rm", {
  x <- c(1, 2, 4, NA)
  expect_equal(GM1(x, na_rm = TRUE), GM1(c(1, 2, 4)))
  expect_true(is.na(GM1(x, na_rm = FALSE)))
})

test_that("GSD1 computes the geometric standard deviation correctly", {
  x <- c(1, 2, 4, 8)
  expect_equal(GSD1(x), exp(sd(log(x))))
  # GSD of identical values is 1 (no spread)
  expect_equal(GSD1(c(3, 3, 3)), 1)
})

test_that("GSD1 handles NA according to na_rm", {
  x <- c(1, 2, 4, NA)
  expect_equal(GSD1(x, na_rm = TRUE), GSD1(c(1, 2, 4)))
  expect_true(is.na(GSD1(x, na_rm = FALSE)))
})

test_that("CV1 computes the coefficient of variation correctly", {
  x <- c(2, 4, 4, 4, 5, 5, 7, 9)
  expect_equal(CV1(x), sd(x) / mean(x))
  expect_equal(CV1(c(5, 5, 5, 5)), 0)
})

test_that("CV1 handles NA according to na_rm", {
  x <- c(2, 4, 4, NA)
  expect_equal(CV1(x, na_rm = TRUE), CV1(c(2, 4, 4)))
  expect_true(is.na(CV1(x, na_rm = FALSE)))
})

test_that("stderr computes the standard error of the mean", {
  x <- c(2, 4, 4, 4, 5, 5, 7, 9)
  expect_equal(stderr(x), sqrt(var(x) / length(x)))
})

test_that("stderr respects na.rm", {
  x <- c(1, 2, 3, NA)
  expect_true(is.na(stderr(x, na.rm = FALSE)))
  expect_equal(stderr(x, na.rm = TRUE), sqrt(var(c(1, 2, 3)) / 3))
})

test_that("sig_star assigns the correct stars by p-value threshold", {
  expect_equal(sig_star(0.0001), "***")
  expect_equal(sig_star(0.0009), "***")
  expect_equal(sig_star(0.005), "**")
  expect_equal(sig_star(0.02), "*")
  expect_equal(sig_star(0.08), ".")
  expect_equal(sig_star(0.5), " ")
  expect_equal(sig_star(1), " ")
})

test_that("sig_star boundary values fall into the stricter (lower) bucket", {
  expect_equal(sig_star(0.001), "**")
  expect_equal(sig_star(0.01), "*")
  expect_equal(sig_star(0.05), ".")
  expect_equal(sig_star(0.1), " ")
})

test_that("sig_star's documented NA/NULL/blank handling does not match its behavior", {
  ## The roxygen docs claim NA, NULL, and " " all return " ", but the
  ## `x < 0.001` comparison is evaluated before any of the is.na()/is.null()
  ## checks are reached, so these inputs behave unexpectedly in practice.
  ## These tests document the *actual* current behavior so a future fix
  ## is caught by a red test rather than silently changing behavior.
  expect_error(sig_star(NA))     ## NA < 0.001 is NA -> if(NA) errors
  expect_error(sig_star(NULL))   ## logical(0) -> if(logical(0)) errors
  ## " " < 0.001 coerces 0.001 to character and does a *string* comparison,
  ## so this incorrectly returns "***" instead of the documented " ".
  expect_equal(sig_star(" "), "***")
})

test_that(":= assigns the whole rhs list, unwrapped, in the single-variable case", {
  ## Note: when only one variable is named on the left, the current
  ## implementation assigns the raw `rhs` (still a list) rather than its
  ## single element - `a` ends up as `list(42)`, not `42`.
  c(a) := list(42)
  expect_type(a, "list")
  expect_equal(a, list(42))
  expect_equal(a[[1]], 42)
})

test_that(":= assigns multiple values from a list in order", {
  c(a, b, c_) := list(1, "two", TRUE)
  expect_equal(a, 1)
  expect_equal(b, "two")
  expect_equal(c_, TRUE)
})

test_that(":= pads missing right-hand values with NULL", {
  c(a, b, c_) := list(1, 2)
  expect_equal(a, 1)
  expect_equal(b, 2)
  expect_null(c_)
})

make_model_data <- function(n = 30, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)                       # noise predictor, weak/no relationship
  y  <- 2 + 3 * x1 - 1.5 * x2 + rnorm(n, sd = 0.2)
  data.frame(
    CODE = seq_len(n),
    PM2.5 = y,
    tree_cover_buffer_500 = x1,
    industry_buffer_500 = x2,
    noise_buffer_500 = x3
  )
}

test_that("create_model builds an lm using all predictors except CODE/response/predicted*", {
  df <- make_model_data()
  m <- create_model(df, "PM2.5")
  expect_s3_class(m, "lm")
  expect_setequal(
    all.vars(formula(m))[-1],
    c("tree_cover_buffer_500", "industry_buffer_500", "noise_buffer_500")
  )
})

test_that("create_model excludes columns containing 'predicted'", {
  df <- make_model_data()
  df$predicted_loocv <- 1
  m <- create_model(df, "PM2.5")
  expect_false("predicted_loocv" %in% all.vars(formula(m)))
})

test_that("create_model builds a single-predictor formula correctly", {
  df <- make_model_data()
  df$industry_buffer_500 <- NULL
  df$noise_buffer_500 <- NULL
  m <- create_model(df, "PM2.5")
  expect_equal(all.vars(formula(m)), c("PM2.5", "tree_cover_buffer_500"))
})

test_that("create_model validates its inputs", {
  df <- make_model_data()
  expect_error(create_model(as.list(df), "PM2.5"), "data.frame")
  expect_error(create_model(df, "not_a_column"), "not a column")
})

test_that("create_model errors when no predictors remain", {
  df <- data.frame(CODE = 1:5, PM2.5 = rnorm(5))
  expect_error(create_model(df, "PM2.5"), "No predictor variables remain")
})

test_that("remove_p_value drops columns whose predictor p-value exceeds the threshold", {
  df <- make_model_data()
  m <- create_model(df, "PM2.5")
  pvals <- summary(m)$coefficients[, 4]
  expect_gt(pvals[["noise_buffer_500"]], 0.05)
  cleaned <- remove_p_value(m, df, p_threshold = 0.05)
  expect_false("noise_buffer_500" %in% names(cleaned))
  expect_true(all(c("tree_cover_buffer_500", "industry_buffer_500") %in% names(cleaned)))
  expect_equal(nrow(cleaned), nrow(df))
})

test_that("remove_p_value leaves data untouched when nothing exceeds the threshold", {
  df <- make_model_data()
  m <- create_model(df, "PM2.5")
  expect_message(
    cleaned <- remove_p_value(m, df, p_threshold = 0.999),
    "No predictors exceeded"
  )
  expect_equal(names(cleaned), names(df))
})

test_that("remove_p_value validates its inputs", {
  df <- make_model_data()
  m <- create_model(df, "PM2.5")
  expect_error(remove_p_value("not a model", df, 0.05), "lm object")
  expect_error(remove_p_value(m, list(), 0.05), "data frame")
  expect_error(remove_p_value(m, df, 1.5), "between 0 and 1")
  expect_error(remove_p_value(m, df, -0.1), "between 0 and 1")
})

test_that("vif_function drops predictors above the VIF threshold", {
  set.seed(2)
  n <- 50
  x1 <- rnorm(n)
  x2 <- x1 * 0.98 + rnorm(n, sd = 0.01)  ## near-collinear with x1 -> high VIF
  x3 <- rnorm(n)
  df <- data.frame(
    CODE = seq_len(n),
    PM2.5 = 2 + x1 + x3 + rnorm(n, sd = 0.1),
    a_buffer_500 = x1,
    b_buffer_500 = x2,
    c_buffer_500 = x3
  )
  m <- create_model(df, "PM2.5")
  cleaned <- suppressMessages(vif_function(m, df, vif_threshold = 5))
  expect_true(all(c("a_buffer_500", "b_buffer_500") %in% names(df)))
  ## at least one of the collinear pair should be removed
  expect_lt(length(intersect(c("a_buffer_500", "b_buffer_500"), names(cleaned))), 2)
  expect_true("c_buffer_500" %in% names(cleaned))
})

test_that("vif_function keeps data unchanged when no VIF exceeds the threshold", {
  df <- make_model_data()
  m <- create_model(df, "PM2.5")
  expect_message(
    cleaned <- vif_function(m, df, vif_threshold = 100),
    "No predictors exceeded"
  )
  expect_equal(names(cleaned), names(df))
})

test_that("vif_function validates its inputs", {
  df <- make_model_data()
  m <- create_model(df, "PM2.5")
  expect_error(vif_function("not a model", df, 5), "lm object")
  expect_error(vif_function(m, list(), 5), "data frame")
  expect_error(vif_function(m, df, -1), "positive numeric")
  # single-predictor model: VIF is undefined
  df2 <- df[, c("CODE", "PM2.5", "tree_cover_buffer_500")]
  m2 <- create_model(df2, "PM2.5")
  expect_error(vif_function(m2, df2, 5), "at least two predictors")
})

test_that("check_sign_from_table returns TRUE when val is NA (no constraint)", {
  expect_true(check_sign_from_table("<", slope = 5, val = NA))
})

test_that("check_sign_from_table evaluates the named comparison correctly", {
  expect_true(check_sign_from_table("<", slope = -2, val = 0))
  expect_false(check_sign_from_table("<", slope = 2, val = 0))
  expect_true(check_sign_from_table(">", slope = 2, val = 0))
  expect_false(check_sign_from_table(">", slope = -2, val = 0))
  expect_true(check_sign_from_table(">=", slope = 0, val = 0))
  expect_true(check_sign_from_table("<=", slope = 0, val = 0))
})

test_that("check_sign_from_table rejects an unrecognized comparison operator", {
  expect_error(check_sign_from_table("!=", slope = 1, val = 0), "must be one of")
})

test_that("add_sign_check derives param from value and flags direction correctly", {
  data_extracted <- tibble::tibble(
    value = c("tree_cover_buffer_500", "industry_buffer_1000"),
    slope = c("-1.2e+00", "3.4e+00")
  )
  doe <- tibble::tibble(
    param = c("tree_cover", "industry"),
    sign  = c("<", ">"),
    val   = c(0, 0)
  )
  out <- add_sign_check(data_extracted, doe)
  expect_equal(out$param, c("tree_cover", "industry"))
  expect_true(all(out$obj))
})

test_that("add_sign_check flags a violated direction of effect as FALSE", {
  data_extracted <- tibble::tibble(
    value = c("tree_cover_buffer_500"),
    slope = c("2.0e+00")
  )
  doe <- tibble::tibble(param = "tree_cover", sign = "<", val = 0)
  out <- add_sign_check(data_extracted, doe)
  expect_false(out$obj)
})

test_that("add_sign_check validates required columns", {
  doe <- tibble::tibble(param = "x", sign = "<", val = 0)
  expect_error(add_sign_check(data.frame(value = "x"), doe), "missing required columns")
  expect_error(add_sign_check(data.frame(value = "x", slope = "1"), data.frame(param = "x")),
               "missing required")
})

test_that("extract_model_data returns coefficient stats joined with direction-of-effect flags", {
  df <- make_model_data()
  m <- create_model(df, "PM2.5")
  doe <- tibble::tibble(
    param = c("tree_cover", "industry", "noise"),
    sign  = c("<", ">", "is.na"),
    val   = c(0, 0, NA)
  )
  out <- extract_model_data(m, sig_star, doe)
  expect_false("(Intercept)" %in% out$value)
  expect_equal(nrow(out), 3)
  expect_true(all(c("slope", "stde", "tval", "prob", "sig", "obj") %in% names(out)))
  ## tree_cover's true slope is positive (+3) but doe expects "<" 0 -> should be FALSE
  expect_false(out$obj[out$value == "tree_cover_buffer_500"])
  ## industry's true slope is negative (-1.5) but doe expects ">" 0 -> should be FALSE
  expect_false(out$obj[out$value == "industry_buffer_500"])
  ## noise has no constraint (is.na is not evaluated because val is NA) -> TRUE
  expect_true(out$obj[out$value == "noise_buffer_500"])
})

test_that("extract_model_data validates its inputs", {
  df <- make_model_data()
  m <- create_model(df, "PM2.5")
  doe <- tibble::tibble(param = "x", sign = "<", val = 0)
  expect_error(extract_model_data("not a model", sig_star, doe), "lm object")
  expect_error(extract_model_data(m, "not a function", doe), "must be a function")
  expect_error(extract_model_data(m, sig_star, data.frame(param = "x")), "must contain columns")
})

test_that("loop_loocv adds one held-out prediction per row", {
  df <- make_model_data(n = 12)
  out <- loop_loocv(df, "PM2.5")
  expect_equal(nrow(out), nrow(df))
  expect_true(all(c("predicted_loocv", "predicted_loocv_r2", "predicted_loocv_r2_adj") %in% names(out)))
  expect_false(any(is.na(out$predicted_loocv)))
})

test_that("loop_loocv validates its inputs", {
  df <- make_model_data(n = 12)
  expect_error(loop_loocv(list(), "PM2.5"), "data frame")
  expect_error(loop_loocv(df, "not_a_col"), "must be a column")
  expect_error(loop_loocv(df[1, ], "PM2.5"), "at least two observations")
})

test_that("loop_kfold returns one row per observation with fold assignments", {
  df <- make_model_data(n = 20)
  out <- loop_kfold(df, "PM2.5", k = 5, seed = 42)
  expect_equal(nrow(out), nrow(df))
  expect_true(all(c("predicted_10fold", "fold") %in% names(out)))
  expect_equal(sort(unique(out$fold)), 1:5)
})

test_that("loop_kfold is reproducible when given the same seed", {
  df <- make_model_data(n = 20)
  out1 <- loop_kfold(df, "PM2.5", k = 5, seed = 123)
  out2 <- loop_kfold(df, "PM2.5", k = 5, seed = 123)
  expect_equal(out1[order(out1$CODE), "predicted_10fold"],
               out2[order(out2$CODE), "predicted_10fold"])
})

test_that("loop_kfold validates k", {
  df <- make_model_data(n = 10)
  expect_error(loop_kfold(df, "PM2.5", k = 1), "greater than or equal to 2")
  expect_error(loop_kfold(df, "PM2.5", k = 100), "cannot be greater than")
})

wgs <- "+proj=longlat +datum=WGS84 +no_defs"
utm <- "+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs"

################################################################################
# convert_sf()
################################################################################
test_that("convert_sf returns an sf object", {
  df <- data.frame(
    lon = c(77, 77.1),
    lat = c(12, 12.1),
    id = 1:2
  )
  x <- convert_sf(df, wgs, "lon", "lat")
  expect_s3_class(x, "sf")
  expect_equal(nrow(x), 2)
})

test_that("convert_sf errors if coordinate columns are missing", {
  df <- data.frame(
    lon = 1,
    value = 2
  )
  expect_error(
    convert_sf(df, wgs, "lon", "lat"),
    "column\\(s\\) not found"
  )
})

test_that("convert_sf errors if coordinates contain NA", {
  df <- data.frame(
    lon = c(77, NA),
    lat = c(12, 13)
  )
  expect_error(
    convert_sf(df, wgs, "lon", "lat"),
    "row\\(s\\) have NA"
  )
})

################################################################################
# convert_sf_proj()
################################################################################
test_that("convert_sf_proj projects coordinates", {
  df <- data.frame(
    lon = c(77, 77.1),
    lat = c(12, 12.1)
  )
  x <- convert_sf_proj(df, wgs, utm, "lon", "lat")
  expect_s3_class(x, "sf")
  expect_false(st_is_longlat(x))
})

################################################################################
# lur_buffer_maker()
################################################################################
test_that("lur_buffer_maker creates named buffers", {
  b <- lur_buffer_maker("rail_buffer_", c(500, 1000))
  expect_equal(
    names(b),
    c("rail_buffer_500m", "rail_buffer_1000m")
  )
  expect_equal(unname(b), c(500, 1000))
})

################################################################################
# buffer_points()
################################################################################
test_that("buffer_points returns one buffer per distance", {
  pts <- st_as_sf(
    data.frame(x = 0, y = 0),
    coords = c("x","y"),
    crs = 32643
  )
  buffers <- c(
    rail_buffer_500m = 500,
    rail_buffer_1000m =1000
  )
  out <- buffer_points(buffers, pts)
  expect_length(out, 2)
  expect_true(all(sapply(out, inherits, "sf")))
})

test_that("buffer_points requires numeric buffers", {
  pts <- st_as_sf(
    data.frame(x = 0, y = 0),
    coords = c("x","y"),
    crs = 32643
  )
  expect_error(
    buffer_points(c("A"),pts),
    "numeric"
  )
})

test_that("buffer_points rejects negative buffers", {
  pts <- st_as_sf(
    data.frame(x = 0, y = 0),
    coords=c("x","y"),
    crs = 32643
  )
  expect_error(
    buffer_points(c(a = -1), pts),
    "negative"
  )
})

test_that("buffer_points requires names", {
  pts <- st_as_sf(
    data.frame(x = 0, y = 0),
    coords = c("x","y"),
    crs = 32643
  )
  expect_error(
    buffer_points(c(100), pts),
    "fully named"
  )
})

################################################################################
# extract_dist_variable()
################################################################################
test_that("extract_dist_variable adds inverse-distance variables", {
  pts <- st_as_sf(
    data.frame(id = 1, x = 0, y = 0),
    coords=c("x","y"),
    crs = 32643
  )
  airport <- st_as_sf(
    data.frame(id = 1, x = 100, y = 0),
    coords = c("x","y"),
    crs = 32643
  )
  industry <- st_as_sf(
    data.frame(id = 1, x = 50, y = 0),
    coords = c("x","y"),
    crs = 32643
  )
  out <- extract_dist_variable(
    pts,
    airport,
    industry
  )
  expect_true(
    all(c(
      "inverse_distance_airport",
      "inverse_distance_industries"
    ) %in% names(out))
  )
})

test_that("extract_dist_variable rejects non-sf input", {
  expect_error(
    extract_dist_variable(data.frame()),
    "sf object"
  )
})

################################################################################
# extract_raster()
################################################################################
test_that("extract_raster extracts raster values", {
  r1 <- raster(
    nrows = 10,
    ncols = 10,
    xmn = 76,
    xmx = 78,
    ymn = 11,
    ymx = 13,
    crs = wgs
  )
  values(r1) <- 1:100
  r2 <- r1
  values(r2) <- 100:1
  pts <- st_as_sf(
    data.frame(
      lon = 77,
      lat = 12
    ),
    coords = c("lon","lat"),
    crs = wgs
  )
  pts <- st_transform(pts,utm)
  out <- extract_raster(
    pts, r1, r2, wgs, utm
  )
  expect_true("elevation" %in% names(out))
  expect_true("aod" %in% names(out))
})

################################################################################
# extract_railway()
################################################################################
test_that("extract_railway returns railway lengths", {
  pts <- st_as_sf(
    data.frame(
      CODE = "A",
      x = 0,
      y = 0
    ),
    coords = c("x","y"),
    crs = 32643
  )
  railway <- st_sf(
    CODE = "rail",
    geometry = st_sfc(
      st_linestring(
        matrix(
          c(
            -500,0, 500,0
          ),
          byrow = TRUE, ncol = 2
        )
      ),
      crs = 32643
    )
  )
  buffers <- lur_buffer_maker(
    "rail_buffer_",
    500
  )
  out <- extract_railway(
    buffers,
    pts,
    railway
  )
  expect_true("CODE" %in% names(out))
  expect_true("rail_buffer_500" %in% names(out))
})

################################################################################
# proj_sf()
################################################################################
test_that("proj_sf errors for missing file", {
  expect_error(
    proj_sf("missing.shp", 32643),
    "file not found"
  )
})

################################################################################
# derive_lulc_as_df()
################################################################################
test_that("derive_lulc_as_df combines shrubland classes", {
  x <- tibble::tibble(
    CODE = "A",
    `20_buffer_500m` = 1,
    `30_buffer_500m` = 2,
    `90_buffer_500m` = 3
  )
  out <- derive_lulc_as_df(x)
  expect_equal(
    out$shrubland_buffer_500,
    6
  )
})

test_that("derive_lulc_as_df renames variables", {
  x <- tibble::tibble(
    CODE = "A",
    `10_buffer_500m` = 4,
    `40_buffer_500m` = 8
  )
  out <- derive_lulc_as_df(x)
  expect_true("tree_cover_buffer_500" %in% names(out))
  expect_true("cropland_buffer_500" %in% names(out))
})

test_that("derive_lulc_as_df handles multiple rows", {
  x <- tibble::tibble(
    CODE = c("A", "B"),
    `10_buffer_500m` = c(1, 2)
  )
  out <- derive_lulc_as_df(x)
  expect_equal(nrow(out), 2)
})

test_that("derive_lulc_as_df ignores NA", {
  x <- tibble::tibble(
    CODE = "A",
    `20_buffer_500m` = NA,
    `30_buffer_500m` = 5
  )
  out <- derive_lulc_as_df(x)
  expect_equal(out$shrubland_buffer_500, 5)
})

test_that("derive_lulc_as_df requires data.frame", {
  expect_error(
    derive_lulc_as_df(1),
    "must be a data frame"
  )
})

test_that("derive_lulc_as_df requires CODE column", {
  expect_error(
    derive_lulc_as_df(data.frame(x = 1)),
    "CODE"
  )
})

################################################################################
# write_lulc_vars()
################################################################################
test_that("write_lulc_vars returns invisibly", {
  tmp <- tempfile(fileext = ".csv")
  dat <- tibble::tibble(
    CODE = "A",
    variable = "10_buffer_500m",
    Area = 1)
  expect_invisible(
    write_lulc_vars(dat, tmp)
  )
})

test_that("write_lulc_vars checks columns", {
  expect_error(
    write_lulc_vars(data.frame(a = 1), "a.csv"),
    "Missing required columns"
  )
})

test_that("write_lulc_vars checks filename", {
  dat <- tibble::tibble(
    CODE = "A",
    variable = "x",
    Area = 1
  )
  expect_error(
    write_lulc_vars(dat, 1),
    "single file path"
  )
})

################################################################################
# extract_road_filtered()
################################################################################
test_that("extract_road_filtered sums road lengths", {
  dat <- tibble::tibble(
    CODE = c("A", "A", "B"),
    buffer_m_1 = c("500", "500", "100"),
    len = c(10, 20, 5)
  )
  out <- extract_road_filtered(dat)
  expect_equal(out$`500`[1], 30)
  expect_equal(out$`100`[2], 5)
})

test_that("extract_road_filtered ignores NA", {
  dat <- tibble::tibble(
    CODE = "A",
    buffer_m_1 = "500",
    len = NA
  )
  out <- extract_road_filtered(dat)
  expect_equal(out$`500`, 0)
})

################################################################################
# summarise_roads()
################################################################################
test_that("summarise_roads filters classes", {
  dat <- tibble::tibble(
    CODE = c("A", "A"),
    buffer_m_1 = c("500", "500"),
    len = c(10, 20),
    fclass = c("primary", "residential")
  )
  out <- summarise_roads(
    dat,
    classes = "primary",
    prefix = "roads_"
  )
  expect_equal(out$roads_500, 10)
})

test_that("summarise_roads prefixes names", {
  dat <- tibble::tibble(
    CODE = "A",
    buffer_m_1 = "500",
    len = 4,
    fclass = "primary"
  )
  out <- summarise_roads(
    dat,
    prefix = "roads_"
  )
  expect_true(
    "roads_500" %in% names(out)
  )
})

test_that("summarise_roads checks required columns", {
  expect_error(
    summarise_roads(
      data.frame(),
      prefix = "roads_"
    ),
    "Missing required columns"
  )
})

################################################################################
# tune_rf()
################################################################################
skip_if_not_installed("ranger")

train <- data.frame(
  y = rnorm(40),
  x1 = rnorm(40),
  x2 = rnorm(40)
)
grid <- data.frame(
  num.trees = 20,
  splitrule = "variance",
  mtry = 2,
  max.depth = 5,
  min.node.size = 1,
  replace = TRUE,
  sample.fraction = .8
)

test_that("tune_rf returns metrics", {
  out <- tune_rf(grid, y ~ x1 + x2, train, "y")
  expect_true("OOB_RMSE" %in% names(out))
  expect_true("OOB_R2" %in% names(out))
  expect_true("valid_RMSE" %in% names(out))
  expect_true("valid_R2" %in% names(out))
})

test_that("character formula works", {
  out <- tune_rf(grid, "y ~ x1 + x2", train, "y")
  expect_s3_class(out, "data.frame")
})

test_that("missing response errors", {
  expect_error(
    tune_rf(grid, y ~ x1, data.frame(x1 = 1), "y")
  )
})

test_that("missing hyperparameter column errors", {
  expect_error(
    tune_rf(
      data.frame(a = 1), y ~ x, data.frame(y = 1, x = 1), "y"
    )
  )
})

test_that("seed is reproducible", {
  out1 <- tune_rf(grid, y ~ x1 + x2, train, "y", seed = 12)
  out2 <- tune_rf(grid, y ~ x1 + x2, train, "y", seed = 12)
  expect_equal(out1,out2)
})