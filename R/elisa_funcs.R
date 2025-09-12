## ----- File containing all the functions to process ELISA data --------

# ======= Fit standard curve functions =======
# standardize data format before analysis
standardize_data <- function(df,
                             id_col = "HCDC_SAMPLE_ID",
                             result_col = "RESULT",
                             dilution_fct_col = "DILUTION_FACTORS",
                             antitoxin_label = "Anti_toxin"
                             ){
  df <- df %>%
    rename(
      sample_id := id_col,
      result := result_col,
      dilution_factors := dilution_fct_col
    ) %>%
    clean_names() %>%
    mutate(
      sample_id = if_else(sample_id == antitoxin_label, "antitoxin", sample_id)
    )

  df
}

# Helper function to do 2 things
# - get antitoxins data
# - compute concentration from reference concentration and dilution factor
get_antitoxins <- function(plate, ref_conc = 10) {
  plate |>
    filter(sample_id == "antitoxin") |>
    mutate(concentration = ref_conc / dilution_factors)
}

# Helper function to look at the data and generate a good initial guess for the parameters
good_guess4PL <- function(x, y, eps = .3) {
  nb_rep <- unique(table(x))
  the_order <- order(x)
  x <- x[the_order]
  y <- y[the_order]
  a <- min(y)
  d <- max(y)
  c <- approx(y, x, (d - a) / 2, ties = "ordered")$y
  list(a = a, c = c, d = d,
       b = (
         approx(x, y, c + eps, ties = "ordered")$y -
           approx(x, y, c - eps, ties = "ordered")$y
       ) / (2 * eps))
}


# function to fit data to a 4PL model
nls4PL <- function(df) {
  nls(result ~ d + (a - d) / (1 + 10^((log10(concentration) - c) * b)),
      data = df,
      start = with(df, good_guess4PL(log10(concentration), result)))
}

rowsplit <- function(df) split(df, 1:nrow(df))

# helper function to generate data to quantify uncertainty (using bootstrapping)
# general flow is:
# - sample different values for parameters (number of samples = nb)
# - use model to compute OD for the new set of parameter values
simulate_nls <- function(object, newdata, nb = 9999) {
  nb |>
    mvtnorm::rmvnorm(mean = coef(object), sigma = vcov(object)) |>
    as.data.frame() |>
    rowsplit() |>
    map(as.list) |>
    map(~ c(.x, newdata)) |>
    map_dfc(eval, expr = parse(text = as.character(formula(object))[3]))
}


#' Function to generate the standard curve and its confidence interval
#'
#' @param df - result for antitoxin
#' @param model
#' @param le
#' @param level
#' @param nb
#'
#' @return
#' @export
#'
#' @examples
standard_curve_data <- function(df, model,
                                le = 512, level = .95, nb = 9999) {

  log_concentration <- log10(df$concentration)
  logc <- seq(min(log_concentration), max(log_concentration), le = le)
  alpha <- (1 - level) / 2
  df |>
    model() |>
    simulate_nls(list(concentration = 10^logc), nb) |>
    apply(1, quantile, c(alpha, .5, 1 - alpha)) |>
    t() |> as.data.frame() |>
    setNames(c("lower", "median", "upper")) |>
    (\(.x) bind_cols(logc = logc, .x))()
}

# Function to turn the standard curve data.frame to a function that coverts OD to LC
data2function <- function(df) {
  with(df, {
    approxfun2 <- function(...) approxfun(y = logc, ...)
    pred_lwr <- approxfun2(upper)
    pred_mdi <- approxfun2(median)
    pred_upp <- approxfun2(lower)
    function(x) c(lower = pred_lwr(x), median = pred_mdi(x), upper = pred_upp(x))
  })
}

# Function to convert samples' OD to LC
process_samples <- function(plate, std_crv) {
  plate %>%
    filter(sample_id != "antitoxin") %>%
    bind_cols(map_dfr(.$result, std_crv))
}


# ===== Quality Control functions ========
# only use dilution factor where negative sample is indeed negative
get_negative_controls <- function(plate, std_crv){
  plate %>%
    select(starts_with("negative")) %>%
    unique() %>%
    tidyr::pivot_longer(everything(), names_to = "dilution", values_to = "od") %>%
    mutate(across(dilution, ~ stringr::str_remove(.x, "negative_") %>%  as.integer()))
}

label_positive <- function(plate, std_curve, positive_threshold = 0.1){
  plate %>%
    bind_cols(map_dfr(.$od, std_curve)) %>%
    mutate(
      positive = upper >= 0.1
    )
}

# ======== Plot functions ==========
plot_standard_curve <- function(scdf, data=NULL, ylim=NULL, datapoint_size = 3){
  if (is.null(ylim)) ylim <- c(0, max(scdf$upper, data$result))

  plot <- ggplot() +
    geom_ribbon(
      aes(
        x = logc, y = median, ymin = lower, ymax = upper
      ),
      fill = adjustcolor("red", alpha.f = 0.2), data = scdf
    ) +
    geom_line(
      aes(
        x = logc, y = median
      ),
      color = "red", lwd = 0.5, data = scdf
    ) +
    labs(
      x = "log10(concentration)",
      y = "Optical density"
    )

  if (!is.null(data)){
    plot <- plot +
      geom_point(
        aes(
          x = log10(concentration),
          y = result
        ),
        shape = 3, size = datapoint_size,
        color = "blue",
        data = data
      )
  }
  plot
}

# a ggplot layer for plotting positive threshold at each dilution
add_dilutions <- function(dilution_factors, positive_threshold = 0.1) {
  geom_vline(
    aes(
      xintercept = log10(positive_threshold / c(1, dilution_factors))
    ),
    color = "green")
}

