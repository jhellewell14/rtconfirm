#' Run the back-calculation model
#'
#' @param country String defining country to run model on (default = South Korea)
#'
#' @return List of a plot of observed confirmed cases (black) and estimated incidence of infection (red ribbon) and the fit of the stan model object returned as a stan object
#' @export
#' @useDynLib rtconfirm, .registration=TRUE
#' @importFrom NCoVUtils get_ecdc_cases
#' @importFrom dplyr filter
#' @importFrom rstan sampling extract
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot geom_line geom_ribbon xlab ylab aes
#' @examples
#' results <- run_backcalc()
run_backcalc <- function(country = "South_Korea") {

 sk <- NCoVUtils::get_ecdc_cases(countries = country)

 sk <- sk %>%
  dplyr::filter(date > "2020-01-31")

 dat <- list(t = nrow(sk),
             y = sk$cases,
             inc_loc = 1.621,
             inc_scale = 0.418,
             delay_alpha = 1.3218,
             delay_beta  = 0.2359536,
             tau = 7)


 model <- stanmodels$backcalc

 fit <- rstan::sampling(mod, data = dat, chains = 4)


 res <- rstan::extract(fit)

 p <- data.frame(obs = dat$y,
            date = sk$date,
            # inf_med = res$infections#infections,
            inf_low = apply(res$infections, 2, quantile, prob = 0.025),
            inf_up = apply(res$infections, 2, quantile, prob = 0.975)) %>%
  ggplot2::ggplot(ggplot2::aes(x = date)) +
  ggplot2::geom_line(ggplot2::aes(y = obs)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = inf_low, ymax = inf_up),
                       fill = "red", alpha = 0.4) +
  cowplot::theme_cowplot() +
  ggplot2::ylab("Cases") +
  ggplot2::xlab("")

 return(list("plot" = p, "fit" = res))
}
