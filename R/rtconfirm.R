#' Run the Rt estimate model
#'
#' @param country String defining country to run model on (default = South Korea)
#' @param cut integer for the number of days to shift Rt by confirmation date back by to get Rt by infection date (default = 10 days)
#' @return List of a plot of observed confirmed cases (black) and estimated infections (red), a plot of the reproduction number by date of infection (blue), and the fit of the stan model object returned as a stan object
#' @export
#' @useDynLib rtconfirm, .registration=TRUE
#' @importFrom NCoVUtils get_ecdc_cases
#' @importFrom dplyr filter
#' @importFrom rstan sampling extract
#' @importFrom ggplot2 ggplot geom_line geom_ribbon geom_hline xlab ylab scale_x_date coord_cartesian geom_bar
#' @importFrom cowplot theme_cowplot
#'
#' @examples
#' plots <- run_rtconfirm()
run_rtconfirm <- function(country = "South_Korea", cut = 10) {

 sk <- NCoVUtils::get_ecdc_cases(countries = country)

 sk <- sk %>%
  dplyr::filter(date > "2020-01-31")

 dat_test <- list(t = nrow(sk),
                  tau = 7,
                  si_loc = 1.3780732,
                  si_scale = 0.6184616,
                  inc_loc = 1.621,
                  inc_scale = 0.418,
                  delay_alpha = 1.3218,
                  delay_beta  = 0.2359536,
                  obs_local = sk$cases,
                  obs_imported = rep(0, nrow(sk)),
                  cut = cut)

 model <- stanmodels$epimodel

 fit <- rstan::sampling(model,
                        iter = 2000,
                        data = dat_test,
                        chains = 4,
                        control = list(adapt_delta = 0.8))

 res <- rstan::extract(fit)

 res_full <- data.frame(meda = c(apply(res$inf_R,MARGIN = 2, median),rep(NA,dat_test$cut)),
                        UQa = c(apply(res$inf_R, MARGIN = 2,
                                      FUN = function(x){quantile(x, prob=0.975)}),rep(NA,dat_test$cut)),
                        LQa = c(apply(res$inf_R, MARGIN = 2,
                                      FUN = function(x){quantile(x, prob=0.025)}),rep(NA,dat_test$cut)),
                        med = apply(res$R,MARGIN = 2, median),
                        UQ = apply(res$R, MARGIN = 2,
                                   FUN = function(x){quantile(x, prob=0.975)}),
                        LQ = apply(res$R, MARGIN = 2,
                                   FUN = function(x){quantile(x, prob=0.025)}),
                        date = seq.Date(from = min(sk$date), length.out = dat_test$t, by = "day"),
                        confirm = c(apply(res$inf_cases,MARGIN = 2, median),rep(NA,dat_test$cut))
 )

 min_date <- sk$date[which(sk$cases > 10)][1]

 p2 <- res_full %>%
  dplyr::filter(date >= min_date) %>%
  ggplot2::ggplot(ggplot2::aes(x=date)) +
  ggplot2::geom_line(ggplot2::aes(y = meda), col = "blue") +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = LQa, ymax = UQa), alpha = 0.1, fill= "blue") +
  cowplot::theme_cowplot() +
  ggplot2::geom_hline(yintercept = 1, lty = 2) +
  ggplot2::ylab("Reproduction number") +
  ggplot2::xlab("") +
  ggplot2::scale_x_date(date_labels = "%b %d", date_breaks = "2 weeks",
                        limits = c(min_date,max(sk$date)))


 p1 <- data.frame(date = sk$date, confirm = sk$cases) %>%
  dplyr::filter(date >= min_date) %>%
  ggplot2::ggplot(ggplot2::aes(x = date, y = confirm)) +
  ggplot2::geom_bar(stat="identity") +
  ggplot2::geom_bar(data = res_full,
                    stat = "identity",
                    fill = "red2",
                    alpha = 0.5) +
  ggplot2::scale_x_date(date_labels = "%b %d", date_breaks = "2 weeks",
                        limits = c(min_date,max(sk$date))) +
  cowplot::theme_cowplot() +
  ggplot2::ylab("Daily cases by infection/confirmation date")


 return(list("cases_plot" = p1, "Rt_plot" = p2, "fit" = fit, "obs" = sk))
}
