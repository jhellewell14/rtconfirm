# rtconfirm
Estimate Time Varying Reproduction Number Directly from Confirmed Case Counts

This repository contains the stan code for running a model that estimates the time-varying reproduction number directly from confirmed case counts.

To get started, try running:
```
library(rstan)
library(patchwork)
# If you want to run in parallel run
# options(mc.cores = parallel::detectCores())
res <- rtconfirm::run_rtconfirm(country = "Poland", cut = 10)
# Output plot
res[[1]] / res[[2]]
```

Country names are matched to data provided by the European Centre for Disease Control, some countries with spaces in use underscores like `South_Korea`, `United_States_of_America`, or `United_Kingdom`.

To get a full list of countries for which there are data available, run `sort(unique(NCoVUtils::get_ecdc_cases()$country))`


I'm starting to put together slides of the theory behind the stan model here: https://docs.google.com/presentation/d/1wg_BDQkhYTJqAylzQVK_wCnw2v6OlFf-sUmdJOx1xvg/edit?usp=sharing
