get_model_info_hcw_pre_full <- function(study) {
    antigenic_map <- get_antigenic_map()
    par_tab <- get_par_tab()
    par_tab_vac <-  rbind(par_tab,
            list("mu_vac",       0.04, 0, 0.1, 0, 1,  0, 1, 1),
            list("mu_short_vac", 0, 0, 0.1, 0, 1,  0, 1, 1),
            list("wane_vac",     0, 0, 0.1, 0, 0.1,    0.001,  0.01, 1),
            list("tau_prev_vac", 0, 1, 0.1,  0, 1,  0.01, 0.1, 1),
            list("sigma1_vac",   1, 1, 0.1, 0, 1,  0.01, 0.1, 1),
            list("sigma2_vac",   0, 0, 0.1, 0, 1,  0.01, 0.1, 1),
            list("rho_boost",   1, 1, 0.1, 0, 10,  0.5, 2, 1),
            list("rho_wane",   1, 1, 0.1, 0, 10,  0.5, 2, 1)
            )

    par_tab_vac <- par_tab_vac %>%
        mutate(value = replace(values, names == "mu", 2.125))
        mutate(value = replace(values, names == "tau", 0.01))
        mutate(value = replace(values, names == "sigma1", 0.079))

    model_vac_full <- make_model_info(
            study = study,
            par_tab = par_tab_vac,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "base",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "mu_short_vac", "wane_vac", "sigma2_vac"),
            prior_function = function(cur_pars) { 
                require(triangle)
                mu <- cur_pars[["mu"]]
                mu_vac <- cur_pars[["mu_vac"]]
                tau <- cur_pars[["tau"]]
                sigma1 <- cur_pars[["sigma1"]]
                log_p <- log(dtriangle(mu, 1.75, 2.50, 2.125))
                log_p <- log_p + log(dtriangle(mu_vac, 0.02, 0.06, 0.04))
                log_p <- log_p + log(dtriangle(tau, 0.0, 0.04, 0.01))
                log_p <- log_p + log(dtriangle(sigma1, 0.074, 0.084, 0.079))
                return(log_p)
            },
            custom_ab_kin_func = ab_kin_vac_prev_hist,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
    )

    par_tab_m1 <- par_tab_vac %>%
        mutate(fixed = replace(fixed, names == "rho_boost", 0))
    model_vac_m1 <- make_model_info(
            study = study,
            par_tab = par_tab_m1,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "m1",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac",  "mu_short_vac", "wane_vac", "sigma2_vac", "rho_boost"),
            prior_function = function(cur_pars) { 
                require(triangle)
                mu <- cur_pars[["mu"]]
                mu_vac <- cur_pars[["mu_vac"]]
                tau <- cur_pars[["tau"]]
                sigma1 <- cur_pars[["sigma1"]]
                log_p <- log(dtriangle(mu, 1.75, 2.50, 2.125))
                log_p <- log_p + log(dtriangle(mu_vac, 0.2, 0.6, 0.4))
                log_p <- log_p + log(dtriangle(tau, 0.0, 0.04, 0.01))
                log_p <- log_p + log(dtriangle(sigma1, 0.074, 0.084, 0.079))
                return(log_p)
            },
            custom_ab_kin_func = ab_kin_vac_prev_hist,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
    )

    par_tab_m2 <- par_tab_vac %>%
        mutate(fixed = replace(fixed, names == "rho_wane", 0))
    model_vac_m2 <- make_model_info(study = study,
            par_tab = par_tab_m2,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "m2",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "mu_short_vac", "wane_vac", "sigma2_vac", "rho_wane"),
            prior_function = function(cur_pars) { 
                require(triangle)
                mu <- cur_pars[["mu"]]
                mu_vac <- cur_pars[["mu_vac"]]
                tau <- cur_pars[["tau"]]
                sigma1 <- cur_pars[["sigma1"]]
                log_p <- log(dtriangle(mu, 1.75, 2.50, 2.125))
                log_p <- log_p + log(dtriangle(mu_vac, 0.2, 0.6, 0.4))
                log_p <- log_p + log(dtriangle(tau, 0.0, 0.04, 0.01))
                log_p <- log_p + log(dtriangle(sigma1, 0.074, 0.084, 0.079))
                return(log_p)
            },
            custom_ab_kin_func = ab_kin_vac_prev_hist,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
            )

    par_tab_m3 <- par_tab_vac %>%
        mutate(fixed = replace(fixed, names == "rho_boost", 0)) %>%
        mutate(fixed = replace(fixed, names == "rho_wane", 0))
    model_vac_m3 <- make_model_info(study = study,
            par_tab = par_tab_m3,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "m3",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "sigma1_vac", "mu_short_vac", "wane_vac", "sigma2_vac", "rho_boost", "rho_wane"),
            prior_function = function(cur_pars) { 
                require(triangle)
                mu <- cur_pars[["mu"]]
                mu_vac <- cur_pars[["mu_vac"]]
                tau <- cur_pars[["tau"]]
                sigma1 <- cur_pars[["sigma1"]]
                log_p <- log(dtriangle(mu, 1.75, 2.50, 2.125))
                log_p <- log_p + log(dtriangle(mu_vac, 0.2, 0.6, 0.4))
                log_p <- log_p + log(dtriangle(tau, 0.0, 0.04, 0.01))
                log_p <- log_p + log(dtriangle(sigma1, 0.074, 0.084, 0.079))
                return(log_p)
            },
            custom_ab_kin_func = ab_kin_vac_prev_hist,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
            )

    all_models <- list(model_vac_full, model_vac_m1, model_vac_m2, model_vac_m3)
    all_models_hcw_pre_full <- all_models %>% map(~convert_to_resolution(.x, 12))

    save(all_models_hcw_pre_full, file = here::here("data", paste0("modelinfo_full_", study$study_name_short, ".RDS")))
}

logit <- function(x, b) {
    b / (1 + exp(-x))
}

logit.inv <- function(y) {
    log((1 - y) / y)
}