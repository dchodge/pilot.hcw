get_model_info_hcw_pre_full <- function(study) {
    antigenic_map <- get_antigenic_map()
    par_tab <- get_par_tab()
    par_tab_vac <-  rbind(par_tab,
        list("mu_vac",       0, 1, 0.1, 0, 1,  0, 1, 1),
        list("mu_short_vac", 0, 1, 0.1, 0, 8,  1, 3, 1),
        list("wane_vac",     0, 1, 0.1, 0, 1,    0.01,  0.1, 1),
        list("tau_vac",      0, 1, 0.1,  0, 1,  0.01, 0.1, 0),
        list("sigma1_vac",   0, 1, 0.1, 0.1, 10,  0.5, 2.0, 1),
        list("sigma2_vac",   0, 1, 0.1, 0, 1,  0.01, 0.1, 1),
        list("rho_boost",   1, 1, 0.1, 0, 10,  0.5, 2, 1),
        list("rho_wane",   1, 1, 0.1, 0, 10,  0.5, 2, 1)
        )

    par_tab_vac <- par_tab_vac %>%
        mutate(values = replace(values, names == "mu", 2.126)) %>%
        mutate(values = replace(values, names == "tau", 0.01610)) %>%
        mutate(values = replace(values, names == "sigma1", 0.0787)) %>%
        mutate(values = replace(values, names == "mu_vac", 0.06463)) %>%
        mutate(values = replace(values, names == "sigma1_vac", 0.33883))
    par_tab_vac <- par_tab_vac %>%
        mutate(lower_start = replace(lower_start, names == "mu", 1.800965)) %>%
        mutate(lower_start = replace(lower_start, names == "tau", 0.0002655416)) %>%
        mutate(lower_start = replace(lower_start, names == "sigma1", 0.07457849)) %>%
        mutate(lower_start = replace(lower_start, names == "mu_vac", 0.03930816)) %>%
        mutate(lower_start = replace(lower_start, names == "sigma1_vac", 0.1055997))
    par_tab_vac <- par_tab_vac %>%
        mutate(upper_start = replace(upper_start, names == "mu", 2.420486)) %>%
        mutate(upper_start = replace(upper_start, names == "tau", 0.0450771230)) %>%
        mutate(upper_start = replace(upper_start, names == "sigma1", 0.08247699)) %>%
        mutate(upper_start = replace(upper_start, names == "mu_vac", 0.11783472)) %>%
        mutate(upper_start = replace(upper_start, names == "sigma1_vac", 1.0052480))
    par_tab_vac <- par_tab_vac %>%
        mutate(steps = replace(steps, names == "mu", 0.1)) %>%
        mutate(steps = replace(steps, names == "tau", 0.004)) %>%
        mutate(steps = replace(steps, names == "sigma1", 0.01)) %>%
        mutate(steps = replace(steps, names == "mu_vac", 0.01)) %>%
        mutate(steps = replace(steps, names == "sigma1_vac", 0.01))
    par_tab_vac$steps <- par_tab_vac$steps / 10

    model_vac_full <- make_model_info(
            study = study,
            par_tab = par_tab_vac,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "base",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "sigma1_vac", "mu_short_vac", "wane_vac", "sigma2_vac"),
            prior_function = function(cur_pars) { 
                require(triangle)
                mu <- cur_pars[["mu"]]
                tau <- cur_pars[["tau"]]
                sigma1 <- cur_pars[["sigma1"]]
                mu_vac <- cur_pars[["mu_vac"]]
                sigma1_vac <- cur_pars[["sigma1_vac"]]
                log_p <- log(dtriangle(mu, 1.800965, 2.420486, 2.126))
                log_p <- log_p + log(dtriangle(tau, 0.0002655416, 0.0450771230, 0.01610))
                log_p <- log_p + log(dtriangle(sigma1, 0.07457849, 0.08247699, 0.0787))
                log_p <- log_p + log(dtriangle(mu_vac, 0.03930816, 0.11783472, 0.06463))
                log_p <- log_p + log(dtriangle(sigma1_vac, 0.1055997, 1.0052480, 0.33883))
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
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "sigma1_vac",  "mu_short_vac", "wane_vac", "sigma2_vac", "rho_boost"),
            prior_function = function(cur_pars) { 
                require(triangle)
                mu <- cur_pars[["mu"]]
                mu_vac <- cur_pars[["mu_vac"]]
                tau <- cur_pars[["tau"]]
                sigma1 <- cur_pars[["sigma1"]]
                sigma1_vac <- cur_pars[["sigma1_vac"]]
                log_p <- log(dtriangle(mu, 1.800965, 2.420486, 2.126))
                log_p <- log_p + log(dtriangle(tau, 0.0002655416, 0.0450771230, 0.01610))
                log_p <- log_p + log(dtriangle(sigma1, 0.07457849, 0.08247699, 0.0787))
                log_p <- log_p + log(dtriangle(mu_vac, 0.03930816, .11783472, 0.06463))
                log_p <- log_p + log(dtriangle(sigma1_vac, 0.1055997, 1.0052480, 0.33883))
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
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "sigma1_vac", "mu_short_vac", "wane_vac", "sigma2_vac", "rho_wane"),
            prior_function = function(cur_pars) { 
                require(triangle)
                 mu <- cur_pars[["mu"]]
                mu_vac <- cur_pars[["mu_vac"]]
                tau <- cur_pars[["tau"]]
                sigma1 <- cur_pars[["sigma1"]]
                sigma1_vac <- cur_pars[["sigma1_vac"]]
                log_p <- log(dtriangle(mu, 1.800965, 2.420486, 2.126))
                log_p <- log_p + log(dtriangle(tau, 0.0002655416, 0.0450771230, 0.01610))
                log_p <- log_p + log(dtriangle(sigma1, 0.07457849, 0.08247699, 0.0787))
                log_p <- log_p + log(dtriangle(mu_vac, 0.03930816, .11783472, 0.06463))
                log_p <- log_p + log(dtriangle(sigma1_vac, 0.1055997, 1.0052480, 0.33883))
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
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "sigma1_vac", "sigma1_vac", "mu_short_vac", "wane_vac", "sigma2_vac", "rho_boost", "rho_wane"),
            prior_function = function(cur_pars) { 
                require(triangle)
                mu <- cur_pars[["mu"]]
                mu_vac <- cur_pars[["mu_vac"]]
                tau <- cur_pars[["tau"]]
                sigma1 <- cur_pars[["sigma1"]]
                sigma1_vac <- cur_pars[["sigma1_vac"]]
                log_p <- log(dtriangle(mu, 1.800965, 2.420486, 2.126))
                log_p <- log_p + log(dtriangle(tau, 0.0002655416, 0.0450771230, 0.01610))
                log_p <- log_p + log(dtriangle(sigma1, 0.07457849, 0.08247699, 0.0787))
                log_p <- log_p + log(dtriangle(mu_vac, 0.03930816, .11783472, 0.06463))
                log_p <- log_p + log(dtriangle(sigma1_vac, 0.1055997, 1.0052480, 0.33883))
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