get_model_info_hanam <- function(study) {
    antigenic_map <- get_antigenic_map()
    par_tab <- get_par_tab()

    model_cross <- make_model_info_hanam(
            study = study,
            par_tab = par_tab,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "novac",
            model_name = "cross",
            pars_plot = c("mu", "sigma1", "tau"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_normal,
            custom_antigenic_maps_func = make_antigenic_maps_default
    )


     par_tab_full <-  rbind(par_tab,
        list("mu_vac",       1, 1, 0.1, 0, 1,  0, 1, 1),
        list("mu_short_vac", 0, 1, 0.1, 0, 4,  1, 3, 1),
        list("wane_vac",     0, 1, 0.1, 0, 1,    0.01,  0.1, 1),
        list("tau_vac",      0, 1, 0.1,  0, 1,  0.01, 0.1, 0),
        list("sigma1_vac",   1, 1, 0.1, 0.1, 10,  0.5, 2.0, 1),
        list("sigma2_vac",   0, 1, 0.1, 0, 0.3,  0.01, 0.1, 1),
        list("rho_boost",   1, 1, 0.1, 0, 10,  0.5, 2, 1),
        list("rho_wane",   1, 1, 0.1, 0, 10,  0.5, 2, 1)
        )
    
    par_tab_full <- par_tab_full %>% # mu, tau, sigma1 also done
        mutate(fixed = replace(fixed, names == "mu", 0)) %>%
        mutate(fixed = replace(fixed, names == "tau", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma1", 0)) %>%
        mutate(fixed = replace(fixed, names == "mu_vac", 0)) %>%
        mutate(fixed = replace(fixed, names == "mu_short_vac", 0)) %>%
        mutate(fixed = replace(fixed, names == "wane_vac", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma1_vac", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma2_vac", 0))

    par_tab_full$steps <- par_tab_full$steps / 10

    model_full <- make_model_info_hanam(
            study = study,
            par_tab = par_tab_full,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "full",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "mu_short_vac", "wane_vac", "sigma1_vac", "sigma2_vac"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_vac_general,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
    )


    all_models <- list(model_cross, model_full)
    all_models_hanam_cross <- all_models %>% map(~convert_to_resolution(.x, 52))

    save(all_models_hanam_cross, file = here::here("data", paste0("modelinfo_cross_", study$study_name_short, ".RDS")))
}

logit <- function(x, b) {
    b / (1 + exp(-x))
}

logit.inv <- function(y) {
    log((1 - y) / y)
}