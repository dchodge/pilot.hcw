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

    par_tab_LN <- rbind(par_tab,
        list("mu_LN",        5, 1, 0.1,   0, 100,  5, 10, 1),
        list("sigma_LN",     5, 1, 0.1,   0, 100,  2, 4, 1),
        list("alpha_LN",      100, 1, 0.1, 0, 1000,  100, 500, 1),
        list("sigma1_vac",   1, 1, 0.1,   0.1, 10,  0.5, 2.0, 1),
        list("sigma2_vac",   0, 1, 0.1,    0,  0.3,  0.01, 0.1, 1)
    )
    par_tab_LN <- par_tab_LN %>% # mu, tau, sigma1 also done
        mutate(fixed = replace(fixed, names == "mu", 0)) %>%
        mutate(fixed = replace(fixed, names == "tau", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma1", 0)) %>%
        mutate(fixed = replace(fixed, names == "mu_LN", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma_LN", 0)) %>%
        mutate(fixed = replace(fixed, names == "alpha_LN", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma1_vac", 0)) 

    model_full_LN <- make_model_info_hanam(
            study = study,
            par_tab = par_tab_LN,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "LN",
            pars_plot = c("mu", "sigma1", "tau", "mu_LN", "sigma_LN", "alpha_LN", "sigma1_vac"),
            prior_function = function(cur_pars) { 
                return(0)
            },
            custom_ab_kin_func = ab_kin_vac_log_normal,
            custom_antigenic_maps_func = make_antigenic_maps_vac2
    )

    par_tab_SP <- rbind(par_tab,
        list("th1_SP",         0.5, 1, 0.01,   0, 1,  0, 1, 1),
        list("th2_SP",         0.5, 1, 0.01,   0, 1,  0, 1, 1),
        list("th3_SP",         0.5, 1, 0.01,   0, 1,  0, 1, 1),
        list("alpha_SP",      10, 1,  0.1,   0, 1000,  0, 10, 1),
        list("sigma1_vac",     1, 1, 0.1,   0.1, 10,  0.5, 2.0, 1),
        list("sigma2_vac",     0, 1, 0.1,    0,  0.3,  0.01, 0.1, 1)
    )

    par_tab_SP <- par_tab_SP %>% # mu, tau, sigma1 also done
        mutate(fixed = replace(fixed, names == "mu", 0)) %>%
        mutate(fixed = replace(fixed, names == "tau", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma1", 0)) %>%
        mutate(fixed = replace(fixed, names == "th1_SP", 0)) %>%
        mutate(fixed = replace(fixed, names == "th2_SP", 0)) %>%
        mutate(fixed = replace(fixed, names == "th3_SP", 0)) %>%
        mutate(fixed = replace(fixed, names == "alpha_SP", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma1_vac", 0))

    model_full_SP <- make_model_info_hanam(
            study = study,
            par_tab = par_tab_SP,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "SP",
            pars_plot = c("mu", "sigma1", "tau", "th1_SP", "th2_SP", "th3_SP", "alpha_SP", "sigma1_vac"),
            prior_function = function(cur_pars) { 
                return(0)
            },
            custom_ab_kin_func = ab_kin_vac_set_point,
            custom_antigenic_maps_func = make_antigenic_maps_vac2
    )


    all_models <- list(model_cross, model_full_LN, model_full_SP)
    all_models_hanam <- all_models %>% map(~convert_to_resolution(.x, 52))

    save(all_models_hanam, file = here::here("data", paste0("modelinfo_cross_", study$study_name_short, ".RDS")))
}

logit <- function(x, b) {
    b / (1 + exp(-x))
}

logit.inv <- function(y) {
    log((1 - y) / y)
}