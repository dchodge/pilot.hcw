get_model_info_hcw_pre_cross <- function(study) {
    antigenic_map <- get_antigenic_map()
    par_tab <- get_par_tab()
    par_tab_vac <-  rbind(par_tab,
            list("mu_vac",       1, 1, 0.1, 0, 1,  0, 1, 1),
            list("mu_short_vac", 0, 1, 0.1, 0, 8,  1, 3, 1),
            list("wane_vac",     0, 1, 0.1, 0, 1,    0.01,  0.1, 1),
            list("tau_prev_vac", 0, 1, 0.1,  0, 1,  0.01, 0.1, 1),
            list("sigma1_vac",   0, 1, 0.1, 0, 1,  0.01, 0.1, 1),
            list("sigma2_vac",   0, 1, 0.1, 0, 1,  0.01, 0.1, 1)
            )
    par_tab$steps <- par_tab$steps / 10
    model_novac <- make_model_info(
            study = study,
            par_tab = par_tab,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "novac",
            model_name = "novac",
            pars_plot = c("mu", "sigma1", "tau"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_normal,
            custom_antigenic_maps_func = make_antigenic_maps_default
    )

    par_tab_vac_ <- par_tab_vac %>%
        mutate(values = replace(values, names == "vac_flag", 1))
    model_vac_ <- make_model_info(
            study = study,
            par_tab = par_tab_vac_,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "vac_base",
            pars_plot = c("mu", "sigma1", "tau"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_vac,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
    )

    par_tab_vac_m <- par_tab_vac_ %>%
        mutate(fixed = replace(fixed, names == "mu_vac", 0))
    model_vac_m <- make_model_info(study = study,
            par_tab = par_tab_vac_m,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "vac_m",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_vac,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
            )

    par_tab_vac_t <- par_tab_vac_ %>%
        mutate(fixed = replace(fixed, names == "tau_prev_vac", 0))
    model_vac_t <- make_model_info(study = study,
            par_tab = par_tab_vac_t,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "vac_t",
            pars_plot = c("mu", "sigma1", "tau", "tau_prev_vac"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_vac,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
            )


    par_tab_vac_s <- par_tab_vac_ %>%
        mutate(fixed = replace(fixed, names == "sigma1_vac", 0))
    model_vac_s <- make_model_info(study = study,
            par_tab = par_tab_vac_s,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "vac_s",
            pars_plot = c("mu", "sigma1", "tau", "sigma1_vac"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_vac,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
            )

    par_tab_vac_mt <- par_tab_vac_ %>%
        mutate(fixed = replace(fixed, names == "mu_vac", 0)) %>%
        mutate(fixed = replace(fixed, names == "tau_prev_vac", 0))
    model_vac_mt <- make_model_info(study = study,
            par_tab = par_tab_vac_mt,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "vac_mt",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "tau_prev_vac"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_vac,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
    )

    par_tab_vac_ms <- par_tab_vac_ %>%
        mutate(fixed = replace(fixed, names == "mu_vac", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma1_vac", 0))
    model_vac_ms <- make_model_info(study = study,
            par_tab = par_tab_vac_ms,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "vac_ms",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "sigma1_vac"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_vac,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
    )

    par_tab_vac_ts <- par_tab_vac_ %>%
        mutate(fixed = replace(fixed, names == "tau_prev_vac", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma1_vac", 0))
    model_vac_ts <- make_model_info(study = study,
            par_tab = par_tab_vac_ts,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "vac_ts",
            pars_plot = c("mu", "sigma1", "tau", "tau_prev_vac", "sigma1_vac"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_vac,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
    )


    par_tab_vac_mts <- par_tab_vac_ %>%
        mutate(fixed = replace(fixed, names == "tau_prev_vac", 0)) %>%
        mutate(fixed = replace(fixed, names == "sigma1_vac", 0)) %>% 
        mutate(fixed = replace(fixed, names == "mu_vac", 0))
    model_vac_mts <- make_model_info(study = study,
            par_tab = par_tab_vac_mts,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "vac_mts",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "tau_prev_vac", "sigma1_vac"),
            prior_function = function(cur_pars) { return(0)},
            custom_ab_kin_func = ab_kin_vac,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
    )

    all_models <- list(model_novac, model_vac_, model_vac_m, model_vac_s, model_vac_t, model_vac_ms, model_vac_mt, model_vac_ts, model_vac_mts)
    all_models_hcw_pre <- all_models %>% map(~convert_to_resolution(.x, 12))

    save(all_models_hcw_pre, file = here::here("data", paste0("modelinfo_cross_", study$study_name_short, ".RDS")))
}

logit <- function(x, b) {
    b / (1 + exp(-x))
}

logit.inv <- function(y) {
    log((1 - y) / y)
}