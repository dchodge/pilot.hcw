get_model_info_hcw_pre_full <- function(study) {
    antigenic_map <- get_antigenic_map()
    par_tab <- get_par_tab()
     par_tab_full <-  rbind(par_tab,
        list("mu_vac",       1, 1, 0.1, 0, 1,  0, 1, 1),
        list("mu_short_vac", 0, 1, 0.1, 0, 4,  1, 3, 1),
        list("wane_vac",     0, 1, 0.1, 0, 1,  0.5 / 12, 4 / 12, 1),
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

 #   par_tab_vac <- par_tab_vac %>%
  #      mutate(values = replace(values, names == "mu", 2.117883)) %>%
   #     mutate(values = replace(values, names == "tau", 0.01017651)) %>%
 #       mutate(values = replace(values, names == "sigma1", 0.07627891)) %>%
#        mutate(values = replace(values, names == "mu_vac", 0.1246219))
#    par_tab_vac <- par_tab_vac %>%
   #     mutate(lower_start = replace(lower_start, names == "mu", 1.891569)) %>%
     #   mutate(lower_start = replace(lower_start, names == "tau", 0.0002803596)) %>%
     #   mutate(lower_start = replace(lower_start, names == "sigma1", 0.07226025)) %>%
    #    mutate(lower_start = replace(lower_start, names == "mu_vac", 0.09169635)) 
 #   par_tab_vac <- par_tab_vac %>%
 #       mutate(upper_start = replace(upper_start, names == "mu", 2.374630)) %>%
  #      mutate(upper_start = replace(upper_start, names == "tau", 0.0284054135)) %>%
 #       mutate(upper_start = replace(upper_start, names == "sigma1", 0.08012938 )) %>%
 #       mutate(upper_start = replace(upper_start, names == "mu_vac", 0.15619417)) 
  #  par_tab_vac <- par_tab_vac %>%
  #      mutate(steps = replace(steps, names == "mu", 0.1)) %>%
  #      mutate(steps = replace(steps, names == "tau", 0.004)) %>%
  #      mutate(steps = replace(steps, names == "sigma1", 0.01)) %>%
  #      mutate(steps = replace(steps, names == "mu_vac", 0.01))



    par_tab_full$steps <- par_tab_full$steps / 10

    model_vac_full <- make_model_info_hcwpre(
            study = study,
            par_tab = par_tab_full,
            antigenic_map = antigenic_map,
            sample_yr = 2016,
            prior_version = 2,
            output_file = "hpc/outputs/",
            vacc_type = "vac",
            model_name = "base",
            pars_plot = c("mu", "sigma1", "tau", "mu_vac", "mu_short_vac", "wane_vac", "sigma1_vac", "sigma2_vac"),
            prior_function = function(cur_pars) { 
                    log_p <- 0
                    mu_short_vac <- cur_pars[["mu_short_vac"]]
                    wane_vac <- cur_pars[["wane_vac"]]
                    log_p <- log(dnorm(mu_short_vac, 2, 0.5))
                    log_p <- log_p + log(dunif(wane_vac, 0.5 / 12, 4 / 12))
                return(log_p)
            },
            custom_ab_kin_func = ab_kin_vac_general,
            custom_antigenic_maps_func = make_antigenic_maps_vac1
    )

    all_models <- list(model_vac_full)
    all_models_hcw_pre_full <- all_models %>% map(~convert_to_resolution(.x, 12))

    save(all_models_hcw_pre_full, file = here::here("data", paste0("modelinfo_full_", study$study_name_short, ".RDS")))
}

logit <- function(x, b) {
    b / (1 + exp(-x))
}

logit.inv <- function(y) {
    log((1 - y) / y)
}