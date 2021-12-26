setup_run_serosolver <- function(
    model_info,
    chains,
    iterations,
    pt = FALSE,
    filename,
    cross_sectional = TRUE,
    continue_run = FALSE
  ) {
    #setup parallel backend to use many processors "longterm, shortterm"
    #cores = detectCores()
    cl <- makeCluster(4, outfile = "") #not to overload your computer
    registerDoParallel(cl)
    # get output files
    no_chains <- chains

    if (cross_sectional) {
      output_file <- "cross_sec"
      model_info$create_post$titre_data <- model_info$create_post$titre_data %>% filter(sample_time == 0)
    } else {
      output_file <- "full"
    }

    dir.create(here::here("outputs",  model_info$others$study$study_name_short, "fits", filename, output_file, model_info$others$model_name))
    filename <- here::here("outputs",  model_info$others$study$study_name_short, "fits", filename, output_file, model_info$others$model_name)
    filenames <- paste0(filename, "/post", "_", 1:chains)
    # Create the posterior solving function that will be used in the MCMC framework 

    model_func <- create_posterior_func(
      par_tab = model_info$create_post$par_tab,
      vaccination_histories = model_info$create_post$vac_history,
      titre_dat = model_info$create_post$titre_data,
      antigenic_map = model_info$create_post$antigenic_map,
      version = model_info$create_post$prior_version,
      custom_ab_kin_func = model_info$custom_ab_kin_func,
      custom_antigenic_maps_func = model_info$custom_antigenic_maps_func)
      
    temp_ladder_len <- 10
    mcmc_sampling_pars <- list("iterations" = iterations,"popt "= 0.44,"popt_hist"=0.44,
                          "opt_freq"=1000, "thin" = 1000, "adaptive_period" = iterations / 10,
                          "save_block"=1000, "thin_hist" = 1000, "hist_sample_prob" = 1,
                          "switch_sample"=2, "burnin" = iterations / 5, "inf_propn"=0.5,
                          "move_size"=3,"hist_opt"=1,"swap_propn"=0.5,
                          "hist_switch_prob"=0.5,"year_swap_propn" = 1,
                          "temperature" = seq(1, temp_ladder_len, by = 1),
                          "parallel_tempering_iter" = 5, "min_temp_diff" = 0.1)

                                    # function in post
    ## Not all random starting conditions return finite likelihood, so for each chain generate random
    ## conditions until we get one with a finite likelihood#
   # for (f in filenames) {
  # for (f in filenames) {
    res <- foreach(f = filenames, .packages = c('rcppfunchcw', 'magrittr', 'serosolver', 'tidyr', 'data.table', 'plyr', 'dplyr')) %dopar% {

        prior_on_mu_vac <- function(par_tab) {
            model_info$prior_function
        }
        start_prob <- -Inf
        while(!is.finite(start_prob)) {
            ## Generating starting antibody kinetics parameters
            start_tab <- generate_start_tab(model_info$create_post$par_tab)
            ## Generate starting infection history
            start_inf <- setup_infection_histories_titre(model_info$create_post$titre_data,
                                                        model_info$create_post$strain_isolation_times,
                                                        space = 3, titre_cutoff = 5) # 50 x XX
            start_prob <- sum(model_func(start_tab$values, start_inf)[[1]])
        }
        if (pt) {
          start_tab_list <- rep(list(start_tab), temp_ladder_len)
          start_inf_list <- rep(list(start_inf), temp_ladder_len)

          res <- run_MCMC_pt(
                  par_tab = start_tab_list,
                  titre_dat = model_info$create_post$titre_data,
                  vaccination_histories = model_info$create_post$vac_history,
                  antigenic_map = model_info$create_post$antigenic_map,
                  start_inf_hist = start_inf_list,
                  mcmc_pars = mcmc_sampling_pars,
                  filename = f,
                  CREATE_POSTERIOR_FUNC = create_posterior_func,
                  CREATE_PRIOR_FUNC = prior_on_mu_vac,
                  continue_run = continue_run,
                  custom_ab_kin_func = model_info$custom_ab_kin_func,
                  custom_antigenic_maps_func = model_info$custom_antigenic_maps_func
                  )
        } else {
          res <- run_MCMC(par_tab = start_tab,
                                  titre_dat = model_info$create_post$titre_data,
                                  vaccination_histories = model_info$create_post$vac_history,
                                  vaccination_histories_mat = model_info$create_post$vaccination_histories_mat,
                                  antigenic_map = model_info$create_post$antigenic_map,
                                  start_inf_hist = start_inf,
                                  mcmc_pars = list("iterations"=2000,"popt"=0.44,"popt_hist"=0.44,
                                          "opt_freq"=100,"thin"=10, "adaptive_period" = 200, 
                                          "save_block"=100, "thin_hist"=100,"hist_sample_prob"=1,
                                          "switch_sample"=2, "burnin"=0, "inf_propn"=0.5, 
                                          "move_size"=3,"hist_opt"=1,"swap_propn"=0.5,
                                          "hist_switch_prob"=0.5,"year_swap_propn" = 1),
                                  filename = f,
                                  CREATE_POSTERIOR_FUNC = create_posterior_func,
                                  CREATE_PRIOR_FUNC = prior_on_mu_vac,
                                  version = model_info$create_post$prior_version,
                                  continue_run = continue_run)
        }
  }

}
