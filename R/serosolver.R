#  hcw_ss <- convert_to_ss_df(hcw_sero_data, hcw_part_data, 2020)
#  hcw_pre_ss <- convert_to_ss_df(hcw_sero_data, hcw_part_data, 2020)

setup_run_serosolver <- function(model_info, study, sample_yr, chains) {
    #setup parallel backend to use many processors
    cores = detectCores()
    cl <- makeCluster(3) #not to overload your computer
    registerDoParallel(cl)

    # study data
    sero_data <- study$sero_data
    part_data <- study$part_data
    vac_data <- study$part_data
    titre_data <- convert_to_ss_df(sero_data, part_data, sample_yr)
    vac_history <- convert_to_ss_vac(sero_data, vac_data)
    # model data
    antigenic_map <- model_info$antigenic_map
    par_tab <- model_info$par_tab
    prior_version <- model_info$prior_version
    output_file <- model_info$output_file
    # trim antigenic map according to titre data
    antigenic_map[antigenic_map$inf_times >= min(titre_data$virus) & antigenic_map$inf_times <= max(titre_data$virus), ]
    strain_isolation_times <- unique(antigenic_map$inf_times)
    # get output files
    no_chains <- chains
    filename <- paste0(output_file, study$study_name_short, "/", "original", "/")
    filenames <- paste0(filename, "original", "_", 1:chains)
    # Create the posterior solving function that will be used in the MCMC framework 
    model_func <- create_posterior_func(par_tab = par_tab,
                                    titre_dat = titre_data,
                                    antigenic_map = antigenic_map,
                                    version = prior_version) # function in post
    ## Not all random starting conditions return finite likelihood, so for each chain generate random
    ## conditions until we get one with a finite likelihood

    res <- foreach(f = filenames, .packages = c('serosolver','data.table','plyr')) %dopar% {
        start_prob <- -Inf
        while(!is.finite(start_prob)) {
            ## Generating starting antibody kinetics parameters
            start_tab <- generate_start_tab(par_tab)
            ## Generate starting infection history
            start_inf <- setup_infection_histories_titre(titre_data, strain_isolation_times,
                                                        space = 3, titre_cutoff = 4)
            start_prob <- sum(model_func(start_tab$values, start_inf)[[1]])
        }

        res <- run_MCMC(par_tab = start_tab,
                        titre_dat = as.data.frame(titre_data),
                        antigenic_map = antigenic_map,
                        start_inf_hist = start_inf,
                        mcmc_pars = c("iterations"=200000,"popt"=0.44,"popt_hist"=0.44,
                                "opt_freq"=1000,"thin"=1,"adaptive_period"=50000, 
                                "save_block"=1000, "thin_hist"=100,"hist_sample_prob"=1,
                                "switch_sample"=2, "burnin"=0, "inf_propn"=0.5, 
                                "move_size"=3,"hist_opt"=1,"swap_propn"=0.5,
                                "hist_switch_prob"=0.5,"year_swap_propn"=1),
                        filename = f,
                        CREATE_POSTERIOR_FUNC = create_posterior_func,
                        version = prior_version)
  }
}