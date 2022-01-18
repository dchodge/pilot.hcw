convert_to_ss_vac <- function(short_name, sero_data, vac_data, sample_yr) {
    vac_df <- sero_data %>% select(pid) %>% unique %>% left_join(vac_data, by = "pid") %>% select(!comorbidites)
    vac_df$individual <- factor(vac_df$pid, levels = unique(vac_df$pid)) %>% as.integer
    vac_df$virus <- vac_df$vac_yr
    vac_df$time <- vac_df$vac_yr
    vac_df <- vac_df %>% mutate(vac_flag =
        case_when(response == "No" ~ 0, response == "Don't know" ~ 0, response == "Yes" ~ 1)) %>%
        select(individual, vac_flag, virus, time)

    vac_df$virus <- as.numeric(as.character(vac_df$virus))
    vac_df$time <- as.numeric(as.character(vac_df$time))

    vac_df <- vac_df %>% complete(individual, virus = min(vac_df$virus):2018,
            time = min(vac_df$virus):2018, fill = list(vac_flag = 0)) %>%
            filter(virus == time) %>% as.data.frame %>% filter(time <= sample_yr)
    vac_df
}

convert_to_ss_sero_hcwpre <- function(shortname, sero_data, part_data, sample_yr) {

    data_b4_vac <- sero_data
    data_b4_vac <- data_b4_vac %>% left_join(unique(part_data), by = "pid")
    data_b4_vac$individual <- factor(data_b4_vac$pid, levels = unique(data_b4_vac$pid)) %>% as.integer

    data_b4_vac <- data_b4_vac %>% mutate(samples =
        case_when(
            sample_time == 0~sample_yr,
            sample_time == 30~sample_yr + 1 / 12,
            sample_time == 60~sample_yr + 2 / 12,
            sample_time == 180~sample_yr + 6 / 12
        ))

    data_b4_vac$virus <- data_b4_vac$virus_isol_yr
    data_b4_vac$titre <- as.integer(round(data_b4_vac$titre_value, 0))
    data_b4_vac$run <- as.integer(1)
    data_b4_vac$DOB <- as.integer(data_b4_vac$yob)
    data_b4_vac$group <- as.integer(1)
    data_b4_vac %>% select(individual, DOB, virus, titre, samples, run, group, sample_time) %>%
        group_by(individual, DOB, virus, samples, run, group, sample_time) %>%
        dplyr::summarise(titre = as.integer(mean(titre))) %>%
        ungroup %>% na.omit %>% as.data.frame
}

convert_to_ss_sero_hanam <- function(shortname, sero_data, part_data, sample_yr) {

    data_b4_vac <- sero_data
    part_data_alt <- unique(part_data) %>% select(pid, yob)
    data_b4_vac <- data_b4_vac %>% left_join(unique(part_data_alt), by = "pid")
    data_b4_vac$individual <- factor(data_b4_vac$pid, levels = unique(data_b4_vac$pid)) %>% as.integer

    data_b4_vac <- data_b4_vac %>% mutate(samples =
        case_when(
            sample_time == 0~sample_yr,
            sample_time == 4~sample_yr + 0.5 / 52,
            sample_time == 7~sample_yr + 1 / 52,
            sample_time == 14~sample_yr + 2 / 52,
            sample_time == 21~sample_yr + 3 / 52,
            sample_time == 280~sample_yr + 40 / 52
        ))

    data_b4_vac$virus <- data_b4_vac$virus_isol_yr
    data_b4_vac$titre <- as.integer(round(data_b4_vac$titre_value, 0))
    data_b4_vac$run <- as.integer(1)
    data_b4_vac$DOB <- as.integer(data_b4_vac$yob)
    data_b4_vac$group <- as.integer(1)
    data_b4_vac %>% select(individual, DOB, virus, titre, samples, run, group, sample_time) %>%
        group_by(individual, DOB, virus, samples, run, group, sample_time) %>%
        dplyr::summarise(titre = as.integer(mean(titre))) %>%
        ungroup %>% na.omit %>% as.data.frame
    data_b4_vac$titre <- data_b4_vac$titre - min(data_b4_vac$titre)
    data_b4_vac %>% filter(sample_time %in% c(0, 7, 14, 21, 280))
    data_b4_vac
}

get_antigenic_map <- function() {

    antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
    antigenic_coords <- read.csv(file = antigenic_coords_path, stringsAsFactors = FALSE)
    antigenic_map <- generate_antigenic_map(antigenic_coords, 1)
    antigenic_map <- bind_rows(antigenic_map, data.frame(
        x_coord = c(43.3, 44.3, 45.3, 46.3, 47.3),
        y_coord = c(14.3, 15.3, 16.3, 17.3, 18.3),
        inf_times = c(2016, 2017, 2018, 2019, 2020)
    ))
    antigenic_map
}

get_par_tab <- function() {

    par_tab_path <- system.file("extdata", "par_tab_base.csv", package = "serosolver")
    par_tab <- read.csv(file = par_tab_path, stringsAsFactors = FALSE)

    ## Set parameters for Beta prior on infection histories
    beta_pars <- find_beta_prior_mode(0.15, 4)
    par_tab[par_tab$names == "alpha", "values"] <- beta_pars$alpha
    par_tab[par_tab$names == "beta", "values"] <- beta_pars$beta

    ## Remove phi parameters, as these are integrated out under prior version 2
    par_tab <- par_tab[par_tab$names != "phi", ]

    ## Fix all short term parameters to 0
    par_tab[par_tab$names %in% c("mu_short", "sigma2", "wane"), "fixed"] <- 1 # mu_short, waning and sigma2 are fixed
    par_tab[par_tab$names %in% c("mu_short", "sigma2", "wane"), "values"] <- 0 # set these values to 0

    par_tab
}

make_model_info_hcwpre <- function(prior_version, output_file, antigenic_map, par_tab, study, sample_yr, vacc_type, model_name, pars_plot,
        prior_function, custom_ab_kin_func, custom_antigenic_maps_func) { 
    
    sero_data <- study$sero_data
    part_data <- study$part_data
    vac_data <- study$vac_data

    titre_data <- convert_to_ss_sero_hcwpre(study$study_name_short, sero_data, part_data, sample_yr) %>% filter(virus <= sample_yr) %>% as.data.frame

    if (vacc_type == "novac") {
        vac_history <- NULL
    } else {
        vac_history <- convert_to_ss_vac(study$study_name_short, sero_data, vac_data, sample_yr) %>% data.frame
        vac_history <- vac_history %>% mutate(vac_flag = case_when(time == sample_yr~1, time != sample_yr~vac_flag))
        vac_history <- vac_history %>% group_by(individual) %>% mutate(prev_vac = cumsum(vac_flag)) %>% as.data.frame
    }

    antigenic_map <- antigenic_map %>% filter(inf_times <= sample_yr)
    antigenic_map <- antigenic_map[antigenic_map$inf_times >= min(titre_data$virus) & antigenic_map$inf_times <= max(titre_data$virus), ]
    strain_isolation_times <- unique(antigenic_map$inf_times)
    ## Maximum recordable log titre in these data is 8
    par_tab[par_tab$names == "MAX_TITRE", "values"] <- max(titre_data$titre)
    
    vaccination_histories_mat <- make_vac_matrix(vac_history, titre_data, strain_isolation_times)

    create_post <- list(par_tab = par_tab,
        vac_history = vac_history,
        vaccination_histories_mat = vaccination_histories_mat,
        titre_data = titre_data,
        antigenic_map = antigenic_map,
        strain_isolation_times = strain_isolation_times,
        prior_version = prior_version
    )

    others <- list(sample_yr = sample_yr, study = study, output_file = output_file, model_name = model_name,
        pars_plot = pars_plot)
    
    list(create_post = create_post, others = others, prior_function = prior_function,
        custom_ab_kin_func = custom_ab_kin_func, custom_antigenic_maps_func = custom_antigenic_maps_func)
}


make_model_info_hanam <- function(prior_version, output_file, antigenic_map, par_tab, study, sample_yr, vacc_type, model_name, pars_plot,
        prior_function, custom_ab_kin_func, custom_antigenic_maps_func) { 
    
    sero_data <- study$sero_data
    part_data <- study$part_data
    vac_data <- study$vac_data

    titre_data <- convert_to_ss_sero_hanam(study$study_name_short, sero_data, part_data, sample_yr) %>% filter(virus <= sample_yr) %>% as.data.frame

    if (vacc_type == "novac") {
        vac_history <- NULL
    } else {
        sit <- unique(antigenic_map$inf_times)
        vac_history_pre <- data.frame(
            individual = 1:length(unique(sero_data$pid)),
            virus = sample_yr,
            time = sample_yr,
            vac_flag = 1
        )

        vac_history <- vac_history_pre  %>% complete(individual, nesting(virus = sit, time = sit), fill = list(vac_flag = 0)) %>% as.data.frame
        vac_history <- vac_history %>% group_by(individual) %>% mutate(prev_vac = cumsum(vac_flag)) %>% filter(time <= sample_yr)
    }

    antigenic_map <- antigenic_map %>% filter(inf_times <= sample_yr)
    antigenic_map <- antigenic_map[antigenic_map$inf_times >= min(titre_data$virus) & antigenic_map$inf_times <= max(titre_data$virus), ]
    strain_isolation_times <- unique(antigenic_map$inf_times)
    ## Maximum recordable log titre in these data is 8
    par_tab[par_tab$names == "MAX_TITRE", "values"] <- max(titre_data$titre)
    
    vaccination_histories_mat <- make_vac_matrix(vac_history, titre_data, strain_isolation_times)

    create_post <- list(par_tab = par_tab,
        vac_history = vac_history,
        vaccination_histories_mat = vaccination_histories_mat,
        titre_data = titre_data,
        antigenic_map = antigenic_map,
        strain_isolation_times = strain_isolation_times,
        prior_version = prior_version
    )

    others <- list(sample_yr = sample_yr, study = study, output_file = output_file, model_name = model_name,
        pars_plot = pars_plot)
    
    list(create_post = create_post, others = others, prior_function = prior_function,
        custom_ab_kin_func = custom_ab_kin_func, custom_antigenic_maps_func = custom_antigenic_maps_func)
}

convert_to_resolution <- function(model_info, resolution) {
  if(!is.null(model_info$create_post$vac_history)) {
    model_info$create_post$vac_history <- model_info$create_post$vac_history %>%
      mutate(virus = virus * resolution, time = time * resolution)
  }
  model_info$create_post$titre_data <- model_info$create_post$titre_data %>%
    mutate(DOB = DOB * resolution, virus = virus * resolution, samples = samples * resolution)
  model_info$create_post$antigenic_map <- model_info$create_post$antigenic_map %>%
    mutate(inf_times = inf_times * resolution)
  model_info$create_post$strain_isolation_times <- model_info$create_post$strain_isolation_times * resolution

  model_info
}

make_vac_matrix <- function(vac_history, titre_data, strain_isolation_times) {
  if (sum(vac_history$vac_flag) == 0) {
    vac_mat <- matrix(0, titre_data$ind %>% unique %>% length, titre_data$ind %>% unique %>% length, strain_isolation_times %>% length)
  } else {
    vac_mat <- vac_history %>%
            select(vac_flag, virus, individual) %>%
            pivot_wider(values_from = vac_flag, names_from = virus) %>%
            select(!individual) %>%
            as.matrix
  }
  vac_mat
}