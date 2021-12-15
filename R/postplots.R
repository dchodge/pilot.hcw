load_mcmc_post <- function(model_info, fit_place, folder_name_1, cross_sectional, burnin, thin) {
    if (cross_sectional) {
        folder_name_2 <- "cross_sec"
    } else {
        folder_name_2 <- "full"
    }
    if (fit_place == "hpc") {
        output_file <- here::here("hpc", "outputs", model_info$other$study$study_name_short, model_info$other$model_name)
    } else if (fit_place == "local") {
        output_file <- here::here("outputs", model_info$other$study$study_name_short, "fits", folder_name_1, folder_name_2, model_info$other$model_name)
    }
    posteriors <- load_mcmc_chains(location = output_file, thin = thin, burnin = burnin,
                               par_tab = model_info$create_post$par_tab, unfixed = FALSE, convert_mcmc = TRUE)
    posteriors
}

# get log_prop and make a matrix in which rows are for each sample and column for each observation
plot_traces <- function(postfull, model_info, file) {

    list_chains <- postfull$theta_list_chains
    j <- 1

    L <- nrow(list_chains[[1]])
    B <- round(L / 2)
    list_chains1 <- lapply(list_chains, 
        function(x) x[B:L, 
            c(model_info$other$pars_plot, "total_infections", "lnlike", "prior_prob")])
  #  for (t_par in model_info$others$pars_plot_trans) {
  #      list_chains1 <- seq_len(length(list_chains1)) %>% map(~(list_chains1[[.x]] %>%
  #          as.data.frame %>% mutate(!!t_par := logit(list_chains1[[.x]][, t_par], model_info$other$upper_boundary[j]))  %>% as.matrix))
   #         j <- j + 1
  #  }

    p_theta_trace <- mcmc_trace(list_chains1)
    savepath <- here::here("outputs", model_info$other$study$study_name_short, "fits", file, model_info$other$model_name)
    dir.create(savepath)
    ggsave(p_theta_trace, filename = paste0(savepath, "/", "traceplot.pdf"))
}

    # get log_prop and make a matrix in which rows are for each sample and column for each observation
plot_histogram <- function(postfull, model_info, file) {

    list_chains <- postfull$theta_list_chains
    j <- 1

    L <- nrow(list_chains[[1]])
    B <- round(L / 2)
    list_chains1 <- lapply(list_chains,
        function(x) x[B:L, 
            c(model_info$other$pars_plot, "total_infections", "lnlike", "prior_prob")])
 #   for (t_par in model_info$others$pars_plot_trans) {
 #       list_chains1 <- seq_len(length(list_chains1)) %>% map(~(list_chains1[[.x]] %>%
  #          as.data.frame %>% mutate(!!t_par := logit(list_chains1[[.x]][, t_par], model_info$other$upper_boundary[j]))  %>% as.matrix))
  #          j <- j + 1
  #  }

    p_theta_trace <- mcmc_dens_overlay(list_chains1)
    savepath <- here::here("outputs", model_info$other$study$study_name_short, "fits", file, model_info$other$model_name)
    dir.create(savepath)
    ggsave(p_theta_trace, filename = paste0(savepath, "/", "histplot.pdf"))
}

plot_post_compare <- function(postfull, model_info, file) {
    df_s <- list()
    for (s in 2:9) {
        j <- 1
     #   for (t_par in model_info[[s]]$others$pars_plot_trans) {
      #      list_chains1 <- seq_len(length(list_chains1)) %>% map(~(list_chains1[[.x]] %>%
      #          as.data.frame %>% mutate(!!t_par := logit(list_chains1[[.x]][, t_par], model_info[[s]]$other$upper_boundary[j]))))
      #          j <- j + 1
       # }
        df_s[[s]] <- postfull[[s]] %>% bind_rows %>% pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
            arrange(parameter) %>%
            mutate(model = model_info[[s]]$other$model_name)
    }

    df_s_final <- df_s %>% bind_rows
  #  df_s_final_trim <- df_s_final[seq(1, nrow(df_s_final), 10), ]

    pars <- df_s_final %>% pull(parameter) %>% unique
    models <- c("vac_base", "vac_m", "vac_s", "vac_t", "vac_mt", "vac_ms", "vac_ts", "vac_mts")

    df_s_final <- df_s_final %>% mutate(model = factor(model, levels = models))

    p1 <- list()
    legend_pos <- c(rep("right", 3), rep("none", 3))
    cols <- c("grey", "#1261A0", "#3895D3", "#58CCED", "#900D09", "#D21404", "#BC544B", "#4C9A2A")
    positions <- list(1:8, 1:8, 1:8, c(2, 5, 6, 8), c(3, 6, 7, 8), c(4, 5, 7, 8))
    j <- 1
    for(par in pars) {
        p1[[j]] <- df_s_final %>%
            filter(parameter == par) %>%
            ggplot() +
                geom_boxplot(aes(x = parameter, y = value, fill = model)) +
                theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
                theme(legend.position = legend_pos[j]) +
                scale_fill_manual(values = cols[positions[[j]]])
        j <- j + 1
    }
    wrap_plots(p1) + plot_annotation(tag_levels = c("A")) + plot_layout(guides = "collect")

    savepath <- here::here("outputs", model_info[[1]]$other$study$study_name_short, "fits", file, "sum_figs")
    ggsave(file = paste0(savepath, "/", "posteriors", ".pdf"))
}


plot_ar <- function(postfull, model_info) {

    inf_chain <- postfull$inf_chain

    ## Look at inferred attack rates
    p_ar <- plot_attack_rates(inf_chain,
        model_info$create_post$titre_data,
        model_info$create_post$strain_isolation_times,
        pad_chain=TRUE, plot_den = TRUE,
        prior_pars = list(prior_version = model_info$create_post$prior_version,
        alpha = model_info$create_post$par_tab[model_info$create_post$par_tab$names == "alpha", "values"],
        beta = model_info$create_post$par_tab[model_info$create_post$par_tab$names == "beta", "values"]))
    ggsave(p_ar, filename = paste0(model_info$other$output_file, model_info$other$study$study_name_short, "/", model_info$other$model_name, "/ar_plot.pdf"))
}

calc_mpsrf <- function(postfull, model_info) {
    minlen <- (postfull %>% map(~nrow(.x))) %>% unlist %>% min
    postfulltrim <- postfull %>% map(~.x[1:minlen, ])
    postfullmcmc <- postfulltrim %>% map(~mcmc(.x))
    data.frame(
        model_name = model_info$other$model_name,
        mpsrf = gelman.diag(postfullmcmc)$mpsrf
    )
}

save_plot_mpsrf <- function(mpsrf_df, model_info, file) {
    save(mpsrf_df, file = here::here("outputs", model_info$other$study$study_name_short, "fits", file, "sum_figs", "mpsrf.RDS"))
    models <- c("vac_base", "vac_m", "vac_s", "vac_t", "vac_mt", "vac_ms", "vac_ts", "vac_mts")
    cols <- c("grey", "#1261A0", "#3895D3", "#58CCED", "#900D09", "#D21404", "#BC544B", "#4C9A2A")
    mpsrf_df %>%
        mutate(model_name = factor(model_name, levels = models)) %>%
        ggplot() +
            geom_col(aes(x = mpsrf, y = model_name, fill = model_name)) +
            scale_fill_manual(values = cols) + 
            theme_minimal() + theme(legend.position = "none", axis.title.y = element_blank()) +
            labs(x = "MPSRF") + 
            geom_vline(xintercept = 1.1)

    ggsave(filename = here::here("outputs", model_info$other$study$study_name_short, "fits", file, "sum_figs", "mpsrf.pdf"))

}


calc_waic <- function(postfull, samp_no_list, model_info) {

    cal_post <- create_posterior_func(par_tab = model_info$create_post$par_tab,
                                    vaccination_histories = model_info$create_post$vac_history,
                                    titre_dat = model_info$create_post$titre_data,
                                    antigenic_map = model_info$create_post$antigenic_map,
                                    version = model_info$create_post$prior_version)

    inf_chain_thin <- postfull$inf_chain %>% filter(sampno %in% samp_no_list[[1]])
    pad_inf_chain_1 <- pad_inf_chain(inf_chain_thin)
    sample_nos <- pad_inf_chain_1 %>% pull(sampno) %>% unique
    NS <- length(sample_nos)
    N <- model_info$create_post$titre_data %>% pull(individual) %>% unique %>% length
    Np <- model_info$create_post$par_tab %>% nrow
    log_lik_ind <- matrix(, ncol = N, nrow = NS)
    log_lik_ind_list <- rep(list(log_lik_ind), 3)
    pad_inf_chain_array <- split(pad_inf_chain_1, pad_inf_chain_1$sampno)
    for (c in 1:4) { 
        for (s in 1:NS) {
            theta_chain_thin <- as.data.frame(postfull$theta_chain) %>% filter(sampno %in% samp_no_list[[1]])
            theta_vals_post <- theta_chain_thin[s, 3:(3 + Np)] %>% as.numeric

            inf_chain_post <- pad_inf_chain_array[[s]] %>%
                filter(chain_no == c) %>%
                arrange(i, j) %>%
                select(i, j, x) %>%
                pivot_wider(values_from = x, names_from = j) %>%
                select(!i) %>%
                as.matrix
            log_lik_ind[s, ] <- cal_post(theta_vals_post, inf_chain_post)[[1]]
        }
        log_lik_ind_list[[c]] <- log_lik_ind
    }

    df_output <- data.frame()
    pWAIC_list <- list()
    for (c in 1:4) {
        log_lik_ind <- log_lik_ind_list[[c]]
        log_lik_ind_temp <- lapply(1:N, function(i) log_lik_ind[, i][log_lik_ind[, i] > -Inf])
        llpd <- sapply(1:N, function(i) log(sum(exp(log_lik_ind_temp[[i]]))) - log(NS))
        pWAIC <- sapply(1:N, function(i) var(log_lik_ind_temp[[i]], na.rm = TRUE))
        WAIC <- -2 * (sum(llpd) - sum(pWAIC, na.rm = TRUE))
        WAICse <- sqrt(N * var(-2 * (llpd - pWAIC), na.rm = TRUE))
        df_output <- rbind(df_output, data.frame(WAIC = WAIC, WAICse = WAICse, chain = c, model = model_info$other$model_name ))
        pWAIC_list[[c]] <- pWAIC
    }

    list(df_output = df_output, pWAIC_list = pWAIC_list)
}


plot_titre_post <- function(postfull, model_info, file) {

    s <- 1
    n <- model_info$create_post$antigenic_map %>% nrow
    dist_mat <- matrix(, nrow = n, ncol = n)
    for (i in 1:n)
        for (j in 1:n)
            dist_mat[i, j] <- sqrt(sum((model_info$create_post$antigenic_map[i, 1:2] - model_info$create_post$antigenic_map[j, 1:2])^2))
    
    theta_vals_post <- postfull$theta_chain %>% apply(2, median)

    model_dist <- pmax(0, 1 - theta_vals_post["sigma1"] * dist_mat) %>% matrix(nrow = n, ncol = n)# matrix
    seninority <- sapply(1:4, function(x) max(0, 1.0 - theta_vals_post["tau"]*(x - 1))) # close to 1
    inf1 <- seninority[1] * theta_vals_post["mu"] * model_dist[3, 1:49]
    inf2 <- seninority[2] * theta_vals_post["mu"] * model_dist[10, 1:49]
    inf3 <- seninority[3] * theta_vals_post["mu"] * model_dist[45, 1:49]


    df_plot <- bind_rows(
        data.frame(titre = inf1, x = 1968:2016, inf_no = "1"),
        data.frame(titre = inf2, x = 1968:2016, inf_no = "2"),
        data.frame(titre = inf3, x = 1968:2016, inf_no = "3")

    )
    inf_heights <- c(3, 10, 45) %>% map(~(df_plot %>% filter(x == 1967 + .x) %>% pull(titre) %>% sum)) %>% unlist

    p1 <- df_plot %>%
        ggplot() +
            geom_col(aes(x = x, y = titre, fill = inf_no), alpha = 0.5) +
            coord_cartesian(ylim = c(0, 5)) +
            geom_segment(data = data.frame(x = 1968 + c(2, 9, 44), y = inf_heights),
                aes(x = x, xend = x ,
                    y = 0, yend = y), color = "red", size = 3) +
            scale_fill_manual(values = wes_palette("Darjeeling2")) + 
            labs(x = "Year", y = "Titre value", title = "Example of antigenic landscape") +
            theme(legend.position = "none") + 
            annotate("label", x = 1968 + 2, y = 4.7, size = 3.8, color = "gray20", lineheight = .9,
                label = "First infection") + 
            annotate("label", x = 1968 + 9, y = 2.7, size = 3.8, color = "gray20", lineheight = .9,
                label = "Second infection")  + 
            annotate("label", x = 1968 + 44, y = 2.7, size = 3.8, color = "gray20", lineheight = .9,
                label = "Third infection")

    p2 <- data.frame(x = 0:9,
        y = sapply(1:10, function(x) max(0, 1.0 - theta_vals_post["tau"]*(x - 1))))  %>%
        ggplot() + 
            geom_line(aes(x, y), size = 3) +
            labs(x = "Number of previous infections", y =  "Relative reduction in titre value") +
            theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))

    p3 <- p1 +
        annotation_custom(ggplotGrob(p2), xmin = 1990, xmax = 2016, ymin = 3, ymax = 5)

    inf_chain <- postfull$inf_chain

    sample_nos <- inf_chain %>% pull(sampno) %>% unique %>% length 
    sample_df <- inf_chain %>%
            group_by(year = j, sampno) %>%
            dplyr::summarise(x = sum(x / (3 * 50))) %>% as.data.frame %>% 
            filter(year < 49)

    p4 <- sample_df %>% 
            ggplot() + 
                    geom_boxplot(aes(x = as.character(year + 1967), y = x)) + 
                    scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
                    labs(x = "Year",  y = "Attack rate", title = "Estimated attack rate")
    
    p3 / p4
    savepath <- here::here("outputs", model_info$other$study$study_name_short, "fits", file, model_info$other$model_name)
    dir.create(savepath)
    ggsave(filename = paste0(savepath, "/", "schematic.pdf"))
}

plot_post_combine <- function(model_info, study, sample_yr) {

    sero_data <- study$sero_data
    part_data <- study$part_data
    vac_data <- study$vac_data
    titre_data <- convert_to_ss_df(sero_data, part_data, sample_yr)
 #   vac_history <- convert_to_ss_vac(sero_data, vac_data)
    # model data
    antigenic_map <- model_info$antigenic_map
    par_tab <- model_info$par_tab
    prior_version <- model_info$prior_version
    output_file <- model_info$output_file
    # trim antigenic map according to titre data
    antigenic_map[antigenic_map$inf_times >= min(titre_data$virus) & antigenic_map$inf_times <= max(titre_data$virus), ]
    strain_isolation_times <- unique(antigenic_map$inf_times)
    ## Maximum recordable log titre in these data is 8
    par_tab[par_tab$names == "MAX_TITRE", "values"] <- max(titre_data$titre)

    output_file_novac <- paste0(output_file, study$study_name_short, "/", "novac", "/")
    all_chains_novac <- load_mcmc_chains(location = output_file_novac, thin = 10, burnin = 1000000,
                               par_tab = par_tab, unfixed = FALSE, convert_mcmc = TRUE)
    list_chains_novac <- all_chains_novac$theta_list_chains

    output_file_vac <- paste0(output_file, study$study_name_short, "/", "vac", "/")
    all_chains_vac <- load_mcmc_chains(location = output_file_vac, thin = 10, burnin = 1000000,
                               par_tab = par_tab, unfixed = FALSE, convert_mcmc = TRUE)
    list_chains_vac <- all_chains_vac$theta_list_chains

    list_chains <- c(list_chains_novac[1:3], list_chains_vac[1:3])
    pars_plot <- c("mu", "sigma1", "tau")
## Look at diagnostics for the free parameters
    list_chains1 <- lapply(list_chains, function(x) x[, c(pars_plot, "total_infections",
                                                     "lnlike", "prior_prob")])
    p_theta_dens <- mcmc_trace(list_chains1, color_chains = TRUE)
    ggsave(p_theta_dens, filename = paste0(output_file, study$study_name_short,  "/post_plot.pdf"))
}


plot_ar_alt <- function(model_info, study, sample_yr, file) {

    sero_data <- study$sero_data
    part_data <- study$part_data
    vac_data <- study$vac_data
    titre_data <- convert_to_ss_df(sero_data, part_data, sample_yr)
#    vac_history <- convert_to_ss_vac(sero_data, vac_data)
    # model data
    antigenic_map <- model_info$antigenic_map
    par_tab <- model_info$par_tab
    prior_version <- model_info$prior_version
    output_file <- model_info$output_file
        output_file_novac <- paste0(output_file, study$study_name_short, "/", "novac", "/")
    all_chains_novac <- load_mcmc_chains(location = output_file_novac, thin = 10, burnin = 1000000,
                               par_tab = par_tab, unfixed = FALSE, convert_mcmc = TRUE)
    list_chains_novac <- all_chains_novac$theta_list_chains

    output_file_vac <- paste0(output_file, study$study_name_short, "/", "vac", "/")
    all_chains_vac <- load_mcmc_chains(location = output_file_vac, thin = 10, burnin = 1000000,
                               par_tab = par_tab, unfixed = FALSE, convert_mcmc = TRUE)
    list_chains_vac <- all_chains_vac$theta_list_chains
    # trim antigenic map according to titre data
    antigenic_map[antigenic_map$inf_times >= min(titre_data$virus) & antigenic_map$inf_times <= max(titre_data$virus), ]
    strain_isolation_times <- unique(antigenic_map$inf_times)
    ## Maximum recordable log titre in these data is 8
    par_tab[par_tab$names == "MAX_TITRE","values"] <- max(titre_data$titre)

    output_file <- paste0(output_file, study$study_name_short, "/", file, "/")
    all_chains <- load_mcmc_chains(location=output_file, thin=10, burnin=100000,
                               par_tab=par_tab,unfixed=FALSE,convert_mcmc=TRUE)
    inf_chain <- all_chains$inf_chain

    ## Look at inferred attack rates
    inf_chain1 <- all_chains_novac$inf_chain
    inf_chain1 <- pad_inf_chain(inf_chain1)
    n_alive1 <- get_n_alive(titre_data, strain_isolation_times)
    data.table::setkey(inf_chain1, "sampno", "j","chain_no")
    tmp <- inf_chain1[, list(V1 = sum(x)), by = key(inf_chain1)]

    quantiles <- ddply(tmp, ~j, function(x) quantile(x$V1, c(0.025, 0.25, 0.5, 0.75,  0.975)))
    colnames(quantiles) <- c("j", "lower", "lower_50","median","upper_50","upper")
    quantiles[c("lower", "lower_50", "median","upper_50","upper")] <- quantiles[c("lower", "lower_50","median","upper_50","upper")] / n_alive1[quantiles$j]
    quantiles$year <- strain_isolation_times[quantiles$j]
    quantiles$taken <- quantiles$year %in% unique(titre_data$samples)
    quantiles$vac_status <- c(rep("Model without vaccinated persons removed",dim(quantiles)[1]))

    # vaccinated
    inf_chain2 <- all_chains_vac$inf_chain
    inf_chain2 <- inf_chain2[inf_chain2$chain_no == 1,]
    inf_chain2 <- pad_inf_chain(inf_chain2)
    n_alive2 <- get_n_alive(titre_data, strain_isolation_times)
    data.table::setkey(inf_chain2, "sampno", "j","chain_no")
    tmp <- inf_chain2[, list(V1 = sum(x)), by = key(inf_chain2)]

    quantiles2 <- ddply(tmp, ~j, function(x) quantile(x$V1, c(0.025, 0.1, 0.5, 0.9,  0.975)))
    colnames(quantiles2) <- c("j", "lower", "lower_50","median","upper_50","upper")
    quantiles2[c("lower", "lower_50","median","upper_50","upper")] <- quantiles2[c("lower", "lower_50","median","upper_50","upper")] / n_alive2[quantiles2$j]
    quantiles2$year <- strain_isolation_times[quantiles2$j]
    quantiles2$taken <- quantiles2$year %in% unique(titre_data$samples)
    quantiles2$vac_status <- c(rep("Model with vaccinated persons removed", dim(quantiles2)[1]))


    quantiles_all <- rbind(quantiles, quantiles2)
    ## Colour depending on vac_status
    colour_fills_unvac <- c("#E69F00","#0072B2")
    colour_fills_age <- c("#CC79A7","#009E73","#56B4E9")

    strain_isolation_times1 <- strain_isolation_times

    ymax <- 0.5
    quantiles_all$vac_status <- factor(quantiles_all$vac_status, levels=c("Model without vaccinated persons removed","Model with vaccinated persons removed"))
    quantiles_all <- quantiles_all %>% na.omit
    p <- ggplot(quantiles_all) +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = vac_status), alpha = 0.25) +
    geom_ribbon(aes(x = year, ymin = lower_50, ymax = upper_50, fill = vac_status), alpha = 0.5) +
    geom_line(aes(x = year, y = median, colour = vac_status), size=0.75) +
    geom_point(aes(x = year, y = median, colour = vac_status), size = 0.75) +
   # scale_y_continuous(limits = c(-0.005, ymax), expand = c(0, 0),breaks=seq(0,ymax,by=0.05)) +
   # scale_x_continuous(expand = c(0, 0), breaks = strain_isolation_times1, labels = labels2,
   #                     limits=c(min(strain_isolation_times-0.1),max(strain_isolation_times+0.1))) +
    theme_pubr() +
    scale_fill_manual(values = colour_fills_unvac) +
    scale_color_manual(values = colour_fills_unvac) +
    ylab("Estimated per capita\n incidence (per quarter)") +
    xlab("Time of virus circulation") +
    theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.title = element_text(size = 10),
            legend.title = element_blank(),
            legend.text = element_text(size = 8,family = "sans"),
            legend.position = c(0.7, 0.99),
            legend.direction = "horizontal",
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.key = element_rect(color = NA),
            legend.background = element_blank(),
            legend.margin = margin(6, 6, 6, 6),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.margin = margin(l = 10,r = 5,t = 5))
    p


}
