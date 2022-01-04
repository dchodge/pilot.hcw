load_point_ll <- function(model_info, folder_name_1, cross_sec = TRUE, burnin = TRUE) {
    if (cross_sec) {
        folder_name_2 <- "cross_sec"
    } else {
        folder_name_2 <- "full"
    }
    location <- here::here("outputs", model_info$other$study$study_name_short, "fits", folder_name_1, folder_name_2, model_info$other$model_name)
    chains <- Sys.glob(file.path(location, "*_indiv_ll.csv"))
    message(cat("Chains detected: ", length(chains), sep = "\t"))
    if (length(chains) < 1) {
        message("Error - no chains found")
        return(NULL)
    }

    read_chains <- lapply(chains, function(x) read.csv(x, header = FALSE))
    L <- (read_chains[[1]] %>% nrow)
    if (burnin) {
        B <- round(L / 2)
    } else {
        B <- 1
    }
    read_chains_trim <- lapply(read_chains, function(x) x[B:L, ])
    ll_point_mat <- read_chains_trim %>% bind_rows %>% as.matrix
    ll_point_mat
}

save_plot_waic_cross <- function(all_models_hcw_pre_cross, file.path = "dis_dep", cross_sec = TRUE) {
    waics <- vector(mode = "list", length = length(all_models_hcw_pre_cross))
    for (i in 1:length(all_models_hcw_pre_cross)) {
        model_info <- all_models_hcw_pre_cross[[i]]
        ll_point_mat <- load_point_ll(model_info, file.path, cross_sec, burnin = TRUE)
        waics[[i]] <- waic(ll_point_mat)
    }

    if (cross_sec) {
        file.path.2 <- "cross_sec"
    } else {
        file.path.2 <- "full"
    }

    save(waics, file =  here::here("outputs", "hcw_pre", "fits", file.path, file.path.2, "sum_figs", "waics.RDS"))

    data.frame(
        y = 1:8,
        name = c("vac_base", "vac_m", "vac_s", "vac_t", "vac_ms", "vac_mt", "vac_ts", "vac_mts"),
        waic_mean = waics[2:9] %>% map(~.x$estimates[3, 1]) %>% unlist,
        waic_sd = waics[2:9] %>% map(~.x$estimates[3, 2]) %>% unlist
    ) %>% 
        ggplot(aes(y = name)) + 
            geom_linerange(aes(xmax = waic_mean + 2 * waic_sd, xmin = waic_mean - 2 * waic_sd)) + 
            geom_point(aes(x = waic_mean), shape = 21, size = 3, fill = "red", color = "black") + 
            labs(x = "WAIC", y = "Model structure")

    ggsave(filename = here::here("outputs", "hcw_pre", "fits",  file.path, file.path.2,  "sum_figs", "waics.pdf"))
}

save_plot_waic_full <- function(all_models_hcw_pre_cross, file.path = "dis_dep", cross_sec = TRUE) {
    waics <- vector(mode = "list", length = length(all_models_hcw_pre_cross))
    for (i in 1:length(all_models_hcw_pre_cross)) {
        model_info <- all_models_hcw_pre_cross[[i]]
        ll_point_mat <- load_point_ll(model_info, file.path, cross_sec, burnin = TRUE)
        waics[[i]] <- waic(ll_point_mat)
    }

    if (cross_sec) {
        file.path.2 <- "cross_sec"
    } else {
        file.path.2 <- "full"
    }

    save(waics, file =  here::here("outputs", "hcw_pre", "fits", file.path, file.path.2, "sum_figs", "waics.RDS"))

    data.frame(
        y = 1:10,
        name = c("base", "m1", "m2", "m3", "m4"),
        waic_mean = waics[1:5] %>% map(~.x$estimates[3, 1]) %>% unlist,
        waic_sd = waics[1:5] %>% map(~.x$estimates[3, 2]) %>% unlist
    ) %>% 
        ggplot(aes(y = name)) + 
            geom_linerange(aes(xmax = waic_mean + 2 * waic_sd, xmin = waic_mean - 2 * waic_sd)) + 
            geom_point(aes(x = waic_mean), shape = 21, size = 3, fill = "red", color = "black") + 
            labs(x = "WAIC", y = "Model structure")

    ggsave(filename = here::here("outputs", "hcw_pre", "fits",  file.path, file.path.2,  "sum_figs", "waics.pdf"))
}