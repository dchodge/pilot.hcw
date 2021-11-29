make_study <- function(study_name_full, study_name_short, part_data = NA, sero_data = NA, vac_data = NA) {
    study <- list(
        study_name_full = study_name_full,
        study_name_short = study_name_short,
        part_data = part_data,
        sero_data = sero_data,
        vac_data = vac_data
    )
    class(study) <- append(class(study), "study")
    return(study)
}

plot_vac_hist <- function(study_class) {
    vac_data <- study_class$vac_data
    p1 <- vac_data %>%
        ggplot(aes(x = vac_yr, y = pid, fill = response)) + 
        geom_raster() +
        scale_fill_manual(values = c("black", "gray", "red")) + 
        theme(aspect.ratio = 0.5, axis.text.y = element_text(size = 1)) +
        labs(x = "Influenza season", y = "PID", fill = "Vaccinated?", 
            subtitle = "Each participant")
    
    p2 <- vac_data %>%
        ggplot(aes(x = vac_yr)) +
            geom_bar(aes(fill = response), position = "fill") +
            scale_fill_manual(values = c("black", "gray", "red")) +
            theme_bw() +
            theme(aspect.ratio = 0.2, legend.position = "none") + 
            scale_y_continuous(labels = scales::percent) + 
            labs(x = "Influenza season", y = "Proportion", fill = "Vaccinated?", subtitle = "Across all participants")

    p1 / p2 + plot_layout(guides = "collect") + plot_annotation(title = study_class$study_name_full, tag_levels = "A")

    ggsave(here::here("inst", "figs", "vac", paste0(study_class$study_name_short, "_1.pdf")))

}

plot_sero <- function(study_class) {
    sero_data <- study_class$sero_data
    virus_order <- sero_data %>%
        select(virus_name, virus_isol_yr) %>%
        arrange(virus_isol_yr) %>%
        pull(virus_name) %>%
        unique
    sero_data$virus_name <- factor(sero_data$virus_name, levels =  virus_order)
    sero_data$sample_time <- factor(sero_data$sample_time, levels =  sero_data$sample_time %>% unique %>% sort)

    sero_data_pre_vac <- sero_data %>%
        group_by(virus_name, titre_value, sample_time) %>%
        mutate(freq = n())

    sero_data_pre_vac_wmean <- sero_data_pre_vac %>%
        ungroup %>%
        group_by(virus_name, sample_time) %>%
        dplyr::summarise(titre_value = weighted.mean(titre_value, freq, rm.na = TRUE))

    p1 <- sero_data_pre_vac %>%
        ggplot(aes(x = virus_name, y = titre_value)) +
            geom_point(aes(fill = freq, size = freq), alpha = 0.1, color = "white", shape = 21) +
            geom_point(data = sero_data_pre_vac_wmean,  color = "black", fill = "red",
                alpha = 0.8, size = 4, shape = 23) +
            scale_size(guide = "none") +
            scale_fill_continuous(type = "gradient") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90)) +
            labs(x = "Strian name", y = "Titre", fill = "Frequency", subtitle = "Pre-vaccination titre values") +
            facet_wrap(vars(sample_time), nrow = (sero_data_pre_vac$sample_time %>% unique %>% length))

    p2 <- sero_data_pre_vac_wmean %>%
        ggplot(aes(x = virus_name)) +
            geom_point(aes(y = titre_value, fill = sample_time), shape = 21, color = "white", size = 3,
                position = position_dodge(width = 0.4)) +
            geom_line(aes(y = titre_value, color = sample_time, group = sample_time), 
                position = position_dodge(width = 0.4)) + 
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90)) +
            labs(x = "Strian name", y = "Mean titre", fill = "Sample time",
            subtitle = "Titre values pre- and post-vaccination by stain")

    p3 <- sero_data_pre_vac_wmean %>%
        ggplot(aes(x = as.numeric(as.character(sample_time)))) +
            geom_line(aes(y = titre_value, color = virus_name), size = 0.2)  +
            geom_point(aes(y = titre_value, color = virus_name), shape = 1, size = 2,
                position = position_dodge(width = 0.4)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90)) +
            labs(x = "Time post-vaccination", y = "Mean titre", color = "Strain name",
            subtitle = "Titre values pre- and post-vaccination by time post-vac")
    p1 + (p2 / p3) + plot_annotation(title = study_class$study_name_full, tag_levels = "A")
    ggsave(here::here("inst", "figs", "sero", paste0(study_class$study_name_short, "_titre_mean.pdf")), 
        height = 12, width = 16, dpi = 50)
    ggsave(here::here("inst", "figs", "sero", paste0(study_class$study_name_short, "_titre_mean.png")), 
        height = 12, width = 16, dpi = 50)

}

plot_part <- function(study_class) {
    part_data <- study_class$part_data  %>% unique

    p1 <- part_data %>%
        ggplot(aes(x = yob)) +
            geom_histogram() +
            theme_bw() +
            theme(aspect.ratio = 0.3) +
            labs(x = "Year of birth", subtitle = "Participant info")
    p1 + plot_annotation(title = study_class$study_name_full)
    ggsave(here::here("inst", "figs", "part", paste0(study_class$study_name_short, "_yob.pdf")))
}

# plot titres agaisnt year of birth, how do people born in different decades respond to vaccination?
plot_yob_titre <- function(study_class) {
    part_data <- study_class$part_data %>% unique
    sero_data <- study_class$sero_data
    combine_data <- left_join(sero_data, part_data)

    df_plot <- combine_data %>% 
        group_by(pid, virus_name) %>%
        mutate(titre_diff = titre_value - lag(titre_value, default = titre_value[1])) %>%
        group_by(virus_id, virus_name, virus_isol_yr, yob, sample_time) %>%
        dplyr::summarise(titre_value = mean(titre_value, na.rm = TRUE), titre_diff = mean(titre_diff, na.rm = TRUE), n = n()) %>%
        ungroup
        
    p1 <- df_plot %>% filter(sample_time == 0) %>%
        ggplot(aes(x = yob, y = titre_value)) +
            geom_point(aes(color = virus_name), size = 0.3, alpha = 0.4) +
            theme_bw() +
            scale_color_discrete(guide = "none") +
            geom_smooth(aes(color = virus_name),  se = FALSE, size = 0.5) +
            labs(x = "Year of birth", y = "Pre-vac titre value", color = "Strain")

    p2 <- df_plot %>% filter(sample_time != 0)  %>%
        ggplot(aes(x = yob, y = titre_diff)) +
            geom_point(aes(color = virus_name), size = 0.5, alpha = 0.4) +
            theme_bw() +
            geom_smooth(aes(color = virus_name),  se = FALSE, size = 0.5) + 
            facet_wrap(vars(sample_time), nrow = 3) + 
            labs(x = "Year of birth", y = "Change in titre value from pre-vac", color = "Strain")
    
    p1 / p2 + plot_layout(heights = c(1, 3), guides = "collect")
    ggsave(here::here("inst", "figs", "sero", paste0(study_class$study_name_short, "_titre_yob.pdf")))

}

plot_yob_titre_spline <- function(study_class) {
    part_data <- study_class$part_data %>% unique
    sero_data <- study_class$sero_data
    combine_data <- left_join(sero_data, part_data, by = "pid") %>% na.omit

    df_plot <- combine_data %>%
        group_by(pid, virus_name) %>%
        mutate(titre_diff = titre_value - lag(titre_value, default = titre_value[1])) %>% ungroup

    df_plot_base <- df_plot %>%
        filter(sample_time == 0) %>%
        select(virus_isol_yr, yob, titre_value) %>%
        arrange(virus_isol_yr, yob, titre_value)

   df_plot_postvac <- df_plot %>%
        filter(sample_time != 0) %>%
        select(virus_isol_yr, yob, sample_time, titre_diff) %>%
        arrange(virus_isol_yr, yob, sample_time, titre_diff)
    
# Thin plate spline regression
    get_smoothed_spline <- function(df_data, titre_val, time) {
        require(fields)
        X <- df_data %>% na.omit %>% select(yob, virus_isol_yr) %>% as.matrix
        y <- df_data %>% na.omit %>% select(!!titre_val) %>% as.matrix()
        fit <- Tps(X, y)
        yob_grid <- c(min(X[, 1]):max(X[, 1]))
        virus_grid <- c(min(X[, 2]):max(X[, 2]))

        grid_points <- expand.grid(x = yob_grid, virus_grid)
        df_surface <- predict(fit, grid_points) %>% as.data.frame %>%
            mutate(x_yob = rep(yob_grid, length(virus_grid)), 
                x_virus = unlist(virus_grid %>% map(~rep(.x, length(yob_grid))  )) ) %>% dplyr::rename(titre = V1) %>% mutate(time = time)
    }

    df_plot_1 <- get_smoothed_spline(df_plot_base, "titre_value", 0)
    df_plot_rest <- unique(df_plot_postvac$sample_time) %>% 
        map(~df_plot_postvac %>% filter(sample_time == .x) %>% get_smoothed_spline("titre_diff", .x)) %>% bind_rows

    require(metR)
    p1 <- df_plot_1 %>%
        ggplot(aes(x = x_yob, y = x_virus)) +
            geom_raster(aes(fill = titre)) +
            geom_contour(aes(z = titre), colour = "white", bins = 20) +
            geom_text_contour(aes(z = titre), stroke = 0.2) +
            scale_fill_gradient(low = "gray", high = "blue") +
            geom_abline(size = 2) + theme_minimal() +
            labs(x = "Year of birth", y = "Year of virus isolation", fill = "Titre")

    p2 <- df_plot_rest %>%
        ggplot(aes(x = x_yob, y = x_virus)) +
            geom_raster(aes(fill = titre)) +
            scale_fill_gradient2(aes(fill = titre)) +
            geom_contour(aes(z = titre), colour = "black", bins = 20) +
            geom_text_contour(aes(z = titre), stroke = 0.2) +
            geom_abline(size = 2) + theme_minimal() +
            labs(x = "Year of birth", y = "Year of virus isolation", fill = "Titre change") + 
            facet_wrap(vars(time))

    p1 / p2 + plot_annotation(tag_levels = "A", title = study_class$study_name_full)

    ggsave(here::here("inst", "figs", "sero", paste0(study_class$study_name_short, "_titre_yob_spline.pdf")))
}

plot_vaccine_impact_spline <- function(study_class) {
    part_data <- study_class$part_data %>% unique
    vac_data <- study_class$vac_data %>% unique
    sero_data <- study_class$sero_data

    vac_data_sum <- vac_data %>%
        filter(response %in% c("Yes", "No")) %>%
        group_by(pid, response) %>%
        dplyr::summarise(no_vac = n()) %>% ungroup %>% filter(response == "Yes") %>% select(pid, no_vac)


    combine_data <- left_join(sero_data, vac_data_sum, by = "pid")

    df_plot <- combine_data %>%
        group_by(pid, virus_name) %>%
        mutate(titre_diff = titre_value - lag(titre_value, default = titre_value[1])) %>% ungroup

    df_plot_base <- df_plot %>%
        filter(sample_time == 0) %>%
        select(virus_isol_yr, no_vac, titre_value) %>%
        arrange(virus_isol_yr, no_vac, titre_value)

   df_plot_postvac <- df_plot %>%
        filter(sample_time != 0) %>%
        select(virus_isol_yr, no_vac, sample_time, titre_diff) %>%
        arrange(virus_isol_yr, no_vac, sample_time, titre_diff)
    
# Thin plate spline regression
    get_smoothed_spline <- function(df_data, titre_val, time) {
        require(fields)
        X <- df_data %>% na.omit %>% select(no_vac, virus_isol_yr) %>% as.matrix
        y <- df_data %>% na.omit %>% select(!!titre_val) %>% as.matrix()
        yob_grid <- c(min(X[, 1]):max(X[, 1]))
        virus_grid <- c(min(X[, 2]):max(X[, 2]))
                if (length(yob_grid) == 1){
            return(data.frame())
        }
        fit <- Tps(X, y)


        grid_points <- expand.grid(x = yob_grid, virus_grid)
        df_surface <- predict(fit, grid_points) %>% as.data.frame %>%
            mutate(x_yob = rep(yob_grid, length(virus_grid)), 
                x_virus = unlist(virus_grid %>% map(~rep(.x, length(yob_grid))  )) ) %>% dplyr::rename(titre = V1) %>% mutate(time = time)
    }

    df_plot_1 <- get_smoothed_spline(df_plot_base, "titre_value", 0)
    df_plot_rest <- unique(df_plot_postvac$sample_time) %>% 
        map(~df_plot_postvac %>% filter(sample_time == .x) %>% get_smoothed_spline("titre_diff", .x)) %>% bind_rows

    require(metR)
    p1 <- df_plot_1 %>%
        ggplot(aes(x = x_yob, y = x_virus)) +
            geom_raster(aes(fill = titre)) +
            geom_contour(aes(z = titre), colour = "white", bins = 20) +
            geom_text_contour(aes(z = titre), stroke = 0.2) + 
            scale_fill_gradient(low = "gray", high = "blue") +
            geom_abline(size = 2) + theme_minimal() +
            labs(x = "Number of previous vaccinations", y = "Year of virus isolation", fill = "Titre")

    p2 <- df_plot_rest %>%
        ggplot(aes(x = x_yob, y = x_virus)) +
            geom_raster(aes(fill = titre)) +
            geom_contour(aes(z = titre), colour = "black", bins = 20) +
            geom_text_contour(aes(z = titre), stroke = 0.2) + 
            scale_fill_gradient2() +
            geom_abline(size = 2) + theme_minimal() +
            labs(x = "Number of previous vaccinations", y = "Year of virus isolation", fill = "Titre change") + 
            facet_wrap(vars(time))

    p1 / p2 + plot_annotation(tag_levels = "A", title = study_class$study_name_full)

    ggsave(here::here("inst", "figs", "vac", paste0(study_class$study_name_short, "_titre_vac_spline.pdf")))
}