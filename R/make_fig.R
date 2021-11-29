plot_schematic_mu <- function(mu_vac, subtitle1) {
    colors_plt <- wes_palette("Royal1", n = 2) %>% as.character

    df_plt <- data.frame(
        type = c("Natural infection", "Vaccination"),
        height = c(1, 1 * mu_vac)
    )

    if (mu_vac == 1) {
        breaks_y <- c(0, 1)
        labels_y <- c("", TeX(r'($\mu$)'))
    } else {
        breaks_y <- c(0, 1, 1 * mu_vac)
        labels_y <- c("", TeX(r'($\mu$)'), TeX(r'($\mu\times\mu_{vac}$)'))
    }
    p1 <- df_plt %>%
        ggplot() +
            geom_col(aes(x = type, y = height, fill = type), width = 0.4) +
            theme_minimal() +
            scale_fill_manual(values = colors_plt) +
            theme(  panel.grid.minor = element_blank(),
                    legend.position = "none") +
            labs(x = "Type of exposure", y = "Titre boost to infecting strain", subtitle = subtitle1) +
            scale_y_continuous(breaks = breaks_y, labels = labels_y, limits = c(0, 1.5))
    p1
}

plot_schematic_mu_both <- function() {
    p1 <- plot_schematic_mu(1, "Exposure independent")
    p2 <- plot_schematic_mu(0.5, "Exposure dependent")
    p1 + p2
    ggsave(here::here("inst", "figs", "schematic", "mu.pdf"), width = 10, height = 5)
}

plot_schematic_tau <- function(tau_vac, subtitle1) {
    tau <- 0.04
    colors <- wes_palette("Royal1", n = 2) %>% as.character
    colors_plt <- c(colors[1], colors[1], colors[2], colors[1], colors[1], colors[2], colors[1], colors[1],
            colors[2], colors[1])
    data_plt <- data.frame(
        infection_no = factor(as.character(1:10), levels = as.character(1:10)),
        type = c("Natural infection", "Natural infection", "Vaccination", "Natural infection", "Natural infection",
            "Vaccination", "Natural infection", "Natural infection",   "Vaccination", "Natural infection"),
        height = c(1, (1 - 1 * tau), (1 - 2 * tau * tau_vac), (1 - 3 * tau), (1 - 4 * tau), (1 - 5 * tau * tau_vac),
            (1 - 6 * tau), (1 - 7 * tau), (1 - 8 * tau * tau_vac), (1 - 9 * tau)),
        colors = c(colors[1], colors[1], colors[2], colors[1], colors[1], colors[2], colors[1], colors[1],
            colors[2], colors[1])
    )
    labels_plt <- data.frame(
        x = c(9, 9),
        y = c(1.1, 0.5),
        #labels = c(TeX(r'($\tau$)'), TeX(r'($\tau_{vac}$)')),
        color = c(colors[1], colors[2])
    )

    p2 <- data_plt %>%
        ggplot() +
            geom_abline(slope = -0.04 * tau_vac, intercept = 1.3, color = colors[2]) +
            geom_text(data = labels_plt, aes(x = x, y = y),  label = c(TeX(r'($\tau$)'), TeX(r'($\tau\times\tau_{vac}$)'))) +
            geom_abline(slope = -0.04, intercept = 1.3, color = colors[1]) +
            geom_col(aes(x = infection_no, y = height, fill = infection_no)) +
            theme_minimal() +
            scale_fill_manual(values = colors_plt) +
            theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                panel.grid.minor = element_blank(), legend.position = "none") +
            scale_y_continuous(limits = c(0, 1.5)) +
            scale_x_discrete(breaks = seq(2, 10, 2)) +
            labs(x = "Number of exposures", subtitle = subtitle1)

    p2
}

plot_schematic_tau_both <- function() {
    p1 <- plot_schematic_tau(1, "Exposure independent")
    p2 <- plot_schematic_tau(2, "Exposure dependent")
    p1 + p2
    ggsave(here::here("inst", "figs", "schematic", "tau.pdf"), width = 10, height = 5)
}

plot_schematic_sigma <- function(sigma_vac, title_text) {
    colors_plt <- wes_palette("Royal1", n = 2) %>% as.character
    df_plot_1 <- data.frame(
        y = dnorm(seq(-1, 1, 0.2), 0, 1),
        x = seq(-1, 1, 0.2),
        color = "Natural infection"
    )
    df_plot_1$y <- (df_plot_1$y) / (max(df_plot_1$y))
    df_plot_2 <- data.frame(
        y = dnorm(seq(-1, 1, 0.2), 0, 1 * sigma_vac) * peak,
        x = seq(-1, 1, 0.2) + 3,
        color = "Vaccination"
    )
    df_plot_2$y <- (df_plot_2$y) / (max(df_plot_2$y))
    df_plot <- bind_rows(df_plot_1, df_plot_2)
    if (sigma_vac == 1) {
        label_text <- c(TeX(r'($\sigma$)'), TeX(r'($\sigma$)'))
    } else {
        label_text <- c(TeX(r'($\sigma$)'), TeX(r'($\sigma\times\sigma_{vac}$)'))
    }
    labels_plt <- data.frame(x = c(0.5, 3 + sigma_vac / 2), y = c(1.15, 1.15))

    df_plot %>%
        ggplot() +
            geom_col(aes(x = x, y = y, fill = color)) +
            scale_fill_manual(values = colors_plt) +
            geom_segment(aes(x = 0, y = 1.1, xend = 1, yend = 1.1), arrow = arrow(length = unit(0.5, "cm"))) +
            geom_segment(aes(x = 3, y = 1.1, xend = 3 + sigma_vac, yend = 1.1), arrow = arrow(length = unit(0.5, "cm"))) +
            geom_text(data = labels_plt, aes(x = x, y = y), label = label_text) +
            theme(axis.title.y = element_blank(),axis.title.x = element_blank(), axis.text.y = element_blank(),
                panel.grid.minor = element_blank(), legend.position = "none") +
                scale_y_continuous(limits = c(0, 1.5)) +
                scale_x_continuous(breaks = c(0, 3), labels = c("Natural infection", "Vaccination")) + 
                labs(title = title_text)
            
}

plot_schematic_sigma_both <- function() {
    p1 <- plot_schematic_sigma(1, "Exposure independent")
    p2 <- plot_schematic_sigma(0.3, "Exposure dependent")
    p1 + p2
    ggsave(here::here("inst", "figs", "schematic", "sigma.pdf"), width = 10, height = 5)
}

plot_schematic <- function() {

    figA <- plot_schematic_1(1, 1, "Exposure independent", "Exposure independent", "Boosting", "Antigenic seniority")
    figB <- plot_schematic_1(0.5, 1, "Exposure independent", "Exposure dependent")
    figC <- plot_schematic_1(1, 3, "Exposure dependent", "Exposure independent")
    figD <- plot_schematic_1(0.5, 3, "Exposure dependent", "Exposure dependent")

    (figA) / (figB) / (figC) / (figD) +
        plot_annotation(tag_levels = 'A')
    ggsave(filename = here::here("inst", "figs", "schematic", "fig1.pdf"), height = 10, width = 10)
}

plot_schematic_new_build <- function() { 
    data.frame(
        type = c(rep("Existing antibody response", 21), rep("De novo antibody response", 21)),
        y = c(-10:10 %>% logit(1), 1 - (-10:10 %>% logit(1))),
        x = c(0:20, 0:20)
    ) %>%
        ggplot() +
            geom_line(aes(x = x, y = y, color = type)) + 
            labs(y = "Probability", x = "Titre value", color = "Type of immune response") + 
            theme(axis.text.x = element_blank(), legend.position = "top")
    ggsave(file = here::here("inst", "figs", "schematic", "fig_denovo.pdf"), height = 5, width = 5)
}