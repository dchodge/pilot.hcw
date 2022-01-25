#' Function to make the antigenic maps (no vaccination invovled)
#' @export
make_antigenic_maps_default_n <- function(antigenic_map_melted, theta) {

    antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
    antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"] * theta["sigma1_vac"])

    return(list(long = antigenic_map_long, short = antigenic_map_short))
}

#' Function to make the antigenic maps (including vaccination)
#' @export
make_antigenic_maps_vac1_n <- function(antigenic_map_melted, theta) {

    antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
    antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])
    antigenic_map_long_vac <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"] * theta["sigma1_vac"])
    antigenic_map_short_vac <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2_vac"])

    return(list(long = antigenic_map_long, short = antigenic_map_short, long_vac = antigenic_map_long_vac, short_vac = antigenic_map_short_vac))
}

#' Function to make the antigenic maps (including vaccination)
#' @export
make_antigenic_maps_vac2_n <- function(antigenic_map_melted, theta) {
    require(serosolver)
    antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
    antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, 0.0)
    antigenic_map_long_vac <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"] * theta["sigma1_vac"])
    antigenic_map_short_vac <- create_cross_reactivity_vector(antigenic_map_melted, 0.0)

    return(list(long = antigenic_map_long, short = antigenic_map_short, long_vac = antigenic_map_long_vac, short_vac = antigenic_map_short_vac))
}