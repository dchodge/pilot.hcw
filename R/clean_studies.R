library(here)
library(magrittr)
library(tidyverse)
library(lubridate)
library(readxl)

import_hcw_pre <- function() {
    hcw_pre_vac <- read.csv(here::here("inst", "extdata", "hcw_pre_vac.csv"))
    hcw_pre_vac_clean <- hcw_pre_vac %>%
        pivot_longer(!c(Participant.ID, comorbidites), names_to = "vac_yr", values_to = "response")
    hcw_pre_vac_clean$comorbidites <- factor(hcw_pre_vac_clean$comorbidites, levels = c("None", "Yes"))
    hcw_pre_vac_clean$vac_yr <- factor(hcw_pre_vac_clean$vac_yr,
        levels = c("vacc2006", "vacc2007", "vacc2008", "vacc2009", "vacc2010", "vacc2011",
            "vacc2012", "vacc2013", "vacc2014", "vacc2015"),
        labels = c("2006", "2007", "2008", "2009", "2010", "2011",
            "2012", "2013", "2014", "2015")
            )
    hcw_pre_vac_clean$response <- factor(hcw_pre_vac_clean$response,
        levels = c("No", "Don't know", "Yes"))
    hcw_pre_vac_data <- hcw_pre_vac_clean %>% dplyr::rename(pid = Participant.ID)
    save(hcw_pre_vac_data, file = here::here("data", "hcw_pre", "vac_data.RDS"))

    hcw_pre_data <- read.csv(here::here("inst", "extdata", "hcw_pre_data.csv"))
    hcw_pre_part_data <- hcw_pre_data %>% select(pid = PID, yob = year_of_birth)
    save(hcw_pre_part_data, file = here::here("data", "hcw_pre", "part_data.RDS"))

    hcw_pre_sero_data <- hcw_pre_data %>% select(pid = PID, virus_id = virus, virus_name = Short_Name, virus_isol_date = IsolDate, 
        virus_isol_yr = Year, titre_pre_vac = L2preHI, titre_post_vac_1 = L2PostVHI, titre_post_vac_2 = L2PostSHI) %>% 
        pivot_longer(c(titre_pre_vac, titre_post_vac_1, titre_post_vac_2),
        names_to = "sample_time", values_to = "titre_value") %>%
        mutate(sample_time = case_when(
            sample_time == "titre_pre_vac"~0, 
            sample_time == "titre_post_vac_1"~30,
            sample_time == "titre_post_vac_2"~180))

    save(hcw_pre_sero_data, file = here::here("data", "hcw_pre", "sero_data.RDS"))
}

makeclass <- function() { 
    import_hcw_pre()

    load(here::here("data", "hcw_pre", "part_data.RDS"))
    load(here::here("data", "hcw_pre", "sero_data.RDS"))
    load(here::here("data", "hcw_pre", "vac_data.RDS"))

    hcwpre <- make_study("HCW_pilot",
        "hcw_pre",
        part_data = hcw_pre_part_data,
        sero_data = hcw_pre_sero_data,
        vac_data = hcw_pre_vac_data)
    
    save(hcwpre, file = here::here("data", "hcwpre_data.RDS"))
}