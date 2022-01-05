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

    hcw_pre_sero_data <- hcw_pre_sero_data %>% filter(pid != "RMH0096") # Remove this person as they have loads of missing data and it effects the WAIC

    save(hcw_pre_sero_data, file = here::here("data", "hcw_pre", "sero_data.RDS"))
}

import_hanam <- function() {
    hanam_data_raw <- read_excel(here::here("inst", "extdata", "hanam_data.xlsx"))

    # get part info
    hanam_part_data <- hanam_data_raw %>%
        select(pid = Subject_ID, dob = DoBS, age = Age, gender = SexS) %>%
        mutate(age = round(age, 0), dob = ymd(dob), yob = year(dob))
    hanam_part_data$gender <- factor(hanam_part_data$gender, levels = c("Male", "Female"))
    save(hanam_part_data, file = here::here("data", "hanam", "part_data.RDS"))

    # 
    hanam_sero_data <- hanam_data_raw %>%
        select(pid = Subject_ID, virus_id = virus, virus_name = Short_Name,
            virus_isol_yr = Year, titre_value = L2titre, d4Diff, d7Diff, d14Diff, d21Diff, d280Diff) %>%
        group_by(pid, virus_name) %>%
        mutate(d0Diff = c(0, rep(NA, n() - 1)), .after = titre_value) %>%
        pivot_longer(c(d0Diff, d4Diff, d7Diff, d14Diff, d21Diff, d280Diff), names_to = "sample_time", values_to = "not") %>%
        na.omit %>%
        mutate(sample_time = case_when(
                sample_time == "d0Diff"~0,
                sample_time == "d4Diff"~4,
                sample_time == "d7Diff"~7,
                sample_time == "d14Diff"~14,
                sample_time == "d21Diff"~21,
                sample_time == "d280Diff"~280)) %>%
        select(!not) %>% select(pid, virus_id, virus_name, virus_isol_yr, sample_time, titre_value)
     save(hanam_sero_data, file = here::here("data", "hanam", "sero_data.RDS"))
}

makeclass_hcwpre <- function() { 
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


makeclass_hanam <- function() { 
    import_hanam()

    load(here::here("data", "hanam", "part_data.RDS"))
    load(here::here("data", "hanam", "sero_data.RDS"))

    hanam <- make_study("Ha Nam study",
        "hanam",
        part_data = hanam_part_data,
        sero_data = hanam_sero_data)
    
    save(hanam, file = here::here("data", "hanam_data.RDS"))
}