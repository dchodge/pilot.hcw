library(tidyverse)
library(plyr)
library(data.table)
library(ggplot2)
library(foreach)
library(doParallel)
library(serosolver)
library(rcppfunchcw)

library(devtools)
library(here)


#install_github("https://github.com/dchodge/serosolver", ref = "custom_ab")
#install(here::here("rcppfunchcw"))

devtools::load_all() # hcwpre
load(here::here("data", "hanam_data.RDS")) # hcwpre
load(file = here::here("data", paste0("modelinfo_cross_", hcwpre$study_name_short, ".RDS"))) # all_models_hcw_pre

setup_run_serosolver(
        all_models_hcw_pre_cross[[2]],
        chains = 4,
        iterations = 1000000,
        pt = TRUE,
        filename = "dis_dep",
        cross_sectional = FALSE,
        continue_run = FALSE)