library(tidyverse)
library(plyr)
library(data.table)
library(ggplot2)
library(foreach)
library(doParallel)
# Must install these two
library(serosolver)
library(rcppfunchcw)

library(devtools)
library(here)


#install_github("https://github.com/dchodge/serosolver", ref = "custom_ab")
#install(here::here("rcppfunchcw"))

devtools::load_all() # hcwpre
load(here::here("data", "hanam_data.RDS")) # hcwpre
load(file = here::here("data", paste0("modelinfo_cross_", hanam$study_name_short, ".RDS"))) # all_models_hcw_pre

setup_run_serosolver(
        all_models_hanam_cross[[1]],
        chains = 4,
        iterations = 100000,
        pt = TRUE,
        filename = "dis_dep",
        cross_sectional = TRUE,
        continue_run = FALSE)
