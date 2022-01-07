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
<<<<<<< HEAD
load(file = here::here("data", paste0("modelinfo_cross_", hanam$study_name_short, ".RDS"))) # all_models_hcw_pre

setup_run_serosolver(
        all_models_hanam_cross[[1]],
=======
load(file = here::here("data", paste0("modelinfo_cross_", hcwpre$study_name_short, ".RDS"))) # all_models_hcw_pre

setup_run_serosolver(
        all_models_hcw_pre_cross[[1]],
>>>>>>> 418e98c3cfcb569d0d803ab7f898932712d6affa
        chains = 4,
        iterations = 1000000,
        pt = TRUE,
        filename = "dis_dep",
        cross_sectional = TRUE,
<<<<<<< HEAD
        continue_run = FALSE)
=======
        continue_run = FALSE)
>>>>>>> 418e98c3cfcb569d0d803ab7f898932712d6affa
