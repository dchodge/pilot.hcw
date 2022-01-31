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

taskIdChar <- Sys.getenv("SGE_TASK_ID")
taskIdInteger <- (as.numeric(taskIdChar))

#install_github("https://github.com/dchodge/serosolver", ref = "custom_ab")
#install(here::here("rcppfunchcw"))

devtools::load_all() # hcwpre
load(here::here("data", "hcwpre_data.RDS")) # hcwpre
load(file = here::here("data", paste0("modelinfo_cross_", hcwpre$study_name_short, ".RDS"))) # all_models_hcw_pre

setup_run_serosolver(
        all_models_hcw_pre_cross[[taskIdInteger]],
        chains = 4,
        iterations = 100000,
        pt = TRUE,
        filename = "dis_dep",
        cross_sectional = TRUE,
        continue_run = FALSE)
