library(tidyverse)
library(plyr)
library(data.table)
library(ggplot2)
library(foreach)
library(doParallel)
# Must install these two
library(serosolver)
library(rcppfunchcw)

taskIdChar <- Sys.getenv("SGE_TASK_ID")
taskIdInteger <- (as.numeric(taskIdChar))

load(here::here("data", "hcwpre_data.RDS")) # hcwpre
load(file = here::here("data", paste0("modelinfo_", hcwpre$study_name_short, ".RDS"))) # all_models_hcw_pre

1:length(all_models_hcw_pre) %>% map(~
    setup_run_serosolver_working_local(
        all_models_hcw_pre[[taskIdInteger]],
        chains = 4,
        pt = TRUE,
        filename = "dis_dep",
        cross_sectional = TRUE,
        continue_run = FALSE)
)