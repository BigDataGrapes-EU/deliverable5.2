if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, hms)

EXPERIMENT <- "PR13-10-23F13" #"PR13-10-07F05"
RUNNING_TIME <- hms(hours = 72) %>% as.numeric() 

read_csv(here::here("data", "raw", "Carignan_k1.csv")) %>%
  select(`code experimentation`, date, `debitCO2 (g)`) %>%
  rename(code_experimentation = `code experimentation`, debitCO2_g = `debitCO2 (g)`) %>%
  group_by(code_experimentation) %>%
  mutate(elapsedtime_s = date - min(date)) %>%
  select(-date) %>%
  write_csv(here::here("data", "interim", "fermentations.csv"))

read_csv(here::here("data", "interim", "fermentations.csv")) %>%
  filter(code_experimentation == EXPERIMENT & elapsedtime_s <= RUNNING_TIME) %>%
  write_csv(here::here("data", "interim", "running_experiment.csv"))

read_csv(here::here("data", "interim", "fermentations.csv")) %>%
  filter(code_experimentation == EXPERIMENT & elapsedtime_s >= RUNNING_TIME) %>%
  write_csv(here::here("data", "interim", "experiment_future.csv"))

read_csv(here::here("data", "interim", "fermentations.csv")) %>%
  filter(code_experimentation != EXPERIMENT) %>%
  write_csv(here::here("data", "interim", "past_experiments.csv"))
