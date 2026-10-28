# Emma Nguyen
# Scratchwork for 4a
# libraries

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(yardstick)

# load data
load('~/Desktop/module1-biomarker-team10/data/biomarker-clean.RData')

# task: simpler panel - fewer proteins that classify as well as in-class panel
# drop the proteins to show metrics barely change 
# reproduce the benchmark in in-class
# eval function

# STEP 1: Reproduce the benchmark from in-class panel
# CODE FROM INCLASS-ANALYSIS.R

test_fn <- function(.df) {t_test(.df, formula = level ~ group, 
                                 order = c("ASD", "TD"),
                                 alternative = "two-sided",
                                 var.equal = FALSE)}

ttests_out <- biomarker_clean %>% select(-ados) %>%
  pivot_longer(-group, names_to = "protein", values_to = "level") %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(),
         hm = log(m) + 1 / (2 * m) - digamma(1),
         rank = row_number(),
         p.adj = m * hm * p_value / rank)

proteins_s1 <- ttests_out %>% slice_min(p.adj, n = 10) %>% pull(protein)

# rf
predictors <- biomarker_clean %>% select(-c(group, ados))
response <- biomarker_clean %>% pull(group) %>% factor()

set.seed(101422)
rf_out <- randomForest(x = predictors,
                       y = response,
                       ntree = 1000,
                       importance = TRUE)

proteins_s2 <- rf_out$importance %>% as_tibble() %>%
                mutate(protein = rownames(rf_out$importance)) %>%
                slice_max(MeanDecreaseGini, n = 10) %>%
                pull(protein)

# In-class panel
proteins_sstar <- intersect(proteins_s1, proteins_s2)
print(proteins_sstar)
length(proteins_sstar)

# ------------------------------------------------------------------
# EVALUATION FUNCTION - 
class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

eval_panel <-function(proteins, label, data = biomarker_clean, seed = 101422) {
  panel_data <- data |>
    mutate(class = factor(group == "ASD", levels = c(FALSE,TRUE))) |>
    select(class, all_of(proteins))

  set.seed(seed)
  
  #split 
  split <- initial_split(panel_data, prop = .80, strata = class)
  
  # fit
  fit <- glm(class ~., data = training(split), family = binomial())
  
  # test split 
  test <- testing(split)
  
  # test 
  test |>
    mutate(pred = predict(fit, newdata = test, type = "response"),
           est = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
    class_metrics(truth = class, estimate = est, pred, event_level = "second") %>%
    select(.metric, .estimate) %>%
    pivot_wider(names_from = .metric, values_from = .estimate) %>%
    mutate(panel = label, n_proteins = length(proteins),
           proteins = paste(proteins, collapse = ", "), .before = 1)
}

# ------------------------------------------------------------------
# BENCHMARK TABLE
results <- eval_panel(proteins_sstar, label = "in-class")
print(results)

# ------------------------------------------------------------------
# METHODS 

# METHOD 1
# top k from the full t-test and RF rankings

rank_t <- ttests_out |>
  arrange(p.adj) |>
  pull(protein)

rank_rf <- rf_out$importance |>
  as_tibble(rownames = "protein") |>
  arrange(desc(MeanDecreaseGini)) |>
  pull(protein)

topk_results <- map_dfr(1:10, function(k) {
  bind_rows(eval_panel(rank_t[1:k], label = paste("t-test top", k)),
            eval_panel(rank_rf[1:k], label = paste("RF top", k)))})

topk_results |> 
  select(panel, n_proteins, proteins, accuracy, roc_auc) |>
  filter(n_proteins < length(proteins_sstar)) |>
  arrange(desc(roc_auc))

# OTHER NOTES:  
# LDA for feature extraction? Combo of existing features to yield a smaller number of features 
# LASSO - automatically performs feature selection, good for irrelevant variables 
# coefficients shrink at different rates
# less useful predictors are shrunk to 0 first  

# BEST SUBSET METHOD 
subsets <- map(1:5, ~ combn(proteins_sstar, .x, simplify = FALSE)) |>
  flatten()

subset_results <- map_dfr(subsets, 
                          ~ eval_panel(.x, label = paste(.x, collapse = ", ")))

subset_table <- subset_results |>
  group_by(n_proteins) |>
  slice_max(roc_auc, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(roc_auc))

print(subset_table)

# ------------------------------------------------------------------
# STABILITY CHECK - multiple splits 
candidates <- c(list(`in-class` = proteins_sstar),
                set_names(map(2:4, ~ rank_rf[1:.x]), paste("RF top", 2:4)),
                set_names(map(2:4, ~ rank_t[1:.x]), paste("t-test top", 2:4)))

stability <- imap_dfr(candidates, function(proteins, label) {
  map_dfr(1:30, ~ eval_panel(proteins, label = label, seed = .x))}) |>
  group_by(panel, n_proteins) |>
  summarize(mean_auc = mean(roc_auc), sd_auc = sd(roc_auc),
            mean_accuracy = mean(accuracy), sd_accuracy = sd(accuracy),
            mean_sensitivity = mean(sensitivity),
            mean_specificity = mean(specificity),
            .groups = "drop") |>
  arrange(desc(mean_auc))

stability

# ------------------------------------------------------------------
# CHOOSING THE SIMPLER PANEL

benchmark <- stability |> filter(panel == "in-class")
benchmark

simpler_panel <- stability |> filter(panel != "in-class",
                                     mean_auc >= benchmark$mean_auc - 0.02,
                                     mean_accuracy >= benchmark$mean_accuracy - 0.02,
                                     n_proteins > 2) |>
  slice_min(n_proteins, n = 1, with_ties = TRUE) |>
  slice_max(mean_auc, n = 1)

simpler_panel
# ------------------------------------------------------------------
# Comparison table of original panel vs simpler panel 

single_split <- bind_rows(results, topk_results) |>
  select(panel, proteins, split_auc = roc_auc, split_accuracy = accuracy)

comparison <- stability |> 
  filter(panel %in% c("in-class", simpler_panel$panel)) |>
  left_join(single_split, by = "panel") |>
  select(panel, proteins, n_proteins,
         split_auc, split_accuracy,
         mean_auc, sd_auc, mean_accuracy, 
         mean_sensitivity, mean_specificity) |>
  arrange(desc(n_proteins))

comparison

