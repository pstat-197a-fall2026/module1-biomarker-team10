# Emma Nguyen
# Scratchwork for 4a
# libraries

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(ggplot2)

# load data
load('data/biomarker-clean.RData')

# task: simpler panel - fewer proteins that classify as well as in-class panel
# drop the proteins to show metrics barely change 
# reproduce the benchmark in in-class
# eval function

# STEP 1
# Reproduce the benchmark from in-class panel
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

# STEP 2
# EVALUATION FUNCTION
# score any panel with same split + model + metrics from class 
class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

eval_panel <-function(proteins, label, seed = 101422) {
  panel_data <- biomarker_clean |>
    mutate(class = factor(group == "ASD", levels = c(FALSE,TRUE))) |>
    select(class, all_of(proteins))

  set.seed(seed)
  split <- initial_split(panel_data, prop = .80, strata = class)
  
  # log reg on training set 
  fit <- glm(class ~., data = training(split), family = binomial())
  
  # test split 
  test <- testing(split)
  
  # predict test children + score, returns one row per panel 
  test |>
    mutate(pred = predict(fit, newdata = test, type = "response"), # probability of ASD 
           est = factor(pred > 0.5, levels = c(FALSE, TRUE))) |>
    class_metrics(truth = class, estimate = est, pred, event_level = "second") |>
    select(.metric, .estimate) |>
    pivot_wider(names_from = .metric, values_from = .estimate) |> # condense to 1 row 
    mutate(panel_type = label, n_proteins = length(proteins),
           proteins_list = paste(proteins, collapse = "|"))
}

# ------------------------------------------------------------------
# BENCHMARK TABLE
results <- eval_panel(proteins_sstar, label = "in-class")
print(results)

# ------------------------------------------------------------------
# # STEP 3
# Method for simpler panel
# CHOSEN METHOD: top-k from the t-test and RF rankings

# best to worst by t test 
rank_t <- ttests_out |> arrange(p.adj) |> pull(protein)

# best to worst by random forest 
rank_rf <- rf_out$importance |>
  as_tibble(rownames = "protein") |>
  arrange(desc(MeanDecreaseGini)) |>
  pull(protein)

# scoring the top k proteins from t-test and rf ranking 
topk <- map(1:4, function(k) {
  bind_rows(eval_panel(rank_t[1:k], label = paste("t-test top", k)),
            eval_panel(rank_rf[1:k], label = paste("RF top", k)))
  }) |>
  list_rbind()

topk |> select(panel_type, n_proteins, proteins_list, accuracy, roc_auc) |>
        arrange(desc(roc_auc))

# OTHER NOTES:  
# LDA for feature extraction? Combo of existing features to yield a smaller number of features 
# LASSO - automatically performs feature selection, good for irrelevant variables 
# coefficients shrink at different rates
# less useful predictors are shrunk to 0 first  

# BEST SUBSET METHOD 
#best_subset <- map(1:5, ~ combn(proteins_sstar, .x, simplify = FALSE)) |>
#  flatten()

#bs_results <- map_dfr(best_subset, 
#                          ~ eval_panel(.x, label = paste(.x, collapse = ", ")))

#bs_table <- bs_results |>
#  group_by(n_proteins) |>
#  slice_max(roc_auc, n = 1, with_ties = FALSE) |>
#  ungroup() |>
#  arrange(desc(roc_auc))

#print(bs_table)

# STEP 4
# STABILITY CHECK - multiple splits (chose 100 random splits)

# includes in-class and top k from each ranking 
candidates <- c(list(`in-class` = proteins_sstar),
                set_names(map(1:4, ~ rank_rf[1:.x]), paste("RF top", 1:4)),
                set_names(map(1:4, ~ rank_t[1:.x]), paste("t-test top", 1:4)))

# run eval_panel() on each and then calculate the averages
stability <- imap(candidates, function(proteins_list, label) {
  map(1:100, ~ eval_panel(proteins_list, label = label, seed = .x)) |>
    list_rbind()}) |>
  list_rbind() |>
  group_by(panel_type, n_proteins, proteins_list) |>
  summarize(mean_auc = mean(roc_auc), 
            sd_auc = sd(roc_auc),
            mean_accuracy = mean(accuracy),
            mean_sensitivity = mean(sensitivity),
            mean_specificity = mean(specificity),
            .groups = "drop") |>
  arrange(desc(mean_auc))

stability

# STEP 5
# CHOOSING THE SIMPLER PANEL
# comparing mean AUC + mean accuracy within 0.025 of the in-class panel

benchmark <- stability |> filter(panel_type == "in-class")
benchmark

simpler_panel <- stability |> filter(panel_type != "in-class",
                                     mean_auc >= benchmark$mean_auc - 0.025,
                                     mean_accuracy >= benchmark$mean_accuracy - 0.025) |>
  # ones with least proteins 
  slice_min(n_proteins, n = 1) |>
  # highest auc of the results from above line 
  slice_max(mean_auc, n = 1)

simpler_panel

# STEP 6
# Comparison table of original in-class panel vs simpler panel 

comparison <- stability |>
  filter(panel_type %in% c("in-class", simpler_panel$panel_type)) |>
  select(panel_type, proteins_list, n_proteins,
  mean_auc, sd_auc, mean_accuracy,
  mean_sensitivity, mean_specificity) |>
  arrange(desc(n_proteins))

comparison

# VISUALIZATION 
comp_long <- comparison |> 
  pivot_longer(starts_with("mean_"), names_to = "metric", values_to = "value")

ggplot(comp_long, aes(x = metric, y = value, fill = panel_type)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("lightblue", "steelblue")) +
  labs(title = "In-class vs simpler panel", x = "Metric", y = "Mean value", fill = "Panel")

results
topk
stability
simpler_panel
comparison

write_csv(comparison, "q4a-simpler-panel.csv")







