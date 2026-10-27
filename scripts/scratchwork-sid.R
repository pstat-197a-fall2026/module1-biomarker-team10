library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

load('../data/biomarker-clean.RData')


test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

cat("Running t-tests across all proteins...\n")
ttests_out <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank_ttest = row_number(),
         p.adj = m * hm * p_value / rank_ttest)



predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

cat("Fitting Random Forest (1000 trees)...\n")
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = TRUE)

rf_ranks <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  arrange(desc(MeanDecreaseGini)) %>%
  mutate(rank_rf = row_number())


class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

evaluate_panel <- function(selected_proteins, panel_name) {
  if (length(selected_proteins) == 0) {
    cat(sprintf("\n[%s] No proteins selected!\n", panel_name))
    return(NULL)
  }
  
  df_subset <- biomarker_clean %>%
    select(group, any_of(selected_proteins)) %>%
    mutate(class = factor(group == 'ASD', levels = c(FALSE, TRUE))) %>%
    select(-group)
  

  set.seed(101422)
  split_obj <- initial_split(df_subset, prop = 0.8)
  
  fit <- glm(class ~ ., 
             data = training(split_obj), 
             family = 'binomial')
  
  # Predict and compute test set metrics
  test_data <- testing(split_obj) %>%
    add_predictions(fit, type = 'response') %>%
    mutate(pred_class = factor(pred > 0.5, levels = c(FALSE, TRUE)))
  
  metrics <- test_data %>%
    class_metrics(truth = class,
                  estimate = pred_class,
                  pred,
                  event_level = 'second') %>%
    mutate(panel = panel_name,
           n_proteins = length(selected_proteins),
           proteins = paste(selected_proteins, collapse = ", "))
  
  return(list(metrics = metrics, fit = fit, proteins = selected_proteins))
}


# Top 10 hard intersection(what we did in class)
proteins_ttest_10 <- ttests_out %>% slice_min(p.adj, n = 10) %>% pull(protein)
proteins_rf_10    <- rf_ranks %>% slice_head(n = 10) %>% pull(protein)
panel_baseline    <- intersect(proteins_ttest_10, proteins_rf_10)

res_baseline <- evaluate_panel(panel_baseline, "Baseline (Top 10 Hard Intersect)")


# Upping the amount of protiens to 30
proteins_ttest_30 <- ttests_out %>% slice_min(p.adj, n = 30) %>% pull(protein)
proteins_rf_30    <- rf_ranks %>% slice_head(n = 30) %>% pull(protein)
panel_top30       <- intersect(proteins_ttest_30, proteins_rf_30)

res_top30 <- evaluate_panel(panel_top30, "Part 3b (Top 30 Hard Intersect)")


# FUZZY INTERSECTIONS

ranking_df <- ttests_out %>%
  select(protein, rank_ttest, p_value, p.adj) %>%
  inner_join(rf_ranks %>% select(protein, rank_rf, MeanDecreaseGini), by = "protein") %>%
  mutate(avg_rank = (rank_ttest + rank_rf) / 2) %>%
  arrange(avg_rank)

# Fuzzy Option 1: Average Rank (Top 10 consensus)
# Selects the top 10 proteins with the best combined ranking across both methods
panel_fuzzy_average <- ranking_df %>%
  slice_head(n = 10) %>%
  pull(protein)

res_fuzzy_average <- evaluate_panel(panel_fuzzy_average, "Part 3c (Fuzzy Average Top 10)")

# Fuzzy Option 2: More breathing room(Top 30 consensus)
# Qualify if top 15 in one method AND top 30 in the other
proteins_ttest_15 <- ttests_out %>% slice_min(p.adj, n = 15) %>% pull(protein)
proteins_rf_15    <- rf_ranks %>% slice_head(n = 15) %>% pull(protein)

panel_fuzzy_slack <- union(
  intersect(proteins_ttest_15, proteins_rf_30),
  intersect(proteins_rf_15, proteins_ttest_30)
)

res_fuzzy_slack <- evaluate_panel(panel_fuzzy_slack, "Part 3c (Fuzzy Slack Top 15/30)")



comparison_table <- bind_rows(
  res_baseline$metrics,
  res_top30$metrics,
  res_fuzzy_average$metrics,
  res_fuzzy_slack$metrics
) %>%
  select(panel, n_proteins, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)

print(as.data.frame(comparison_table), row.names = FALSE)

cat("\n--- Selected Proteins by Panel ---\n")
cat("1. Baseline (Top 10 Hard Intersect) [", length(panel_baseline), "proteins]:\n  ", 
    paste(panel_baseline, collapse = ", "), "\n")
cat("2. Part 3b (Top 30 Hard Intersect) [", length(panel_top30), "proteins]:\n  ", 
    paste(panel_top30, collapse = ", "), "\n")
cat("3. Part 3c (Fuzzy Average Top 10) [", length(panel_fuzzy_average), "proteins]:\n  ", 
    paste(panel_fuzzy_average, collapse = ", "), "\n")
cat("4. Part 3c (Fuzzy Slack Top 15/30) [", length(panel_fuzzy_slack), "proteins]:\n  ", 
    paste(panel_fuzzy_slack, collapse = ", "), "\n")

