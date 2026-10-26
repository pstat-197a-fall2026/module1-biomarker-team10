# Emma Nguyen
# Scratchwork for 4a
# libraries

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

# load data
load('~/Desktop/module1-biomarker-team10/data/biomarker-clean.RData')

# task: simpler panel - fewer proteins that classify as well as in-class panel
# drop the proteins to show metrics barely change 
# reproduce the benchmark in in-class
# eval function

# STEP 1: Reproduce the benchmark from in-class panel
# CODE FROM INCLASS-ANALYSIS.R

test_fn <- function(.df) {
  t_test(
    .df,
    formula = level ~ group,
    order = c("ASD", "TD"),
    alternative = "two-sided",
    var.equal = FALSE
  )
}

ttests_out <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = "protein", values_to = "level") %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(
    m = n(),
    hm = log(m) + 1 / (2 * m) - digamma(1),
    rank = row_number(),
    p.adj = m * hm * p_value / rank
  )

proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

# Random forest
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>%
  pull(group) %>%
  factor()

set.seed(101422)
rf_out <- randomForest(
  x = predictors,
  y = response,
  ntree = 1000,
  importance = TRUE
)

proteins_s2 <- rf_out$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

# In-class panel
proteins_sstar <- intersect(proteins_s1, proteins_s2)
proteins_sstar
length(proteins_sstar)

# ------------------------------------------------------------------
# EVALUATION FUNCTION

# ------------------------------------------------------------------
# METHODS 

# best subset of the existing panel ? 

# top k 

# lasso ?