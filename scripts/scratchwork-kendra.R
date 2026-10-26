#Q1 and Q2 draft

library(tidyverse)
load('data/biomarker-clean.RData')

protein_cols <- setdiff(names(biomarker_clean), c('group', 'ados'))

## Q1: WHY LOG-TRANSFORM?

var_names <- read_csv('data/biomarker-raw.csv',
                      col_names = F, n_max = 2, col_select = -(1:2)) %>%
  t() %>% as_tibble() %>%
  rename(name = V1, abbreviation = V2) %>%
  na.omit()

biomarker_raw <- read_csv('data/biomarker-raw.csv',
                          skip = 2, col_select = -2L,
                          col_names = c('group', 'empty',
                                        pull(var_names, abbreviation), 'ados'),
                          na = c('-', '')) %>%
  filter(!is.na(group))

# random sample of 9 proteins
set.seed(101422)
sample_proteins <- sample(protein_cols, 9)

# beofre: raw levels
biomarker_raw %>%
  select(all_of(sample_proteins)) %>%
  pivot_longer(everything(), names_to = 'protein', values_to = 'level') %>%
  ggplot(aes(x = level)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ protein, scales = 'free') +
  labs(title = 'Raw protein levels (random sample of 9)',
       x = 'Raw level', y = 'Count')

# after: same proteins from biomarker_clean (log10 + scaled)
biomarker_clean %>%
  select(all_of(sample_proteins)) %>%
  pivot_longer(everything(), names_to = 'protein', values_to = 'level') %>%
  ggplot(aes(x = level)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ protein, scales = 'free') +
  labs(title = 'Cleaned protein levels (log10, scaled, trimmed)',
       x = 'Standardized log10 level', y = 'Count')

# skewness across all proteins: raw vs. cleaned
skew <- function(x) {
  x <- x[!is.na(x)]
  mean((x - mean(x))^3) / sd(x)^3
}

skew_tbl <- tibble(
  protein = protein_cols,
  raw   = map_dbl(protein_cols, ~ skew(biomarker_raw[[.x]])),
  clean = map_dbl(protein_cols, ~ skew(biomarker_clean[[.x]]))
)

skew_tbl %>%
  summarize(median_skew_raw   = median(raw),
            median_skew_clean = median(clean),
            pct_skewed_raw    = mean(raw > 1) * 100,
            pct_skewed_clean  = mean(clean > 1) * 100)

skew_tbl %>%
  pivot_longer(c(raw, clean), names_to = 'version', values_to = 'skewness') %>%
  mutate(version = factor(version, levels = c('raw', 'clean'))) %>%
  ggplot(aes(x = skewness, fill = version)) +
  geom_histogram(bins = 50, alpha = 0.6, position = 'identity') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(title = 'Skewness of all proteins: raw vs. cleaned',
       x = 'Skewness', y = 'Number of proteins')

## Q2: OUTLIERS

outlier_long <- biomarker_clean %>%
  mutate(subject_id = row_number()) %>%
  select(-ados) %>%
  pivot_longer(-c(subject_id, group),
               names_to = 'protein', values_to = 'z') %>%
  mutate(outlier = abs(z) == 3)

# overall: how many values did trimming change?
outlier_long %>%
  summarize(n_outlying   = sum(outlier),
            pct_outlying = mean(outlier) * 100)

# number of outlying values PER SUBJECT
outlier_counts <- outlier_long %>%
  group_by(subject_id, group) %>%
  summarize(n_out  = sum(outlier),
            n_high = sum(z == 3),
            n_low  = sum(z == -3),
            .groups = 'drop') %>%
  arrange(desc(n_out))

# subjects with the most outlying values
outlier_counts %>% head(15)

# distribution of counts by group
outlier_counts %>%
  group_by(group) %>%
  summarize(n_subjects = n(),
            mean_out   = mean(n_out),
            median_out = median(n_out),
            max_out    = max(n_out))

ggplot(outlier_counts, aes(x = n_out, fill = group)) +
  geom_histogram(bins = 40, alpha = 0.7, position = 'identity') +
  labs(title = 'Outlying values (|z| > 3) per subject',
       x = 'Number of outlying protein values', y = 'Number of subjects')

ggplot(outlier_counts, aes(x = group, y = n_out)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(title = 'Outlying values per subject by group',
       x = 'Group', y = 'Number of outlying values')

# define an outlier subject with the boxplot rule on the counts
cutoff <- quantile(outlier_counts$n_out, 0.75) + 1.5 * IQR(outlier_counts$n_out)
cutoff

outlier_counts <- outlier_counts %>%
  mutate(outlier_subject = n_out > cutoff)

outlier_counts %>% filter(outlier_subject)

# are outlier subjects more frequent in one group?
outlier_table <- table(outlier_counts$group, outlier_counts$outlier_subject)
outlier_table
prop.table(outlier_table, margin = 1)
fisher.test(outlier_table)

# do counts differ by group overall?
wilcox.test(n_out ~ group, data = outlier_counts)