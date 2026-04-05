# ===========================================================
# Ukol 3 
#   3.1 Strukturalni zlom
#   3.3 Asymetricka beta
# ===========================================================

library(tidyverse)
library(broom)
library(purrr)
library(strucchange)   # Chow test
library(car)           # linearHypothesis
library(sandwich)      # robust vcov
library(lmtest)

# ===========================================================
# Strukturalni zlom v charakteristickych krivkach
# -----------------------------------------------------------
# Model:
#   ex_rj = alpha + delta*D + beta*ex_rm + theta*(D*ex_rm) + u
#
# D = 1 v poslednich 24 mesicich
#
# Interpretace:
#   pred zlomem:
#     alpha_pre = alpha
#     beta_pre  = beta
#   po zlomu:
#     alpha_post = alpha + delta
#     beta_post  = beta + theta
#
# Testujeme:
# - joint dummy/interactions test
# - Chowuv strukturalni test
# - Chowuv predpovedni test
# ===========================================================

# -----------------------------------------------------------
# 1.1 Break date = prvni mesic v poslednich 24 mesicich
# -----------------------------------------------------------
all_dates <- sort(unique(df_long$Datum))
break_date <- all_dates[length(all_dates) - 24 + 1]

df_long_break <- df_long %>%
  mutate(post_break = if_else(Datum >= break_date, 1, 0))

# -----------------------------------------------------------
# 1.2 Chowuv predpovedni test - vlastni funkce
# -----------------------------------------------------------
# Odhad prvniho subvzorku, predikce druheho subvzorku
# -----------------------------------------------------------
chow_forecast_test <- function(df){
  df1 <- df %>% filter(post_break == 0)
  df2 <- df %>% filter(post_break == 1)
  
  m1 <- lm(ex_rj ~ ex_rm, data = df1)
  
  RSS1 <- sum(resid(m1)^2)
  pred2 <- predict(m1, newdata = df2)
  RSSf <- sum((df2$ex_rj - pred2)^2)
  
  n1 <- nrow(df1)
  n2 <- nrow(df2)
  k  <- length(coef(m1))
  
  Fstat <- ((RSSf - RSS1) / n2) / (RSS1 / (n1 - k))
  pval  <- 1 - pf(Fstat, df1 = n2, df2 = (n1 - k))
  
  tibble(
    F_chow_forecast = Fstat,
    p_chow_forecast = pval
  )
}

extract_joint_p <- function(obj){
  as.data.frame(obj)$`Pr(>F)`[2]
}

# -----------------------------------------------------------
# 1.3 Odhad modelu pro kazdy akciovy titul
# -----------------------------------------------------------
break_models <- df_long_break %>%
  group_by(Akcie) %>%
  nest() %>%
  mutate(
    bp_index = map_int(data, ~ which(.x$Datum == break_date)[1]),
    
    # Model s dummy promennou a interakci
    model = map(data, ~ lm(ex_rj ~ post_break + ex_rm + post_break:ex_rm, data = .x)),
    tidy = map(model, tidy),
    
    # Spolecny test: zadny zlom v interceptu ani ve sklonu
    test_joint = map(model, ~ linearHypothesis(
      .x,
      c("post_break = 0", "post_break:ex_rm = 0"),
      vcov. = vcovHC(.x, type = "HC1"),
      test = "F"
    )),
    
    # Chowuv test strukturalniho zlomu
    chow_obj = map2(data, bp_index, ~ sctest(ex_rj ~ ex_rm, type = "Chow", point = .y, data = .x)),
    
    # Predikcni Chowuv test
    chow_forecast = map(data, chow_forecast_test)
  )


# -----------------------------------------------------------
# 1.4 Vysledky 3.1.
# -----------------------------------------------------------
break_results <- break_models %>%
  mutate(
    joint_p = map_dbl(test_joint, extract_joint_p),
    F_chow_break = map_dbl(chow_obj, ~ unname(.x$statistic)),
    p_chow_break = map_dbl(chow_obj, ~ .x$p.value)
  ) %>%
  select(Akcie, tidy, joint_p, F_chow_break, p_chow_break, chow_forecast) %>%
  unnest(tidy) %>%
  select(Akcie, term, estimate, p.value, joint_p, F_chow_break, p_chow_break, chow_forecast) %>%
  pivot_wider(names_from = term, values_from = c(estimate, p.value)) %>%
  unnest(chow_forecast) %>%
  mutate(
    alpha_pre  = `estimate_(Intercept)`,
    alpha_post = `estimate_(Intercept)` + estimate_post_break,
    beta_pre   = estimate_ex_rm,
    beta_post  = estimate_ex_rm + `estimate_post_break:ex_rm`,
    break_joint_signif   = joint_p < 0.05,
    chow_break_signif    = p_chow_break < 0.05,
    chow_forecast_signif = p_chow_forecast < 0.05
  ) %>%
  select(
    Akcie,
    alpha_pre, alpha_post,
    beta_pre, beta_post,
    joint_p, break_joint_signif,
    F_chow_break, p_chow_break, chow_break_signif,
    F_chow_forecast, p_chow_forecast, chow_forecast_signif
  ) %>%
  arrange(p_chow_break)

print(break_results, n = Inf)

break_r <- break_results %>%
  mutate(beta_change = beta_post - beta_pre,
         abs_beta_change = abs(beta_change)) %>%
  arrange(desc(abs_beta_change)) %>%
  select(Akcie, beta_pre, beta_post, beta_change, joint_p)
# -----------------------------------------------------------
# 1.5 Shrnuti a vizualizace
# -----------------------------------------------------------
summary_breaks <- break_results %>%
  ungroup() %>%
  summarise(
    pocet_signif_joint = sum(break_joint_signif, na.rm = TRUE),
    pocet_signif_chow_break = sum(chow_break_signif, na.rm = TRUE),
    pocet_signif_chow_forecast = sum(chow_forecast_signif, na.rm = TRUE)
  )

print(summary_breaks)

graf_breaks <- tibble(
  Test = c("Joint test (dummy + slope shift)",
           "Chow structural break",
           "Chow forecast"),
  Pocet = c(summary_breaks$pocet_signif_joint,
            summary_breaks$pocet_signif_chow_break,
            summary_breaks$pocet_signif_chow_forecast)
) %>%
  ggplot(aes(x = Test, y = Pocet)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Ukol 3.1 - Number of stocks with significant break",
       x = "Test type",
       y = "Number of stocks")

print(graf_breaks)

# -----------------------------------------------------------
# 1.5 Analyza vybranych akciovych titulu
# -----------------------------------------------------------

# Vyznamne akciove tituly
vybrane_akcie <- c("ORCL", "COP", "XOM", "META", "MSFT")

df_plot_break <- df_long_break %>%
  filter(Akcie %in% vybrane_akcie) %>%
  mutate(
    Obdobi = if_else(post_break == 1, "Po zlomu", "Pred zlomem")
  )

graf_break_akcie <- ggplot(df_plot_break,
                           aes(x = ex_rm, y = ex_rj, color = Obdobi)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_wrap(~ Akcie, scales = "free") +
  theme_minimal() +
  labs(
    title = "Charakteristicke krivky pred a po strukturialnim zlomu",
    x = "Nadvynos trhu (rm - rf)",
    y = "Nadvynos akcie (rj - rf)",
    color = "Obdobi"
  )

print(graf_break_akcie)

# ===========================================================
# ASYMETRICKA BETA
# -----------------------------------------------------------
# Model:
#   ex_rj = alpha + beta_up * ex_rm_plus + beta_down * ex_rm_minus + u
#
# kde:
#   ex_rm_plus  = max(ex_rm, 0)
#   ex_rm_minus = min(ex_rm, 0)
#
# Test:
#   H0: beta_up = beta_down
# ===========================================================

df_long_asym <- df_long %>%
  mutate(
    ex_rm_plus  = pmax(ex_rm, 0),
    ex_rm_minus = pmin(ex_rm, 0)
  )

extract_equal_p <- function(obj){
  as.data.frame(obj)$`Pr(>F)`[2]
}

asym_models <- df_long_asym %>%
  group_by(Akcie) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(ex_rj ~ ex_rm_plus + ex_rm_minus, data = .x)),
    tidy = map(model, tidy),
    test_equal = map(model, ~ linearHypothesis(
      .x,
      "ex_rm_plus = ex_rm_minus",
      vcov. = vcovHC(.x, type = "HC1"),
      test = "F"
    ))
  )

asym_results <- asym_models %>%
  mutate(
    p_beta_diff = map_dbl(test_equal, extract_equal_p)
  ) %>%
  select(Akcie, tidy, p_beta_diff) %>%
  unnest(tidy) %>%
  select(Akcie, term, estimate, p.value, p_beta_diff) %>%
  pivot_wider(names_from = term, values_from = c(estimate, p.value)) %>%
  mutate(
    beta_up = estimate_ex_rm_plus,
    beta_down = estimate_ex_rm_minus,
    asym_signif = p_beta_diff < 0.05
  ) %>%
  select(Akcie, beta_up, beta_down, p_beta_diff, asym_signif) %>%
  arrange(p_beta_diff)

print(asym_results, n = Inf)


# -----------------------------------------------------------
# 2.1 Vysledky
# -----------------------------------------------------------
summary_asym <- asym_results %>%
  summarise(
    pocet_asym_signif = sum(asym_signif, na.rm = TRUE),
    podil_procent = round(mean(asym_signif, na.rm = TRUE) * 100, 1)
  )

print(summary_asym)

graf_asym <- ggplot(asym_results, aes(x = beta_down, y = beta_up, color = asym_signif)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "tomato"),
                     labels = c("Not significant", "Significant")) +
  labs(title = "Ukol 3.3 - Beta v down market vs up market",
       x = "Beta v down market",
       y = "Beta v up market",
       color = "Test H0: beta_up = beta_down")

print(graf_asym)
