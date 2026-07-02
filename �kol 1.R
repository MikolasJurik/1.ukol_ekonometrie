
# 1.1 SBĚR DAT A KONVERZE
# Balíček pro data z Yahoo Finance a FRED
library(quantmod)
library(tidyverse)

# 50 TOP S&P 500
tkrs <- c("AAPL","MSFT","NVDA","AMZN","META","GOOGL","BRK-B","LLY","AVGO","JPM",
          "TSLA","UNH","XOM","V","JNJ","MA","PG","HD","COST","MRK",
          "ABBV","CRM","BAC","NFLX","AMD","CVX","PEP","KO","WMT","TMO",
          "MCD","CSCO","INTC","ABT","INTU","WFC","CMCSA","IBM","AMGN","CAT",
          "UNP","TXN","PM","BA","COP","HON","GE","SPG","DIS","ORCL","NKE")

# 1. Akcie  (>= 60M) → r_j
# do.call(merge, ...) = spojení extrahovaných sloupců do 1 matice
# lapply = cyklus přes všechny tickery
ceny_j <- do.call(merge, lapply(tkrs, function(x) {
  # getSymbols = stáhnutí měsíčních dat; auto.assign=FALSE = přímý výstup; Cl = izolace zavírací ceny($P_t$)
  Cl(getSymbols(x, auto.assign=FALSE, from="2018-01-01", periodicity="monthly"))
}))
# ROC = $P_t → výnosy $r_j$; [-1,] = odstranění 1. řádku (`NA`)
r_j <- ROC(ceny_j, type="discrete")[-1,]

# 2. Trh stáhnutí S&P 500 → $r_m$
# auto.assign=FALSE = řeší syntax error s "^GSPC" jelikož stahuje měsíční data trhu přímo do proměnné
trh_ceny <- getSymbols("^GSPC", auto.assign=FALSE, from="2018-01-01", periodicity="monthly")
# zavírací cena → výnosy $r_m$ → smazání `NA`
r_m <- ROC(Cl(trh_ceny), type="discrete")[-1,]

# 3. stáhnutí $r_f$ z FRED (3M T-Bill) → měsíční $r_f$
# stáhnutí roční sazby (% p.a.)
rf_data <- getSymbols("TB3MS", src="FRED", auto.assign=FALSE)
# roční p.a. → efektivní měsíční $r_f$
r_f <- (1 + rf_data/100)^(1/12) - 1

# merge = sjednocení dle data; na.omit = smazání chybějících hodnot
data_capm <- na.omit(merge(r_j, r_m, r_f))
# úhledné názvy sloupců → příprava na OLS
colnames(data_capm) <- c(tkrs, "r_m", "r_f")




# 1.2. VIZUALIZACE
# xts objekt → data.frame (ggplot2 vyžaduje běžný data.frame)
df_capm <- data.frame(Datum = index(data_capm), coredata(data_capm))

# 1. Časové řady: Trh ($r_m$) vs. Risk-free ($r_f$)
# pivot_longer → převede sloupce do dlouhého formátu pro kreslení čar
df_ts <- df_capm %>% 
  select(Datum, r_m, r_f) %>% 
  pivot_longer(-Datum, names_to="Metrika", values_to="Hodnota")

graf_ts <- ggplot(df_ts, aes(x=Datum, y=Hodnota, color=Metrika)) +
  geom_line(linewidth=0.7) +
  theme_minimal() +
  labs(title="Vývoj: Trh (r_m) vs. Bezriziková sazba (r_f)", y="Měsíční míra")

print(graf_ts)


# 2. Akcie: Sdružený boxplot 51 výnosů ($r_j$)
# Odstraníme trh a rf → zbyde 51 akcií → dlouhý formát
df_akcie <- df_capm %>% 
  select(-r_m, -r_f) %>% 
  pivot_longer(-Datum, names_to="Akcie", values_to="Vynos")

# x=reorder → seřadí akcie zleva doprava dle mediánu
graf_box <- ggplot(df_akcie, aes(x=reorder(Akcie, Vynos, FUN=median), y=Vynos)) +
  geom_boxplot(fill="steelblue", outlier.size=0.5, alpha=0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=7)) +
  labs(title="Rozdělení výnosů 51 akcií (dle mediánu)", x="Akciový titul", y="Měsíční výnos (r_j)")

print(graf_box)

# ===========================================================
# 1.3 ODHAD SYSTEMATICKEHO RIZIKA
library(broom)

# 1. Výpočet nadvýnosů (Excess returns)
# transformuje data na čisté vstupy pro OLS regresi →na levé straně rovnice jsou nyní nadvýnosy akcií (ex_AAPL atd.) a na pravé straně nadvýnos trhu (ex_rm).
df_excess <- df_capm %>%
  # across(-c(...)) -> aplikuj na všechny sloupce KROMĚ výjimek
  mutate(across(-c(Datum, r_m, r_f), ~ .x - r_f, .names = "ex_{.col}"),
         ex_rm = r_m - r_f) %>%
  select(Datum, starts_with("ex_"))

# Matice → Dlouhý formát pro hromadné operace
df_long <- df_excess %>%
  pivot_longer(cols = starts_with("ex_") & !matches("ex_rm"), 
               names_to = "Akcie", values_to = "ex_rj") %>%
  mutate(Akcie = str_remove(Akcie, "ex_")) # Smazání prefixu "ex_"


# 2. Hromadný OLS odhad bety pro 51 akcií
# nest → seskupí data dle akcie
# map + lm → 50x regrese
# tidy → extrakce koeficientů do tabulky
vysledky_capm <- df_long %>%
  nest(data = c(Datum, ex_rj, ex_rm)) %>%
  mutate(model = map(data, ~ lm(ex_rj ~ ex_rm, data = .x)),
         odhad = map(model, tidy)) %>%
  unnest(odhad) %>%
  select(Akcie, term, estimate, p.value)

# Výpis jen beta koeficientů
bety <- vysledky_capm %>% 
  filter(term == "ex_rm") %>% 
  rename(Beta = estimate) %>%
  arrange(desc(Beta))
print(bety, n = Inf)


# 3. Vizualizace regresních křivek (51 grafů v jednom) 
graf_krivky <- ggplot(df_long, aes(x = ex_rm, y = ex_rj)) +
  geom_point(alpha = 0.4, size = 0.5, color="darkgray") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.6) +
  facet_wrap(~ Akcie, ncol = 10) + # 51 malých grafů v jedné mřížce
  theme_minimal() +
  theme(strip.text = element_text(size = 7, face="bold"),
        axis.text = element_text(size = 5)) +
  labs(title = "Regresní křivky vybraných 51 akciích", 
       x = "Nadvýnos trhu (rm - rf)", 
       y = "Nadvýnos akcie (rj - rf)")

print(graf_krivky)

# ===============================================================
# 1.4 TEST H0: ALFA = 0

# 1. Extrakce urovńových konstant (Intercept) a p-hodnot
alfy_test <- vysledky_capm %>%
  filter(term == "(Intercept)") %>%
  rename(Alfa = estimate, P_hodnota = p.value) %>%
  # Vytvoření sloupce s výsledkem testu (Hladina významnosti = 0.05)
  mutate(Zaver_Testu = ifelse(P_hodnota < 0.05, 
                              "Zamítáme H0 (Alfa ≠ 0)", 
                              "Nezamítáme H0 (Alfa = 0)"))

# 2.Souhrnná tabulka testu H0
summary_alfy <- alfy_test %>%
  group_by(Zaver_Testu) %>%
  summarise(Pocet_Akcii = n(),
            Podil = paste0(round(n() / 51 * 100, 1), " %"))

print("--- Souhrn testování nulovosti konstanty ---")
print(summary_alfy)

# Výpis akcií, pro které CAPM selhává (alfa != 0)
akcie_s_alfou <- alfy_test %>% filter(P_hodnota < 0.05) %>% pull(Akcie)
cat("\nAkcie se statisticky významnou alfou:", paste(akcie_s_alfou, collapse=", "), "\n")

# 3. Kompaktní přehled (Bar chart)
# Zobrazí p-hodnoty všech akcií s vyznačenou 5% hranicí
graf_alfy <- ggplot(alfy_test, aes(x = reorder(Akcie, P_hodnota), y = P_hodnota, fill = Zaver_Testu)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed", linewidth = 1) +
  coord_flip() + # Otočení pro lepší čitelnost 51 názvů
  scale_fill_manual(values = c("Nezamítáme H0 (Alfa = 0)" = "steelblue", 
                               "Zamítáme H0 (Alfa ≠ 0)" = "tomato")) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5),
        legend.position = "bottom") +
  labs(title = "Test významnosti alfy pro 51 akcií",
       subtitle = "Červená přerušovaná čára = 5% hladina významnosti (0.05)",
       x = "Akciový titul",
       y = "P-hodnota",
       fill = "Výsledek testu: ")

print(graf_alfy)


# PŘÍKLAD 2: TESTOVÁNÍ PLATNOSTI CAPM

library(tidyverse)
library(broom)
library(ggplot2)

# 2.1 SESTAVENÍ PORTFOLIÍ DLE BETA KOEFICIENTŮ

# beta koeficienty z cvičení 1
bety_p2 <- vysledky_capm %>%
  filter(term == "ex_rm") %>%
  select(Akcie, Beta = estimate)

# extrakce rozptylu reziduí pro Fama-MacBeth test
modely_p2 <- df_long %>%
  nest(data = c(Datum, ex_rj, ex_rm)) %>%
  mutate(
    model  = map(data, ~ lm(ex_rj ~ ex_rm, data = .x)),
    sigma2 = map_dbl(model, ~ var(residuals(.x)))   # σ2_εj
  ) %>%
  select(Akcie, sigma2)

# přiřazení akcií do 10 portfolií dle beta (sestupně: P10 = nejvyšší)
params_p2 <- bety_p2 %>%
  left_join(modely_p2, by = "Akcie") %>%
  arrange(Beta) %>%
  mutate(
    Portfolio = paste0("P", sprintf("%02d", ntile(Beta, 10)))
  )

cat("Rozdělení akcií do portfolií:\n")
print(params_p2 %>% group_by(Portfolio) %>%
        summarise(n = n(), Beta_min = min(Beta), Beta_max = max(Beta)))

# výpočet průměrných výnosů portfolií (equal-weighted)
# měsíční výnosy portfolia = průměr nadvýnosů akcií v portfoliu
df_portfolia <- df_long %>%
  left_join(params_p2 %>% select(Akcie, Portfolio), by = "Akcie") %>%
  group_by(Datum, Portfolio) %>%
  summarise(rp = mean(ex_rj), .groups = "drop")

# agregace přes celé období: průměrný měsíční nadvýnos portfolia
rp_mean <- df_portfolia %>%
  group_by(Portfolio) %>%
  summarise(rp_avg = mean(rp), .groups = "drop")

# beta portfolia = aritmetický průměr individuálních bet
beta_p <- params_p2 %>%
  group_by(Portfolio) %>%
  summarise(
    Beta_p   = mean(Beta),       # průměrná beta portfolia
    sigma2_p = mean(sigma2),     # průměrný σ2 pro referenci
    n_akcii  = n(),
    .groups  = "drop"
  )

# finální datová sada pro SML
data_sml <- beta_p %>%
  left_join(rp_mean, by = "Portfolio")

cat("\nPřehled portfolií (beta × průměrný výnos):\n")
print(data_sml)

# 2.2 SECURITY MARKET LINE (SML)

# ols: r̄_pi = γ0 + γ1 * β_pi + ε_pi
# (r̄_pi je nadvýnos portfolia — rf již odečteno v Příkladu 1)
model_sml <- lm(rp_avg ~ Beta_p, data = data_sml)
summ_sml  <- summary(model_sml)
tidy_sml  <- tidy(model_sml)

cat("\n=== Odhad Security Market Line ===\n")
print(tidy_sml)
cat(sprintf("R2 = %.4f\n", summ_sml$r.squared))

# průměrná tržní riziková prémie — benchmark γ1
trzni_premie <- mean(df_long$ex_rm, na.rm = TRUE)
gamma0 <- coef(model_sml)[1]
gamma1 <- coef(model_sml)[2]
p_g0   <- tidy_sml %>% filter(term == "(Intercept)") %>% pull(p.value)
p_g1   <- tidy_sml %>% filter(term == "Beta_p")      %>% pull(p.value)

cat(sprintf("\nBenchmark γ1 (průměrná tržní prémie): %.6f\n", trzni_premie))
cat(sprintf("Odhadnuté γ1:                          %.6f\n", gamma1))
cat(sprintf("\nTest γ0 = 0: p-value = %.4f → %s\n",
            p_g0, ifelse(p_g0 < 0.05, "ZAMÍTÁME H0", "NEZAMÍTÁME H0")))
cat(sprintf("Test γ1 = tržní prémie: γ1 = %.6f vs. benchmark = %.6f\n",
            gamma1, trzni_premie))

# SML vizualizace
graf_sml <- ggplot(data_sml, aes(x = Beta_p, y = rp_avg)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(aes(label = Portfolio), vjust = -0.7, size = 3) +
  # Odhadnutá SML (OLS)
  geom_smooth(method = "lm", se = TRUE,
              color = "tomato", linewidth = 0.8,
              fill = "tomato", alpha = 0.1) +
  # Teoretická SML dle CAPM (γ0 = 0, γ1 = tržní prémie)
  geom_abline(intercept = 0, slope = trzni_premie,
              linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
  theme_minimal() +
  labs(
    title    = "Security Market Line (SML) — 10 portfolií dle beta",
    subtitle = sprintf("y1 = %.5f (p = %.3f)    y2 = %.5f (p = %.3f)",
                       gamma0, p_g0, gamma1, p_g1),
    x       = "Beta koeficient portfolia (β̂p)",
    y       = "Průměrný měsíční nadvýnos portfolia (r̄p − rf)",
    caption = "Červená = odhadnutá SML  |  Zelená přerušovaná = teoretická SML (CAPM: γ0 = 0)"
  )
print(graf_sml)

# 2.3 TEST FAMA-MACBETH / GIBBONS (1982)
# Specifikace cross-sekční regrese na individuálních akciích:
#   r̄j = δ1 + δ2·βj + δ3·βj2 + δ4·σ2εj + ε

# průměrné nadvýnosy individuálních akcií (závisle proměnná)
rj_mean <- df_long %>%
  group_by(Akcie) %>%
  summarise(rj_avg = mean(ex_rj), .groups = "drop")

# Datová sada: beta, beta2, sigma2
data_fm <- params_p2 %>%
  left_join(rj_mean, by = "Akcie") %>%
  mutate(Beta2 = Beta^2)

# Modely porovnání

# Model 1: Základní SML na individuálních akciích (lineární)
model_fm1 <- lm(rj_avg ~ Beta, data = data_fm)

# Model 2: Plná Fama-MacBeth specifikace
model_fm2 <- lm(rj_avg ~ Beta + Beta2 + sigma2, data = data_fm)

summ_fm2 <- summary(model_fm2)
tidy_fm2 <- tidy(model_fm2) %>%
  mutate(
    stars = case_when(
      p.value < 0.01 ~ "***",
      p.value < 0.05 ~ "**",
      p.value < 0.10 ~ "*",
      TRUE           ~ ""
    )
  )

cat("\n=== Fama-MacBeth / Gibbons test ===\n")
print(tidy_fm2 %>% select(term, estimate, std.error, p.value, stars))
cat(sprintf("R2 = %.4f\n", summ_fm2$r.squared))

# výsledky testů klíčových hypotéz
delta3_p <- summ_fm2$coefficients["Beta2",  4]
delta4_p <- summ_fm2$coefficients["sigma2", 4]

cat(sprintf("\n[TEST LINEARITY]        H0: δ3 = 0  →  p-value = %.4f  →  %s\n",
            delta3_p, ifelse(delta3_p < 0.05, "ZAMÍTÁME H0 ✗", "NEZAMÍTÁME H0 ✓")))
cat(sprintf("[TEST SYS. RIZIKA]      H0: δ4 = 0  →  p-value = %.4f  →  %s\n",
            delta4_p, ifelse(delta4_p < 0.05, "ZAMÍTÁME H0 ✗", "NEZAMÍTÁME H0 ✓")))

# vizualizace: koeficienty FM testu s intervalem spolehlivosti
graf_fm <- tidy_fm2 %>%
  mutate(
    term_label = recode(term,
                        "(Intercept)" = "δ1 (intercept)",
                        "Beta"        = "δ2 (β)",
                        "Beta2"       = "δ3 (β2)",
                        "sigma2"      = "δ4 (σ2ε)"
    )
  ) %>%
  ggplot(aes(x = term_label, y = estimate,
             ymin = estimate - 2 * std.error,
             ymax = estimate + 2 * std.error,
             color = p.value < 0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(linewidth = 0.8) +
  geom_text(aes(label = stars), vjust = -1.2, size = 4, color = "black") +
  coord_flip() +
  scale_color_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue"),
                     labels = c("TRUE" = "Sig. (p < 0.05)", "FALSE" = "Nesig."),
                     name = "") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    title    = "Fama-MacBeth test: odhady koeficientů (±2 SE)",
    subtitle = "H0: δ3 = 0 (linearita)  a  δ4 = 0 (pouze systematické riziko)",
    x = NULL,
    y = "Odhad koeficientu"
  )
print(graf_fm)

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


# ===========================================================
# 3.2 Security Market Line - dva subperiody (pred/po zlomu)
# ===========================================================
# Pouzivame beta_pre a beta_post z break_results
# a odpovidajici prumerne nadvynosy portfolii v kazdem obdobi

# -----------------------------------------------------------
# 2.1 Priprava dat - prumerne nadvynosy portfolii v obou obdobich
# -----------------------------------------------------------

# Nadvynosy v pre-break obdobi
portfolio_pre <- df_long_break %>%
  filter(post_break == 0) %>%
  group_by(Akcie) %>%
  summarise(mean_ex_rj = mean(ex_rj, na.rm = TRUE), .groups = "drop") %>%
  left_join(break_results %>% select(Akcie, beta_pre), by = "Akcie") %>%
  rename(beta = beta_pre) %>%
  mutate(Obdobi = "Pred zlomem")

# Nadvynosy v post-break obdobi
portfolio_post <- df_long_break %>%
  filter(post_break == 1) %>%
  group_by(Akcie) %>%
  summarise(mean_ex_rj = mean(ex_rj, na.rm = TRUE), .groups = "drop") %>%
  left_join(break_results %>% select(Akcie, beta_post), by = "Akcie") %>%
  rename(beta = beta_post) %>%
  mutate(Obdobi = "Po zlomu")

portfolio_sml_32 <- bind_rows(portfolio_pre, portfolio_post)

# -----------------------------------------------------------
# 2.2 Tvorba portfolii serazenych dle bety (min 10 portfolii)
# -----------------------------------------------------------
n_portfolios <- 10

make_portfolios <- function(df) {
  df %>%
    arrange(beta) %>%
    mutate(portfolio = ntile(beta, n_portfolios)) %>%
    group_by(portfolio) %>%
    summarise(
      beta_p    = mean(beta, na.rm = TRUE),
      mean_ex_rp = mean(mean_ex_rj, na.rm = TRUE),
      .groups = "drop"
    )
}

portfolios_pre  <- portfolio_pre  %>% make_portfolios() %>% mutate(Obdobi = "Pred zlomem")
portfolios_post <- portfolio_post %>% make_portfolios() %>% mutate(Obdobi = "Po zlomu")
portfolios_32   <- bind_rows(portfolios_pre, portfolios_post)

# -----------------------------------------------------------
# 2.3 Odhad SML pro kazde obdobi
# -----------------------------------------------------------
sml_models_32 <- portfolios_32 %>%
  group_by(Obdobi) %>%
  nest() %>%
  mutate(
    model   = map(data, ~ lm(mean_ex_rp ~ beta_p, data = .x)),
    tidy    = map(model, tidy),
    glance  = map(model, glance)
  )

sml_coefs_32 <- sml_models_32 %>%
  select(Obdobi, tidy) %>%
  unnest(cols = tidy) %>%
  select(Obdobi, term, estimate, std.error, statistic, p.value) %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "gamma0 (intercept)",
    term == "beta_p"      ~ "gamma1 (slope)",
    TRUE                  ~ term
  ))

print(sml_coefs_32)

# Pozorova trzni rizikova premie v obou obdobich (pro porovnani s gamma1)
market_rp_pre <- df_long_break %>%
  filter(post_break == 0) %>%
  summarise(mean_ex_rm = mean(ex_rm, na.rm = TRUE)) %>%
  pull(mean_ex_rm)

market_rp_post <- df_long_break %>%
  filter(post_break == 1) %>%
  summarise(mean_ex_rm = mean(ex_rm, na.rm = TRUE)) %>%
  pull(mean_ex_rm)

cat("Pozorovana trzni rizikova premie pred zlomem:", round(market_rp_pre, 4), "\n")
cat("Pozorovana trzni rizikova premie po zlomu:  ", round(market_rp_post, 4), "\n")

# -----------------------------------------------------------
# 2.4 Vizualizace SML pro obe obdobi
# -----------------------------------------------------------
graf_sml_32 <- ggplot(portfolios_32, aes(x = beta_p, y = mean_ex_rp, color = Obdobi)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  theme_minimal() +
  labs(
    title = "Ukol 3.2 - Security Market Line: pred a po strukturialnim zlomu",
    x     = "Beta portfolia",
    y     = "Prumerny nadvynos portfolia",
    color = "Obdobi"
  )

print(graf_sml_32)


# ===========================================================
# 3.4 Security Market Line - Bull vs. Bear trh
# ===========================================================
# Pouzivame beta_up a beta_down z asym_results
# Bull obdobi: mesice kdy ex_rm >= 0 -> pouzivame beta_up
# Bear obdobi: mesice kdy ex_rm <  0 -> pouzivame beta_down

# -----------------------------------------------------------
# 4.1 Priprava dat - prumerne nadvynosy ve vzestupnych/klesajicich mesicich
# -----------------------------------------------------------
portfolio_bull <- df_long_asym %>%
  filter(ex_rm_plus > 0) %>%
  group_by(Akcie) %>%
  summarise(mean_ex_rj = mean(ex_rj, na.rm = TRUE), .groups = "drop") %>%
  left_join(asym_results %>% select(Akcie, beta_up), by = "Akcie") %>%
  rename(beta = beta_up) %>%
  mutate(Trh = "Bull (rust)")

portfolio_bear <- df_long_asym %>%
  filter(ex_rm_minus < 0) %>%
  group_by(Akcie) %>%
  summarise(mean_ex_rj = mean(ex_rj, na.rm = TRUE), .groups = "drop") %>%
  left_join(asym_results %>% select(Akcie, beta_down), by = "Akcie") %>%
  rename(beta = beta_down) %>%
  mutate(Trh = "Bear (pokles)")

# -----------------------------------------------------------
# 4.2 Tvorba portfolii serazenych dle bety (10 portfolii)
# -----------------------------------------------------------
n_portfolios <- 10

make_portfolios <- function(df) {
  df %>%
    arrange(beta) %>%
    mutate(portfolio = ntile(beta, n_portfolios)) %>%
    group_by(portfolio) %>%
    summarise(
      beta_p     = mean(beta, na.rm = TRUE),
      mean_ex_rp = mean(mean_ex_rj, na.rm = TRUE),
      .groups    = "drop"
    )
}

portfolios_bull <- portfolio_bull %>% make_portfolios() %>% mutate(Trh = "Bull (rust)")
portfolios_bear <- portfolio_bear %>% make_portfolios() %>% mutate(Trh = "Bear (pokles)")
portfolios_34   <- bind_rows(portfolios_bull, portfolios_bear)

# -----------------------------------------------------------
# 4.3 Odhad SML pro bull a bear obdobi
# -----------------------------------------------------------
sml_models_34 <- portfolios_34 %>%
  group_by(Trh) %>%
  nest() %>%
  mutate(
    model  = map(data, ~ lm(mean_ex_rp ~ beta_p, data = .x)),
    tidy   = map(model, tidy),
    glance = map(model, glance)
  )

sml_coefs_34 <- sml_models_34 %>%
  select(Trh, tidy) %>%
  unnest(cols = tidy) %>%
  select(Trh, term, estimate, std.error, statistic, p.value) %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "gamma0 (intercept)",
    term == "beta_p"      ~ "gamma1 (slope)",
    TRUE                  ~ term
  ))

print(sml_coefs_34)

# Pozorovana trzni rizikova premie v bull/bear obdobich (pro porovnani s gamma1)
market_rp_bull <- df_long_asym %>%
  filter(ex_rm_plus > 0) %>%
  summarise(mean_ex_rm = mean(ex_rm, na.rm = TRUE)) %>%
  pull(mean_ex_rm)

market_rp_bear <- df_long_asym %>%
  filter(ex_rm_minus < 0) %>%
  summarise(mean_ex_rm = mean(ex_rm, na.rm = TRUE)) %>%
  pull(mean_ex_rm)

cat("Pozorovana trzni rizikova premie (bull):", round(market_rp_bull, 4), "\n")
cat("Pozorovana trzni rizikova premie (bear):", round(market_rp_bear, 4), "\n")

# -----------------------------------------------------------
# 4.4 Vizualizace SML pro bull a bear obdobi
# -----------------------------------------------------------
graf_sml_34 <- ggplot(portfolios_34, aes(x = beta_p, y = mean_ex_rp, color = Trh)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  theme_minimal() +
  labs(
    title = "Ukol 3.4 - Security Market Line: Bull vs. Bear trh",
    x     = "Beta portfolia",
    y     = "Prumerny nadvynos portfolia",
    color = "Obdobi trhu"
  )

print(graf_sml_34)

