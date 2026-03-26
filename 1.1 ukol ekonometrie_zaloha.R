
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

