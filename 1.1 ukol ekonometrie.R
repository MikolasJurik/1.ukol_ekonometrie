# Balíček pro data z Yahoo Finance a FRED
library(quantmod)
library(tidyverse)

# 50 TOP S&P 500
tkrs <- c("AAPL","MSFT","NVDA","AMZN","META","GOOGL","BRK-B","LLY","AVGO","JPM",
          "TSLA","UNH","XOM","V","JNJ","MA","PG","HD","COST","MRK",
          "ABBV","CRM","BAC","NFLX","AMD","CVX","PEP","KO","WMT","TMO",
          "MCD","CSCO","INTC","ABT","INTU","WFC","CMCSA","IBM","AMGN","CAT",
          "UNP","TXN","PM","BA","COP","HON","GE","SPG","PLTR","ORCL")

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

# 3. $r_f$ $\downarrow$ FRED (3M T-Bill) → měsíční $r_f$
# stáhnutí roční sazba (% p.a.)
rf_data <- getSymbols("TB3MS", src="FRED", auto.assign=FALSE)
# roční p.a. → efektivní měsíční $r_f$
r_f <- (1 + rf_data/100)^(1/12) - 1

# merge = sjednocení dle data; na.omit = smazání chybějících hodnot
data_capm <- na.omit(merge(r_j, r_m, r_f))
# úhledné názvy sloupců $\to$ příprava na OLS
colnames(data_capm) <- c(tkrs, "r_m", "r_f")