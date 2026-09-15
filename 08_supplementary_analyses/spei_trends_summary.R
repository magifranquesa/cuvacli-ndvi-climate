# ==============================================================================
# TENDENCIAS DEL SPEI 1984-2023: RESUMEN POR PARQUE
# ==============================================================================
# Resume la tabla de tendencias Mann-Kendall del SPEI por parque, mes y escala
# y cuenta cuántas son significativas (p < 0,05) en cada parque.
#
# Entrada: resultados/clima/tendencias_mensuales_spei_por_parque_1984_2023.csv
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

suppressMessages({library(dplyr); library(readr)})
options(width = 130)
d <- read_csv(paste0(ROOT, "/articulo_2/resultados/clima/tendencias_mensuales_spei_por_parque_1984_2023.csv"),
              show_col_types = FALSE)
cat("Tendencias SPEI significativas (p < 0,05), 1984-2023\n")
cat("(cada parque: 5 escalas x 12 meses = 60 tests)\n\n")
cat(sprintf("%-16s %7s %8s %7s   %s\n","parque","tests","signif","%","meses con tendencia signif."))
r <- d %>% group_by(Parque) %>%
     summarise(n = n(), s = sum(P.Value < 0.05, na.rm = TRUE),
               meses = paste(sort(unique(Mes[P.Value < 0.05])), collapse=","), .groups="drop") %>%
     arrange(s)
for (i in 1:nrow(r)) cat(sprintf("%-16s %7d %8d %6.1f%%   %s\n", r$Parque[i], r$n[i], r$s[i],
    100*r$s[i]/r$n[i], substr(r$meses[i],1,44)))
cat(sprintf("\n%-16s %7d %8d %6.1f%%\n", "TOTAL", sum(r$n), sum(r$s), 100*sum(r$s)/sum(r$n)))
