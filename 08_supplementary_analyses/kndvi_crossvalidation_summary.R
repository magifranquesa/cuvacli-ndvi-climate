# ==============================================================================
# VALIDACIÓN CRUZADA NDVI - kNDVI: RESUMEN
# ==============================================================================
# A partir de xval.rds (kndvi_crossvalidation.R): las combinaciones positivas
# en que el NDVI da al menos un 20 % de superficie significativa en media
# anual, con sus estadísticos de acuerdo, y cuántas de las 48 combinaciones
# positivas dan menos superficie con kNDVI que con NDVI.
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")
OUT <- file.path(ROOT, "articulo_2/resultados_v2/analisis_complementarios"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

options(width = 130)
res <- readRDS(file.path(OUT, "xval.rds"))
x <- res[res$signo == "pos" & res$med_ndvi >= 20, ]
x <- x[order(-x$med_ndvi), ]
cat("Las 12 combinaciones seleccionadas (positivas, media anual NDVI >= 20 %):\n\n")
cat(sprintf("%-14s %-5s %8s %8s %9s %8s %7s %7s\n",
    "parque","var","NDVI %","kNDVI %","perfil r","pixel r","kappa","esc %"))
for (i in 1:nrow(x)) cat(sprintf("%-14s %-5s %8.1f %8.1f %9.3f %8.3f %7.3f %7.1f\n",
  x$parque[i], x$variable[i], x$med_ndvi[i], x$med_kndvi[i],
  x$perfil_r[i], x$px_r[i], x$kappa[i], x$escala_ac[i]))
cat(sprintf("\n%-20s %8.1f %8.1f\n", "MEDIA de las 12:", mean(x$med_ndvi), mean(x$med_kndvi)))
cat(sprintf("%-20s %26.3f %8.3f %7.3f %7.1f\n", "MEDIANA de las 12:",
    median(x$perfil_r), median(x$px_r), median(x$kappa), median(x$escala_ac)))
cat(sprintf("%-20s %26s %8s %7s\n", "  rango:",
    sprintf("%.2f-%.2f", min(x$perfil_r), max(x$perfil_r)),
    sprintf("%.2f-%.2f", min(x$px_r), max(x$px_r)),
    sprintf("%.2f-%.2f", min(x$kappa), max(x$kappa))))
cat(sprintf("\nCuantas de las 48 positivas dan MENOS superficie con kNDVI: %d\n",
    sum(res$signo=="pos" & res$med_kndvi < res$med_ndvi)))
