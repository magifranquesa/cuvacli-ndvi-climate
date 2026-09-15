# ================================================================================
# UMBRAL DE AÑOS VÁLIDOS Y CONTROL DEL FDR — CORRELACIONES POSITIVAS
# ================================================================================
# Se ejecuta después de los scripts de máximos (02_maxima/). Lee
# maximos_cor_positive/ y escribe maximos_cor_positive_fdr/ con los mismos
# ficheros y nombres de variable, recalculando la significancia en dos pasos:
#
#   A) Umbral de n. Se descartan los píxeles con menos de N_MIN años válidos en
#      el mes (n_obs). Los descartados salen del numerador y del denominador de
#      cualquier porcentaje de superficie posterior.
#
#   B) FDR de Benjamini-Hochberg (q = Q_FDR) sobre los p-valores de cada campo
#      parque-mes, es decir, sobre todos los píxeles válidos de un parque en un
#      mes dado (Wilks, 2016, BAMS 97:2263-2273). El procedimiento sigue
#      controlando el FDR bajo la dependencia positiva que cabe esperar entre
#      píxeles vecinos (Benjamini & Yekutieli, 2001).
#
# El p-valor de cada píxel es el de la escala seleccionada por el script de
# máximos; no se aplica ninguna corrección por la selección de escala.
#
# Codificación de significance en la salida: 2 / -2 significativo positivo /
# negativo tras A y B; 1 / -1 no significativo; NA descartado por el umbral.
#
# Este script y fdr_postproceso_neg.R son idénticos salvo DIR_ENTRADA,
# DIR_SALIDA y SIGNO.
#
# Entrada: correlaciones_v2/maximos_cor_positive/<parque>_max_pos_correlation_<var>_<indice>[_selscales].nc
# Salida : correlaciones_v2/maximos_cor_positive_fdr/  (mismos nombres)
# ================================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(ncdf4)
options(width = 130)

# --------------------------------------------------------------------------------
# CONFIGURACIÓN
# --------------------------------------------------------------------------------
DIR_ENTRADA <- paste0(ROOT, '/articulo_2/correlaciones_v2/maximos_cor_positive/')
DIR_SALIDA  <- paste0(ROOT, '/articulo_2/correlaciones_v2/maximos_cor_positive_fdr/')
SIGNO       <- "pos"

parques <- c("aiguestortes","cabaneros","cabrera","cies","daimiel","donana",
             "guadarrama","monfrague","ordesa","picoseuropa",
             "sierranevada","sierranieves")
indices  <- c("NDVI","KNDVI")
var_list <- c("tmax","tmin","pr","spei")

N_MIN    <- 20      # años válidos mínimos por píxel y mes (A)
Q_FDR    <- 0.05    # nivel del FDR (B)

# --------------------------------------------------------------------------------
cat("salida ->", DIR_SALIDA, "\n")
dir.create(DIR_SALIDA, showWarnings = FALSE, recursive = TRUE)
sufijo <- if (SIGNO == "pos") "max_pos_correlation" else "max_neg_correlation"
signo_sig <- if (SIGNO == "pos") 2 else -2

cat(sprintf("=== FDR sobre correlaciones %s ===\n", toupper(SIGNO)))
cat(sprintf("    n>=%d | q=%.2f\n", N_MIN, Q_FDR))

resumen <- data.frame()

for (parque in parques) for (vi in indices) for (var in var_list) {

  # Dos convenciones de nombre según el script de máximos que los generó:
  #   variables (tmax/tmin/pr) -> {parque}_{sufijo}_{var}_{indice}_selscales.nc
  #   SPEI                     -> {parque}_{sufijo}_spei_{indice}.nc   (sin sufijo)
  cand <- c(file.path(DIR_ENTRADA, sprintf("%s_%s_%s_%s_selscales.nc", parque, sufijo, var, tolower(vi))),
            file.path(DIR_ENTRADA, sprintf("%s_%s_%s_%s.nc",           parque, sufijo, var, tolower(vi))))
  f_in <- cand[file.exists(cand)][1]
  if (is.na(f_in)) { cat("falta:", basename(cand[1]), "\n"); next }

  nc  <- nc_open(f_in)
  cor_max <- ncvar_get(nc, "max_correlation")
  sig_ini <- ncvar_get(nc, "significance")
  p_val   <- ncvar_get(nc, "p_value")
  escala  <- ncvar_get(nc, "scale")
  # Sin n_obs no se puede aplicar el umbral de años: se detiene.
  if (!("n_obs" %in% names(nc$var))) {
    nc_close(nc)
    stop("El fichero ", basename(f_in), " no tiene la variable n_obs, ",
         "asi que NO se puede aplicar el umbral N_MIN = ", N_MIN, " anios.\n",
         "  Regeneralo con los scripts de 02_maxima/, ",
         "o pon N_MIN <- 0 si aceptas explicitamente no filtrar por numero de anios.")
  }
  n_obs <- ncvar_get(nc, "n_obs")
  lon <- nc$dim$lon$vals; lat <- nc$dim$lat$vals
  nc_close(nc)

  # (A) umbral de n
  valido <- !is.na(cor_max) & !is.na(n_obs) & n_obs >= N_MIN
  descartado <- 100 * (1 - sum(valido) / sum(!is.na(cor_max)))

  # (B) FDR por campo parque-mes, y significancia nueva con la misma codificación
  sig_new <- array(NA, dim(cor_max))
  crudo <- fdr <- numeric(12)

  for (m in 1:12) {
    ok <- valido[, , m] & !is.na(p_val[, , m])
    den <- sum(valido[, , m], na.rm = TRUE)
    if (den == 0) next

    q <- rep(NA_real_, length(ok)); dim(q) <- dim(ok)
    q[ok] <- p.adjust(p_val[, , m][ok], method = "BH")

    s <- array(NA, dim(ok))
    positivo <- !is.na(cor_max[, , m]) & cor_max[, , m] >= 0
    s[valido[, , m] &  positivo] <-  1
    s[valido[, , m] & !positivo] <- -1
    s[valido[, , m] &  positivo & !is.na(q) & q <= Q_FDR] <-  2
    s[valido[, , m] & !positivo & !is.na(q) & q <= Q_FDR] <- -2
    sig_new[, , m] <- s

    crudo[m] <- 100 * sum(sig_ini[, , m] == signo_sig, na.rm = TRUE) / den
    fdr[m]   <- 100 * sum(s == signo_sig, na.rm = TRUE) / den
  }

  cat(sprintf("\n%s | %s ~ %s (%s)   [descartado por n<%d: %.1f%% de la superficie]\n",
              parque, vi, var, SIGNO, N_MIN, descartado))
  cat("   mes:  ", paste(sprintf("%5s", month.abb), collapse = " "), "\n")
  cat("   crudo:", paste(sprintf("%5.1f", crudo), collapse = " "),
      sprintf("  media %.1f%%\n", mean(crudo, na.rm = TRUE)))
  cat("   FDR:  ", paste(sprintf("%5.1f", fdr), collapse = " "),
      sprintf("  media %.1f%%\n", mean(fdr, na.rm = TRUE)))

  resumen <- rbind(resumen, data.frame(parque = parque, indice = vi, variable = var,
                                       crudo = mean(crudo, na.rm = TRUE),
                                       fdr = mean(fdr, na.rm = TRUE),
                                       descartado = descartado))

  # --- escribir con los mismos nombres de variable que esperan las figuras ---
  # Se arrastra n_obs para que estos ficheros sean autocontenidos: así se puede
  # reprocesar con otro N_MIN sin volver a los máximos.
  f_out <- file.path(DIR_SALIDA, basename(f_in))
  dim_lon <- ncdim_def("lon", "degrees_east", lon)
  dim_lat <- ncdim_def("lat", "degrees_north", lat)
  dim_month <- ncdim_def("month", "month", 1:12)

  v_cor <- ncvar_def("max_correlation", "double", list(dim_lon, dim_lat, dim_month), -9999)
  v_sig <- ncvar_def("significance",    "double", list(dim_lon, dim_lat, dim_month), -9999)
  v_p   <- ncvar_def("p_value",         "double", list(dim_lon, dim_lat, dim_month), -9999)
  v_esc <- ncvar_def("scale",           "double", list(dim_lon, dim_lat, dim_month), -9999)
  v_n   <- ncvar_def("n_obs",           "double", list(dim_lon, dim_lat, dim_month), -9999)

  out <- nc_create(f_out, list(v_cor, v_sig, v_p, v_esc, v_n))
  ncvar_put(out, v_cor, cor_max)
  ncvar_put(out, v_sig, sig_new)
  ncvar_put(out, v_p,   p_val)
  ncvar_put(out, v_esc, escala)
  ncvar_put(out, v_n,   n_obs)
  ncatt_put(out, 0, "correcciones",
            sprintf("n>=%d valid years; BH FDR q=%.2f per park-month field; no scale correction", N_MIN, Q_FDR))
  ncatt_put(out, 0, "p_value", "uncorrected two-sided p-value at the selected scale; FDR applied to it")
  ncatt_put(out, 0, "referencia_fdr", "Wilks 2016 BAMS 97:2263-2273")
  nc_close(out)
  cat("   guardado:", basename(f_out), "\n")
}

cat("\n\n================ RESUMEN ", toupper(SIGNO), " ================\n", sep = "")
if (nrow(resumen)) {
  agg <- aggregate(cbind(crudo, fdr) ~ parque, data = resumen, FUN = mean)
  agg$retiene <- round(100 * agg$fdr / agg$crudo)
  agg[, 2:3] <- round(agg[, 2:3], 1)
  print(agg[order(-agg$fdr), ], row.names = FALSE)
  cat(sprintf("\nGLOBAL: %.1f%% -> %.1f%% (retiene %.0f%%)\n",
              mean(resumen$crudo), mean(resumen$fdr),
              100 * mean(resumen$fdr) / mean(resumen$crudo)))
}
