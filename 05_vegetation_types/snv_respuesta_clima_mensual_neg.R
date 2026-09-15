# ==============================================================================
# CORRELACIONES NEGATIVAS SIGNIFICATIVAS NDVI-CLIMA POR TIPO DE VEGETACIÓN
# (Tmax, Tmin, Pr)
# ==============================================================================
# Para cada parque, variable, escala de acumulación (1, 3, 6, 9, 12) y mes,
# cuenta los píxeles de cada tipo de vegetación (Id_SNVeg, Sistemas Naturales
# de Vegetación) y cuántos de ellos tienen correlación negativa significativa.
#
# Esta cadena no usa los máximos entre escalas: lee los ficheros de
# correlación por escala (01_correlations/), de modo que un píxel puede contar
# en más de una escala. La significancia se recalcula aquí con el mismo
# criterio que en 03_fdr/: se descartan los píxeles con menos de N_MIN años
# válidos y se aplica el FDR de Benjamini-Hochberg (q = Q_FDR) a todos los
# píxeles válidos del campo parque-mes-escala. La familia del FDR es el campo
# completo del parque, no el tipo de vegetación, para que el resultado de un
# tipo no dependa de cómo se estratifique el resto.
#
# Los píxeles descartados por el umbral de años no cuentan ni en n_sig ni en
# n_total.
#
# Sierra de las Nieves no está en la lista de parques: no existe mapa de
# vegetación (SNV_sierranieves.shp) para ese parque.
#
# Entradas: mapas_vegetacion/sistemasNaturales_vegetacion/SNV_parques/SNV_<parque>.shp (ETRS89 / UTM 30N)
#           correlaciones_v2/<parque>/correlation_<indice>_<var>_scale_<s>_<parque>.nc (WGS84)
# Salida  : resultados_v2/snv_respuesta_clima_mensual_negativo_<indice>_<var>.csv
#           columnas: Id_SNVeg, n_total, n_sig, prop_sig, parque, escala, mes
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(terra)
library(sf)
library(dplyr)

# 📁 Rutas
shp_dir <- paste0(ROOT, "/mapas_vegetacion/sistemasNaturales_vegetacion/SNV_parques")
nc_base_dir <- paste0(ROOT, "/articulo_2/correlaciones_v2")
out_dir <- paste0(ROOT, "/articulo_2/resultados_v2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ⚙️ Parámetros
indice    <- "ndvi"     # "ndvi" o "kndvi"
variables <- c("tmax", "tmin", "pr")
escalas <- c(1, 3, 6, 9, 12)

N_MIN <- 20      # años válidos mínimos por píxel (mismo valor que en 03_fdr/)
Q_FDR <- 0.05    # nivel del FDR (Benjamini-Hochberg)

parques <- c("aiguestortes", "cabaneros", "cabrera", "daimiel", "donana",
             "guadarrama", "cies", "monfrague", "ordesa", "picoseuropa",
             "sierranevada")

cat(sprintf("=== SNV negativas | indice=%s | n>=%d | FDR q=%.2f ===\n", indice, N_MIN, Q_FDR))

# 🔁 Bucle por variable
for (var in variables) {
  cat("Procesando variable:", var, "\n")
  resultados <- data.frame()

  for (parque in parques) {
    cat("  ➤ Parque:", parque, "\n")

    # Leer shapefile
    shp_path <- file.path(shp_dir, paste0("SNV_", parque, ".shp"))
    if (!file.exists(shp_path)) {
      cat("  ⚠️ No existe shapefile para", parque, "\n")
      next
    }

    snv <- st_read(shp_path, quiet = TRUE)
    snv <- st_transform(snv, "EPSG:4326")
    snv$ID <- 1:nrow(snv)
    snv_vect <- vect(snv)

    # Bucle por escalas
    for (escala in escalas) {
      nc_path <- file.path(nc_base_dir, parque,
                           paste0("correlation_", indice, "_", var, "_scale_", escala, "_", parque, ".nc"))

      if (!file.exists(nc_path)) {
        cat("    ⚠️ No existe:", nc_path, "\n")
        next
      }

      r_cor <- rast(nc_path, subds = "correlation")
      r_p <- rast(nc_path, subds = "p_value")
      r_n <- rast(nc_path, subds = "n_obs")

      for (mes in 1:12) {
        cor_m <- r_cor[[mes]]

        # Significancia con umbral de años y FDR sobre el campo parque-mes-escala
        cc <- values(cor_m)[, 1]
        pp <- values(r_p[[mes]])[, 1]
        nn <- values(r_n[[mes]])[, 1]
        ok <- !is.na(cc) & !is.na(pp) & !is.na(nn) & nn >= N_MIN
        s  <- rep(NA_real_, length(cc))
        if (any(ok)) {
          q <- rep(NA_real_, length(pp))
          q[ok] <- p.adjust(pp[ok], method = "BH")
          s[ok & cc >= 0] <-  1
          s[ok & cc <  0] <- -1
          s[ok & cc >= 0 & !is.na(q) & q <= Q_FDR] <-  2
          s[ok & cc <  0 & !is.na(q) & q <= Q_FDR] <- -2
        }
        sig_m <- setValues(cor_m, s)
        rm(cc, pp, nn, ok, s)

        cor_vals <- terra::extract(cor_m, snv_vect, cells = TRUE)
        sig_vals <- terra::extract(sig_m, snv_vect, cells = TRUE)

        if (nrow(cor_vals) == 0) next

        vals <- cor_vals %>%
          rename(cor = 2) %>%
          left_join(sig_vals %>% rename(sig = 2), by = c("ID", "cell")) %>%
          left_join(st_drop_geometry(snv)[, c("ID", "Id_SNVeg")], by = "ID") %>%
          # Se filtra por sig (no por cor): los píxeles descartados por el
          # umbral de años tienen cor válida pero sig NA y no deben contar
          filter(!is.na(sig))

        resumen <- vals %>%
          group_by(Id_SNVeg) %>%
          summarise(
            n_total = n(),
            n_sig = sum(sig == -2, na.rm = TRUE),
            prop_sig = round(100 * n_sig / n_total, 2),
            .groups = "drop"
          ) %>%
          mutate(parque = parque, escala = escala, mes = mes)

        resultados <- bind_rows(resultados, resumen)
      }
      rm(r_cor, r_p, r_n); gc(verbose = FALSE)
    }
  }

  # Guardar CSV por variable
  out_csv <- file.path(out_dir, paste0("snv_respuesta_clima_mensual_negativo_", indice, "_", var, ".csv"))
  write.csv(resultados, out_csv, row.names = FALSE)
  cat("✅ Guardado:", out_csv, "\n")
}
