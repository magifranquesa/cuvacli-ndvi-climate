# ==============================================================================
# CORRELACIONES NEGATIVAS MÁXIMAS POR TIPO DE VEGETACIÓN Y ESCALA (Tmax, Tmin, Pr)
# ==============================================================================
# Para cada variable, parque, mes y tipo de vegetación (Id_SNVeg), cuenta los
# píxeles válidos del tipo y, de ellos, los que tienen correlación negativa
# significativa, desglosados por la escala de acumulación en que se alcanza la
# correlación máxima. A diferencia de la cadena de 05_vegetation_types/, aquí
# cada píxel cuenta una sola vez, en la escala de su máximo.
#
# Lee los máximos con umbral de años válidos y FDR (03_fdr/). Los píxeles
# descartados por el umbral tienen significance NA y no cuentan ni en n ni en
# n_total. Se emite la rejilla completa tipo x mes x escala, con n = 0 donde no
# hay píxeles significativos, para que el denominador de las figuras sea
# siempre el mismo.
#
# Entradas: correlaciones_v2/maximos_cor_negative_fdr/<parque>_max_neg_correlation_<var>_<indice>_selscales.nc
#           mapas_vegetacion/sistemasNaturales_vegetacion/SNV_parques/SNV_<parque>.shp
# Salida  : resultados_v2/snv_respuesta_<var>_<indice>_absoluto_negativo.csv
#           columnas: Id_SNVeg, n_total, mes, escala, n, parque
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

# 🏞️ Parques (Sierra de las Nieves no tiene mapa de vegetación)
parques <- c("aiguestortes", "cabaneros", "cabrera", "daimiel", "donana",
             "guadarrama", "cies", "monfrague", "ordesa", "picoseuropa",
             "sierranevada")

# 🌡️ Variables climáticas a procesar
variables <- c("pr", "tmin", "tmax")

# 📁 Directorios
base_nc <- paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_negative_fdr")
base_snv <- paste0(ROOT, "/mapas_vegetacion/sistemasNaturales_vegetacion/SNV_parques")
indice <- "ndvi"    # "ndvi" o "kndvi"

# 🗓️ Meses y escalas
n_meses <- 12
escalas <- c(1, 3, 6, 9, 12)

# 🔁 Iterar por variable
for (variable in variables) {
  cat("\n🌍 Procesando variable:", variable, "\n")
  todos_resultados <- list()

  # 🔁 Iterar por parque
  for (parque in parques) {
    message("🔄 Parque: ", parque)

    # Archivos
    nc_path <- file.path(base_nc, sprintf("%s_max_neg_correlation_%s_%s_selscales.nc", parque, variable, indice))
    snv_path <- file.path(base_snv, sprintf("SNV_%s.shp", parque))

    # Raster y vector
    r_lag <- rast(nc_path, subds = "scale")
    r_sig <- rast(nc_path, subds = "significance")
    snv <- vect(snv_path)

    # Proyección
    snv_proj <- project(snv, crs(r_lag))
    snv_proj$Id_SNVeg <- as.factor(snv_proj$Id_SNVeg)
    snv_proj$ID <- 1:nrow(snv_proj)

    # 🔁 Por mes
    for (m in 1:n_meses) {
      lag_m <- r_lag[[m]]
      sig_m <- r_sig[[m]]
      names(lag_m) <- "scale"

      # Píxeles válidos: los descartados por el umbral de años tienen
      # significance NA aunque conserven scale
      lag_m_valida <- ifel(is.na(sig_m), NA, lag_m)
      names(lag_m_valida) <- "scale"

      # 1️⃣ Total de píxeles válidos por tipo
      total_pixeles <- terra::extract(lag_m_valida, snv_proj) %>%
        left_join(as.data.frame(snv_proj)[, c("ID", "Id_SNVeg")], by = "ID") %>%
        filter(!is.na(scale), !is.na(Id_SNVeg)) %>%
        count(Id_SNVeg, name = "n_total") %>%
        mutate(mes = m)

      # 2️⃣ Píxeles significativos (ifel es seguro con NA en significance)
      lag_m_filtrada <- ifel(is.na(sig_m) | sig_m != -2, NA, lag_m)
      names(lag_m_filtrada) <- "scale"

      df_extract <- terra::extract(lag_m_filtrada, snv_proj) %>%
        left_join(as.data.frame(snv_proj)[, c("ID", "Id_SNVeg")], by = "ID") %>%
        filter(!is.na(scale), !is.na(Id_SNVeg)) %>%
        mutate(
          escala = as.character(scale),
          mes = m
        )

      # 3️⃣ Contar por escala y completar la rejilla tipo x escala con n = 0
      resumen <- df_extract %>%
        count(Id_SNVeg, escala, mes, name = "n")

      resumen_final <- total_pixeles %>%
        crossing(escala = as.character(escalas)) %>%
        left_join(resumen, by = c("Id_SNVeg", "mes", "escala")) %>%
        mutate(n = ifelse(is.na(n), 0L, n), parque = parque)

      todos_resultados[[length(todos_resultados) + 1]] <- resumen_final
    }
  }

  # 📊 Guardar CSV
  df_final <- bind_rows(todos_resultados)
  dir.create(paste0(ROOT, "/articulo_2/resultados_v2"), showWarnings = FALSE, recursive = TRUE)
  output_csv <- sprintf(paste0(ROOT, "/articulo_2/resultados_v2/snv_respuesta_%s_%s_absoluto_negativo.csv"), variable, indice)
  write_csv(df_final, output_csv)

  cat("✅ CSV generado para", variable, "en:\n", output_csv, "\n")
}
