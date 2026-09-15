# ==============================================================================
# CORRELACIONES POSITIVAS MÁXIMAS POR TIPO DE VEGETACIÓN Y ESCALA (SPEI)
# ==============================================================================
# Equivalente a generar_csv_pos_clima.R para el SPEI. Para cada parque, mes y
# tipo de vegetación (Id_SNVeg), cuenta los píxeles válidos del tipo y, de
# ellos, los que tienen correlación positiva significativa, desglosados por la
# escala del SPEI en que se alcanza la correlación máxima. Cada píxel cuenta
# una sola vez, en la escala de su máximo.
#
# Lee los máximos con umbral de años válidos y FDR (03_fdr/). Los píxeles
# descartados por el umbral tienen significance NA y no cuentan ni en n ni en
# n_total. Se emite la rejilla completa tipo x mes x escala, con n = 0 donde no
# hay píxeles significativos.
#
# Entradas: correlaciones_v2/maximos_cor_positive_fdr/<parque>_max_pos_correlation_spei_<indice>.nc
#           mapas_vegetacion/sistemasNaturales_vegetacion/SNV_parques/SNV_<parque>.shp
# Salida  : resultados_v2/snv_respuesta_spei_<indice>_absoluto_positivo.csv
#           columnas: Id_SNVeg, n_total, mes, escala, n, parque
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

# 🏞️ Parques a procesar (Sierra de las Nieves no tiene mapa de vegetación)
parques <- c("aiguestortes", "cabaneros", "cabrera", "daimiel", "donana",
             "guadarrama", "cies", "monfrague", "ordesa", "picoseuropa",
             "sierranevada")

# 📁 Directorios de entrada
base_nc <- paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_positive_fdr")
base_snv <- paste0(ROOT, "/mapas_vegetacion/sistemasNaturales_vegetacion/SNV_parques")
indice <- "ndvi"    # "ndvi" o "kndvi"

# 🗓️ Meses y escalas
meses <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
n_meses <- length(meses)
escalas <- c(1, 3, 6, 9, 12)

# 📦 Lista para almacenar todos los resultados
todos_resultados <- list()

# 🔁 Iterar por parque
for (parque in parques) {
  message("🔄 Procesando parque: ", parque)

  # Rutas a NetCDF y shapefile
  nc_path <- file.path(base_nc, sprintf("%s_max_pos_correlation_spei_%s.nc", parque, indice))
  snv_path <- file.path(base_snv, sprintf("SNV_%s.shp", parque))

  # Leer raster y shapefile
  r_esc <- rast(nc_path, subds = "scale")          # Raster con escalas SPEI
  r_sig <- rast(nc_path, subds = "significance")   # Raster con significancia
  snv <- vect(snv_path)

  # Reproyectar shapefile al CRS del raster
  snv_proj <- project(snv, crs(r_esc))
  snv_proj$Id_SNVeg <- as.factor(snv_proj$Id_SNVeg)
  snv_proj$ID <- 1:nrow(snv_proj)

  # 🔁 Iterar por mes
  for (m in 1:n_meses) {
    escala_m <- r_esc[[m]]
    sig_m <- r_sig[[m]]
    names(escala_m) <- "scale"

    # Píxeles válidos: los descartados por el umbral de años tienen
    # significance NA aunque conserven scale
    escala_m_valida <- ifel(is.na(sig_m), NA, escala_m)
    names(escala_m_valida) <- "scale"

    # =============================
    # 1️⃣ Total de píxeles válidos por tipo
    # =============================
    total_pixeles <- terra::extract(escala_m_valida, snv_proj) %>%
      left_join(as.data.frame(snv_proj)[, c("ID", "Id_SNVeg")], by = "ID") %>%
      filter(!is.na(scale), !is.na(Id_SNVeg)) %>%
      count(Id_SNVeg, name = "n_total") %>%
      mutate(mes = m)

    # =============================
    # 2️⃣ Píxeles significativos (ifel es seguro con NA en significance)
    # =============================
    escala_m_filtrada <- ifel(is.na(sig_m) | sig_m != 2, NA, escala_m)
    names(escala_m_filtrada) <- "scale"

    df_extract <- terra::extract(escala_m_filtrada, snv_proj)

    df_extract <- df_extract %>%
      left_join(as.data.frame(snv_proj)[, c("ID", "Id_SNVeg")], by = "ID") %>%
      filter(!is.na(scale), !is.na(Id_SNVeg)) %>%
      mutate(
        escala = as.character(scale),
        mes = m
      )

    # =============================
    # 3️⃣ Contar por escala y completar la rejilla tipo x escala con n = 0
    # =============================
    resumen <- df_extract %>%
      count(Id_SNVeg, escala, mes, name = "n")

    resumen_final <- total_pixeles %>%
      crossing(escala = as.character(escalas)) %>%
      left_join(resumen, by = c("Id_SNVeg", "mes", "escala")) %>%
      mutate(n = ifelse(is.na(n), 0L, n), parque = parque)

    todos_resultados[[length(todos_resultados) + 1]] <- resumen_final
  }
}

# 📊 Unir todos los resultados
df_final <- bind_rows(todos_resultados)

# 💾 Guardar como CSV
dir.create(paste0(ROOT, "/articulo_2/resultados_v2"), showWarnings = FALSE, recursive = TRUE)
output_csv <- sprintf(paste0(ROOT, "/articulo_2/resultados_v2/snv_respuesta_spei_%s_absoluto_positivo.csv"), indice)
write_csv(df_final, output_csv)

cat("✅ CSV generado en:\n", output_csv, "\n")
