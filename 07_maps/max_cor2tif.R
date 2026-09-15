# ============================================
# CAPAS GEOTIFF PARA LOS MAPAS (FIGURA 8 Y S1-S8)
# ============================================
# Para cada parque, variable (SPEI, Pr, Tmax, Tmin) y signo, exporta tres
# GeoTIFF a partir de los máximos con umbral de años y FDR (03_fdr/):
#   <parque>_<var>_<pos|neg>_<indice>_cor_max.tif   correlación más fuerte del año
#   <parque>_<var>_<pos|neg>_<indice>_mes_max.tif   mes en que ocurre (1-12)
#   <parque>_<var>_<pos|neg>_<indice>_esc_max.tif   escala de acumulación en ese mes
# Solo se consideran los píxeles-mes significativos del signo correspondiente;
# para las correlaciones negativas "más fuerte" es la más negativa (min).
#
# Los mapas se componen después en ArcGIS Pro a partir de estas capas.
#
# Entrada: correlaciones_v2/maximos_cor_{positive,negative}_fdr/
# Salida : figuras/mapas_v2/fdr/rasters_max_cor_{positive,negative}/<var>/
# ============================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

# 📚 Librerías necesarias
library(ncdf4)
library(terra)

# =============================
# 🧭 Configuración
# =============================
parques <- c("aiguestortes", "cabaneros", "daimiel", "donana", "guadarrama",
             "monfrague", "ordesa", "picoseuropa", "sierranevada",
             "sierranieves", "cabrera", "cies")

indice    <- "ndvi"    # "ndvi" o "kndvi"
variables <- c("spei", "pr", "tmax", "tmin")
tipos     <- c("max_pos", "max_neg")

for (tipo in tipos) {

  base_dir <- sprintf(paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_%s_fdr"),
                      ifelse(tipo == "max_pos", "positive", "negative"))

  output_root_dir <- sprintf(paste0(ROOT, "/articulo_2/figuras/mapas_v2/fdr/rasters_max_cor_%s"),
                             ifelse(tipo == "max_pos", "positive", "negative"))

  # Código de significancia y agregación según el signo: en las negativas la
  # correlación más fuerte es la mínima
  cod_sig <- if (tipo == "max_pos") 2 else -2
  agrega  <- if (tipo == "max_pos") max else min
  cual    <- if (tipo == "max_pos") which.max else which.min

  for (variable in variables) {

    output_raster_dir <- file.path(output_root_dir, variable)
    dir.create(output_raster_dir, showWarnings = FALSE, recursive = TRUE)
    message("\n===== ", tipo, " | ", variable, " -> ", output_raster_dir)

    # =============================
    # 🔁 Procesamiento por parque
    # =============================
    for (parque in parques) {
      message("Procesando: ", parque)

      # Nombre del archivo de entrada (el SPEI no lleva sufijo _selscales)
      filename <- if (variable == "spei") {
        paste0(parque, "_", tipo, "_correlation_spei_", indice, ".nc")
      } else {
        paste0(parque, "_", tipo, "_correlation_", variable, "_", indice, "_selscales.nc")
      }

      file <- file.path(base_dir, filename)

      if (!file.exists(file)) {
        warning("No se encuentra el archivo: ", file)
        next
      }

      # Leer NetCDF
      nc  <- nc_open(file)
      cor <- ncvar_get(nc, "max_correlation")
      sig <- ncvar_get(nc, "significance")
      esc <- ncvar_get(nc, "scale")
      lon <- ncvar_get(nc, "lon")
      lat <- ncvar_get(nc, "lat")
      nc_close(nc)

      # Máscara de significancia a prueba de NA (los píxeles descartados por el
      # umbral de años tienen significance NA)
      cor[is.na(sig) | sig != cod_sig] <- NA

      # Correlación extrema del año y su mes
      dims    <- dim(cor)
      cor_max <- apply(cor, c(1, 2), function(x) if (all(is.na(x))) NA else agrega(x, na.rm = TRUE))
      mes_max <- apply(cor, c(1, 2), function(x) if (all(is.na(x))) NA else cual(x))

      # Escala correspondiente al valor extremo
      escala_max <- matrix(NA, nrow = dims[1], ncol = dims[2])
      for (i in 1:dims[1]) {
        for (j in 1:dims[2]) {
          m <- mes_max[i, j]
          if (!is.na(m)) {
            escala_max[i, j] <- esc[i, j, m]
          }
        }
      }

      # Crear rasters
      ext_r <- ext(min(lon), max(lon), min(lat), max(lat))
      r_cor <- rast(t(cor_max),    extent = ext_r, crs = "EPSG:4326")
      r_mes <- rast(t(mes_max),    extent = ext_r, crs = "EPSG:4326")
      r_esc <- rast(t(escala_max), extent = ext_r, crs = "EPSG:4326")

      # Exportar
      sfx <- sprintf("%s_%s_%s", variable, ifelse(tipo == "max_pos", "pos", "neg"), indice)
      writeRaster(r_cor, file.path(output_raster_dir, paste0(parque, "_", sfx, "_cor_max.tif")), overwrite = TRUE)
      writeRaster(r_mes, file.path(output_raster_dir, paste0(parque, "_", sfx, "_mes_max.tif")), overwrite = TRUE)
      writeRaster(r_esc, file.path(output_raster_dir, paste0(parque, "_", sfx, "_esc_max.tif")), overwrite = TRUE)
    }
  }
}

message("\nListo: capas en figuras/mapas_v2/fdr/")
