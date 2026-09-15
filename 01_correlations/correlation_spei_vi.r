# ================================================================================
# CORRELACIONES MENSUALES ENTRE ÍNDICES DE VEGETACIÓN Y SPEI
# ================================================================================
# Equivalente a correlation_vi_nc.r para el SPEI, que tiene su propia cadena de
# ficheros: un NetCDF por escala (spei01, spei03, spei06, spei09, spei12) en la
# misma rejilla que el resto del clima (1,1 km, EPSG:23030), del que se usan las
# 480 capas del periodo 1984-2023. El SPEI ya es un índice acumulado, de modo que
# la "escala" es la del propio SPEI.
#
# Procedimiento por píxel y mes (idéntico al de las variables climáticas):
#   enmascarado de años sin dato en alguna de las dos series; detrend por mínimos
#   cuadrados contra el año real; Pearson r y p-valor bilateral con n-3 grados de
#   libertad; Spearman rho en una submuestra reproducible de píxeles.
#
# Salida: NetCDF [lon, lat, month] con correlation, p_value, significance
# (2/-2: p<=0.05 sin corregir; 1/-1: no), n_obs, correlation_spearman y
# p_value_spearman. El umbral de años válidos y el FDR se aplican en 03_fdr/.
#
# Entradas: MOSAICOS/mosaicos_all/<parque>_<INDICE>_filled_clean.nc
#           CLIMA/indices/SPEI_monthly/spei<ss>.nc
#           limites_red/.../Limites_PN.shp
# Salida  : articulo_2/correlaciones_v2/<parque>/correlation_<indice>_spei_scale_<s>_<parque>.nc
# ================================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(ncdf4)
library(dplyr)
library(terra)

# --- Clima al grid del índice: reproyección, recorte y remuestreo bilineal ---
extender_clima_al_ndvi <- function(clima_stack, parque_vect, vi_rast) {
  crs(clima_stack) <- "EPSG:23030"
  clima_subset <- clima_stack[[277:756]]                  # enero 1984 - diciembre 2023
  clima_proj <- project(clima_subset, crs(vi_rast))
  clima_crop <- mask(clima_proj, parque_vect, touches = TRUE)
  clima_resampled <- resample(clima_crop, vi_rast, method = "bilinear")
  return(clima_resampled)
}

# --- Igual, con relleno focal previo: parques de borde (islas, costa, límite
#     con Francia) que la rejilla climática no cubre por completo ---
extender_clima_al_ndvi_cies <- function(clima_stack, parque_vect, vi_rast, n_iter = 1, w = 3) {
  crs(clima_stack) <- "EPSG:23030"
  clima_subset <- clima_stack[[277:756]]
  clima_subset <- focal(clima_subset, w = w, fun = mean, na.policy = "only", na.rm = TRUE)
  clima_proj <- project(clima_subset, crs(vi_rast))
  clima_crop <- mask(clima_proj, parque_vect, touches = TRUE)
  clima_filled <- clima_crop

  for (i in 1:n_iter) {
    clima_filled <- focal(clima_filled, w = w, fun = mean, na.policy = "only", na.rm = TRUE)
  }

  clima_resampled <- resample(clima_filled, vi_rast, method = "bilinear")
  return(clima_resampled)
}


# --- Residuos de la recta de mínimos cuadrados z ~ tt (equivale a
#     residuals(lm(z ~ tt)), sin el coste de construir el modelo) ---
quita_tendencia <- function(z, tt) {
  tc <- tt - mean(tt)
  z - mean(z) - tc * (sum(tc * z) / sum(tc * tc))
}


# --- Configuración ---
index_list <- c('NDVI','KNDVI')
path_spei <- paste0(ROOT, '/CLIMA/indices/SPEI_monthly/')
vi_path <- paste0(ROOT, '/MOSAICOS/mosaicos_all/')
out_path_base <- paste0(ROOT, '/articulo_2/correlaciones_v2/')
shp_path <- paste0(ROOT, '/limites_red/desc_Red_PN_LIM_Enero2023/Limites_PN.shp')

parques <- c("aiguestortes","cabaneros","cabrera","cies","daimiel","donana",
             "guadarrama","monfrague","ordesa","picoseuropa",
             "sierranevada","sierranieves")

# Parques que la rejilla climática peninsular no cubre por completo (superficie
# sin dato nativo: Islas Atlánticas 53,9 %, Cabrera 14,1 %, Ordesa 2,2 %,
# Doñana 1,2 %); en ellos se aplica el relleno focal.
parques_gapfill <- c("cabrera", "cies", "ordesa", "donana")

anios <- 1984:2023                        # 480 meses = 40 años

# Spearman en una submuestra de píxeles (la misma para todas las escalas y
# meses del parque). FRAC_SPEARMAN <- 1 lo calcula en todos.
calcular_spearman <- TRUE
FRAC_SPEARMAN <- 0.05
set.seed(42)                              # submuestra reproducible

# --- Leer shapefile general de parques ---
parques_vect <- vect(shp_path)
parques_vect <- project(parques_vect, "EPSG:4326")

for (parque in parques) {
  cat("\n==============================\nProcesando parque:", parque, "\n==============================\n")
  out_path <- file.path(out_path_base, parque)
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  funcion_ext <- if (parque %in% parques_gapfill) extender_clima_al_ndvi_cies else extender_clima_al_ndvi
  cat("   relleno focal del clima:", parque %in% parques_gapfill, "\n")

  for (vi in index_list) {
    vi_file <- paste0(vi_path, parque, '_', vi, '_filled_clean.nc')
    vi_nc <- nc_open(vi_file)
    lon <- ncvar_get(vi_nc, 'lon')
    lat <- ncvar_get(vi_nc, 'lat')
    nc_close(vi_nc)

    vi_rast_full <- rast(vi_file)
    vi_rast <- vi_rast_full[[grep(paste0("^", vi), names(vi_rast_full))]]

    # --- Aplicar máscara del parque con centroides ---
    parque_vect <- parques_vect[parques_vect$name_pn == parque, ]
    parque_vect <- project(parque_vect, crs(vi_rast))
    vi_rast <- mask(vi_rast, parque_vect, touches = FALSE)

    vi_data <- as.array(vi_rast)
    time_dim <- dim(vi_data)[3]
    lon_dim <- dim(vi_data)[2]
    lat_dim <- dim(vi_data)[1]

    # Píxeles de la submuestra de Spearman
    submuestra_sp <- matrix(runif(lon_dim * lat_dim) < FRAC_SPEARMAN, lon_dim, lat_dim)

    for (scale in c(1, 3, 6, 9, 12)) {
      clim_file <- file.path(path_spei, paste0("spei", sprintf("%02d", scale), ".nc"))

      if (!file.exists(clim_file)) {
        cat("⚠️  Archivo no encontrado para SPEI escala", scale, ":", clim_file, "\n")
        next
      }

      cat("\nProcesando SPEI escala:", scale, "mes(es) para índice", vi, "\n")

      clima_stack <- rast(clim_file)
      clima_resampled <- funcion_ext(clima_stack, parque_vect, vi_rast)
      var_data_subset <- as.array(clima_resampled)

      correlation <- array(NA, dim = c(lon_dim, lat_dim, 12))
      significant <- array(NA, dim = c(lon_dim, lat_dim, 12))
      p_value_array <- array(NA, dim = c(lon_dim, lat_dim, 12))
      n_obs <- array(NA, dim = c(lon_dim, lat_dim, 12))
      correlation_s <- array(NA, dim = c(lon_dim, lat_dim, 12))
      p_value_s <- array(NA, dim = c(lon_dim, lat_dim, 12))

      for (m in 1:12) {
        cat("    Procesando mes:", m, "\n")
        vi_month <- vi_data[, , seq(m, time_dim, by = 12)]
        var_month <- var_data_subset[, , seq(m, time_dim, by = 12)]

        for (j in 1:lon_dim) {
          for (i in 1:lat_dim) {
            s_series <- vi_month[i, j, ]
            v_series <- var_month[i, j, ]

            if (sum(!is.na(s_series)) > 1 && sum(!is.na(v_series)) > 1) {

              # Años con dato en las dos series; la tendencia se ajusta
              # después del enmascarado y contra el año real
              mask <- !is.na(s_series) & !is.na(v_series)
              n_valid <- sum(mask)
              n_obs[j, i, m] <- n_valid

              if (n_valid >= 5) {
                x  <- s_series[mask]
                y  <- v_series[mask]
                tt <- anios[mask]

                if (sd(x) > 0 && sd(y) > 0) {
                  x <- quita_tendencia(x, tt)
                  y <- quita_tendencia(y, tt)

                  r <- unname(cor(x, y))
                  correlation[j, i, m] <- r

                  # df = n-3: se ha estimado y eliminado una tendencia
                  p_val <- 2 * pt(-abs(r * sqrt((n_valid - 3) / (1 - r^2))), n_valid - 3)
                  p_value_array[j, i, m] <- p_val

                  significant[j, i, m] <- ifelse(
                    r >= 0,
                    ifelse(p_val <= 0.05, 2, 1),
                    ifelse(p_val <= 0.05, -2, -1)
                  )

                  # Spearman sobre los mismos residuos, solo en la submuestra
                  if (calcular_spearman && submuestra_sp[j, i]) {
                    rs <- unname(cor(x, y, method = "spearman"))
                    correlation_s[j, i, m] <- rs
                    p_value_s[j, i, m] <- 2 * pt(-abs(rs * sqrt((n_valid - 3) / (1 - rs^2))),
                                                 n_valid - 3)
                  }
                }
              }
            }
          }
        }
      }

      outfile <- file.path(out_path, paste0('correlation_', tolower(vi), '_spei_scale_', scale, '_', parque, '.nc'))
      lon_out <- xFromCol(vi_rast, 1:lon_dim)
      lat_out <- yFromRow(vi_rast, 1:lat_dim)

      dim_lon <- ncdim_def("lon", "degrees_east", lon_out)
      dim_lat <- ncdim_def("lat", "degrees_north", lat_out)
      dim_month <- ncdim_def("month", "1=Jan,...,12=Dec", 1:12)

      var_cor <- ncvar_def("correlation", "double", list(dim_lon, dim_lat, dim_month), -9999)
      var_sig <- ncvar_def("significance", "double", list(dim_lon, dim_lat, dim_month), -9999)
      var_p <- ncvar_def("p_value", "double", list(dim_lon, dim_lat, dim_month), -9999)
      var_n <- ncvar_def("n_obs", "double", list(dim_lon, dim_lat, dim_month), -9999)
      var_cor_s <- ncvar_def("correlation_spearman", "double", list(dim_lon, dim_lat, dim_month), -9999)
      var_p_s <- ncvar_def("p_value_spearman", "double", list(dim_lon, dim_lat, dim_month), -9999)

      out_nc <- nc_create(outfile, list(var_cor, var_sig, var_p, var_n, var_cor_s, var_p_s))
      ncvar_put(out_nc, var_cor, correlation)
      ncvar_put(out_nc, var_sig, significant)
      ncvar_put(out_nc, var_p, p_value_array)
      ncvar_put(out_nc, var_n, n_obs)
      ncvar_put(out_nc, var_cor_s, correlation_s)
      ncvar_put(out_nc, var_p_s, p_value_s)
      ncatt_put(out_nc, 0, "version", "v2: OLS detrend after masking, df=n-3")
      ncatt_put(out_nc, 0, "significance", "code 2/-2: p<=0.05 uncorrected; FDR applied downstream")
      nc_close(out_nc)

      cat("  ✅ Guardado:", outfile, "\n")
      rm(var_data_subset); gc()
    }
    cat("✔️ Correlaciones SPEI completadas para índice", vi, "\n")
  }
}
