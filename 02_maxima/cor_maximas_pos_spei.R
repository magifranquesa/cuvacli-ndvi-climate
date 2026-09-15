# ==============================================================================
# CORRELACIÓN POSITIVA MÁXIMA ENTRE ESCALAS DEL SPEI
# ==============================================================================
# Para cada parque, índice (NDVI, kNDVI), píxel y mes, selecciona entre las
# cinco escalas del SPEI (1, 3, 6, 9, 12 meses) la correlación más alta, y
# guarda con ella la escala, su p-valor sin corregir, su código de
# significancia y el número de años válidos.
#
# El SPEI tiene su propia cadena de ficheros: la salida no lleva el sufijo
# "_selscales" y los scripts de figuras de SPEI leen ese nombre.
#
# Entrada: correlaciones_v2/<parque>/correlation_<indice>_spei_scale_<s>_<parque>.nc
# Salida : correlaciones_v2/maximos_cor_positive/<parque>_max_pos_correlation_spei_<indice>.nc
#          variables: max_correlation, significance, p_value, scale, n_obs
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(ncdf4)
library(abind)

# =============================
# ⚙️ Configuración general
# =============================

indices <- c("NDVI", "KNDVI")
scales <- c(1, 3, 6, 9, 12)

parques <- c("aiguestortes","cabaneros","cabrera","cies","daimiel","donana",
             "guadarrama","monfrague","ordesa","picoseuropa",
             "sierranevada","sierranieves")

base_input_dir <- paste0(ROOT, "/articulo_2/correlaciones_v2")
output_dir <- paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_positive")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================
# 🔁 Bucle por índice y parque
# =============================

for (index in indices) {
  for (parque in parques) {
    cat("\n📍 Procesando parque:", parque, "con índice", index, "\n")

    input_dir <- file.path(base_input_dir, parque)

    correlation_list <- significance_list <- p_value_list <- n_obs_list <- list()
    dims <- NULL

    for (s in seq_along(scales)) {
      scale_val <- scales[s]
      file_name <- sprintf("correlation_%s_spei_scale_%d_%s.nc", tolower(index), scale_val, parque)
      file_path <- file.path(input_dir, file_name)

      if (!file.exists(file_path)) {
        warning("❌ Archivo no encontrado: ", file_path)
        next
      }

      nc <- nc_open(file_path)
      if (is.null(dims)) {
        lon <- ncvar_get(nc, "lon")
        lat <- ncvar_get(nc, "lat")
        dims <- c(length(lon), length(lat), 12)
      }

      correlation_list[[s]] <- ncvar_get(nc, "correlation")
      significance_list[[s]] <- ncvar_get(nc, "significance")
      p_value_list[[s]] <- ncvar_get(nc, "p_value")
      n_obs_list[[s]] <- ncvar_get(nc, "n_obs")
      nc_close(nc)
    }

    if (length(Filter(Negate(is.null), correlation_list)) == 0) {
      cat("⚠️ No se encontraron archivos válidos para", parque, "y el índice", index, "\n")
      next
    }

    # Apilar en arrays 4D: [scale, lon, lat, month]
    correlation_array <- abind::abind(correlation_list, along = 0)
    significance_array <- abind::abind(significance_list, along = 0)
    p_value_array <- abind::abind(p_value_list, along = 0)
    n_obs_array <- abind::abind(n_obs_list, along = 0)

    # Arrays de salida
    max_correlation <- array(NA, dims)
    significance <- array(NA, dims)
    p_value <- array(NA, dims)
    scale_array <- array(NA, dims)
    n_obs_sel <- array(NA, dims)

    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        for (m in 1:12) {
          vals <- correlation_array[, i, j, m]
          if (all(is.na(vals))) next
          idx <- which.max(vals)
          max_correlation[i, j, m] <- correlation_array[idx, i, j, m]
          significance[i, j, m] <- significance_array[idx, i, j, m]
          p_value[i, j, m] <- p_value_array[idx, i, j, m]
          scale_array[i, j, m] <- scales[idx]
          n_obs_sel[i, j, m] <- n_obs_array[idx, i, j, m]
        }
      }
    }

    # =============================
    # 💾 Guardar archivo NetCDF
    # =============================

    outfile <- file.path(output_dir, paste0(parque, "_max_pos_correlation_spei_", tolower(index), ".nc"))

    dim_lon <- ncdim_def("lon", "degrees_east", lon)
    dim_lat <- ncdim_def("lat", "degrees_north", lat)
    dim_time <- ncdim_def("month", "1-12", 1:12)

    vars <- list(
      ncvar_def("max_correlation", "double", list(dim_lon, dim_lat, dim_time), -9999),
      ncvar_def("significance", "double", list(dim_lon, dim_lat, dim_time), -9999),
      ncvar_def("p_value", "double", list(dim_lon, dim_lat, dim_time), -9999),
      ncvar_def("scale", "integer", list(dim_lon, dim_lat, dim_time), -9999),
      ncvar_def("n_obs", "double", list(dim_lon, dim_lat, dim_time), -9999)
    )

    nc_out <- nc_create(outfile, vars)
    ncvar_put(nc_out, vars[[1]], max_correlation)
    ncvar_put(nc_out, vars[[2]], significance)
    ncvar_put(nc_out, vars[[3]], p_value)
    ncvar_put(nc_out, vars[[4]], scale_array)
    ncvar_put(nc_out, vars[[5]], n_obs_sel)
    ncatt_put(nc_out, 0, "scales_considered", paste(scales, collapse = ","))
    nc_close(nc_out)

    cat("✅ Guardado:", outfile, "\n")
  }
}
