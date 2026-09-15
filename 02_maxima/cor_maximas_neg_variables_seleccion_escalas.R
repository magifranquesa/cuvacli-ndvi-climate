# ==============================================================================
# CORRELACIÓN NEGATIVA MÁXIMA ENTRE ESCALAS DE ACUMULACIÓN (Tmax, Tmin, Pr)
# ==============================================================================
# Para cada parque, índice (NDVI, kNDVI), variable, píxel y mes, selecciona
# entre las cinco escalas de acumulación (1, 3, 6, 9, 12 meses) la correlación
# más baja (la negativa más fuerte), y guarda con ella la escala en que ocurre,
# su p-valor sin corregir, su código de significancia y el número de años
# válidos.
#
# Entrada: correlaciones_v2/<parque>/correlation_<indice>_<var>_scale_<s>_<parque>.nc
# Salida : correlaciones_v2/maximos_cor_negative/<parque>_max_neg_correlation_<var>_<indice>_selscales.nc
#          variables: max_correlation, significance, p_value, scale, n_obs
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(ncdf4)
library(abind)

# =============================
# ⚙️ Configuración general
# =============================

variables <- c("tmax", "tmin", "pr")

indices <- c("NDVI", "KNDVI")

parques <- c("aiguestortes","cabaneros","cabrera","cies","daimiel","donana",
             "guadarrama","monfrague","ordesa","picoseuropa",
             "sierranevada","sierranieves")

scales <- c("1", "3", "6", "9", "12")

base_input_dir <- paste0(ROOT, "/articulo_2/correlaciones_v2")
output_dir <- paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_negative")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================
# 🔁 Bucle principal
# =============================

for (parque in parques) {
  input_dir <- file.path(base_input_dir, parque)

  for (index in indices) {
    for (var in variables) {
      cat("\n📍 Procesando:", parque, "-", index, "-", var, "\n")

      correlation_list <- significance_list <- p_value_list <- n_obs_list <- list()
      dims <- NULL

      for (scale_val in scales) {
        file_name <- sprintf("correlation_%s_%s_scale_%s_%s.nc", tolower(index), var, scale_val, parque)
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

        correlation_list[[scale_val]] <- ncvar_get(nc, "correlation")
        significance_list[[scale_val]] <- ncvar_get(nc, "significance")
        p_value_list[[scale_val]] <- ncvar_get(nc, "p_value")
        n_obs_list[[scale_val]] <- ncvar_get(nc, "n_obs")

        nc_close(nc)
      }

      if (length(Filter(Negate(is.null), correlation_list)) == 0) {
        cat("⚠️ No se encontraron archivos válidos para", parque, var, index, "\n")
        next
      }

      # Apilar arrays
      correlation_array <- abind::abind(correlation_list[scales], along = 0)
      significance_array <- abind::abind(significance_list[scales], along = 0)
      p_value_array <- abind::abind(p_value_list[scales], along = 0)
      n_obs_array <- abind::abind(n_obs_list[scales], along = 0)

      # Inicializar arrays de salida
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
            idx <- which.min(vals)
            max_correlation[i, j, m] <- correlation_array[idx, i, j, m]
            significance[i, j, m] <- significance_array[idx, i, j, m]
            p_value[i, j, m] <- p_value_array[idx, i, j, m]
            scale_array[i, j, m] <- as.integer(scales[idx])
            n_obs_sel[i, j, m] <- n_obs_array[idx, i, j, m]
          }
        }
      }

      # Crear archivo NetCDF
      outfile <- file.path(output_dir, sprintf("%s_max_neg_correlation_%s_%s_selscales.nc", parque, var, tolower(index)))

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
}
