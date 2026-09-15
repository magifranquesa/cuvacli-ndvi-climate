# ==============================================================================
# COMBINA LOS CSV POR VARIABLE EN UN ÚNICO CSV — CORRELACIONES NEGATIVAS
# ==============================================================================
# Une los CSV que escribe snv_respuesta_clima_mensual_neg.R (pr, tmax, tmin; el
# SPEI negativo no se analiza) en un solo fichero con una columna 'variable', que es
# el que leen los scripts de figuras de esta carpeta. Se detiene si falta
# alguno de los CSV de entrada.
#
# Entrada: resultados_v2/snv_respuesta_clima_mensual_negativo_<indice>_<var>.csv
# Salida : resultados_v2/snv_respuesta_clima_mensual_negativo_<indice>_completo.csv
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(dplyr)
library(readr)
library(purrr)

# 📁 Rutas
dir_base <- paste0(ROOT, "/articulo_2/resultados_v2")

indice <- "ndvi"   # "ndvi" o "kndvi"

# 📄 Archivos y sus etiquetas
vars <- c("pr", "tmax", "tmin")   # el SPEI negativo no se analiza
archivos <- setNames(
  sprintf("snv_respuesta_clima_mensual_negativo_%s_%s.csv", indice, vars),
  vars)


faltan <- archivos[!file.exists(file.path(dir_base, archivos))]
if (length(faltan)) {
  stop("Faltan estos CSV en ", dir_base, ":\n  ", paste(faltan, collapse = "\n  "),
       "\nEjecuta antes snv_respuesta_clima_mensual_neg.R",
       " con indice <- \"", indice, "\".")
}

# 📦 Leer, etiquetar y combinar
df_combinado <- imap_dfr(archivos, function(archivo, var) {
  read_csv(file.path(dir_base, archivo), show_col_types = FALSE) %>%
    mutate(variable = var)
})

# 🧹 Asegurar tipos adecuados
df_combinado <- df_combinado %>%
  mutate(
    Id_SNVeg = as.factor(Id_SNVeg),
    parque   = as.factor(parque),
    escala   = as.integer(escala),
    mes      = as.integer(mes),
    variable = factor(variable, levels = c("pr", "tmax", "tmin")),
    n_total  = as.integer(n_total),
    n_sig    = as.integer(n_sig),
    prop_sig = as.numeric(prop_sig)
  )

# 💾 Guardar CSV final
output_csv <- file.path(dir_base, sprintf("snv_respuesta_clima_mensual_negativo_%s_completo.csv", indice))
write_csv(df_combinado, output_csv)
cat("✅ CSV combinado y limpiado guardado en:\n", output_csv, "\n")
cat(sprintf("   %d filas | %d parques | %d tipos de vegetación\n",
            nrow(df_combinado), n_distinct(df_combinado$parque), n_distinct(df_combinado$Id_SNVeg)))
