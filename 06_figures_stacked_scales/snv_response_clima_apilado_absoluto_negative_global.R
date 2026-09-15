# ==============================================================================
# FIGURAS S80-S82: CORRELACIONES NEGATIVAS MÁXIMAS POR TIPO DE VEGETACIÓN,
# APILADAS POR ESCALA DE ACUMULACIÓN (red completa de parques)
# ==============================================================================
# Equivalente a snv_response_clima_apilado_absoluto_positive_global.R para las
# correlaciones negativas (Pr, Tmax, Tmin). Un panel por tipo de vegetación con,
# para cada mes, el porcentaje de píxeles del tipo con correlación negativa
# significativa, agregado sobre todos los parques y apilado por la escala en
# que se alcanza la correlación máxima.
#
# Entrada: resultados_v2/snv_respuesta_<var>_<indice>_absoluto_negativo.csv
#          (generar_csv_neg_clima.R)
# Salida : figuras/figuras_mensuales_apilados_absolutos_SNV_v2/fdr/<indice>/global_neg/
#            snv_response_<var>_ponderado_global_negativo.{tif,pdf,svg}
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(purrr)
library(svglite)

indice     <- "ndvi"    # "ndvi" o "kndvi"
out_base   <- sprintf(paste0(ROOT, "/articulo_2/figuras/figuras_mensuales_apilados_absolutos_SNV_v2/fdr/%s/global_neg"), indice)
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
message("salida: ", out_base)

# 🎨 Paletas de colores por variable
colores_escala_roja   <- c("1"="#ffffb2","3"="#fecc5c","6"="#fd8d3c","9"="#f03b20","12"="#bd0026")
colores_escala_azul   <- c("1"="#ffffcc","3"="#a1dab4","6"="#41b6c4","9"="#2c7fb8","12"="#253494")
colores_escala_morado <- c("1"="#edf8fb","3"="#b3cde3","6"="#8c96c6","9"="#8856a7","12"="#810f7c")

paletas_variables <- list(
  tmax = colores_escala_roja,
  tmin = colores_escala_morado,
  pr   = colores_escala_azul,
  spei = colores_escala_roja
)

# 🗓️ Etiquetas de meses
meses <- c("Jan","Feb","Mar","Apr","May","Jun",
           "Jul","Aug","Sep","Oct","Nov","Dec")

# 🟢 Tipos de vegetación de interés
tipos_interes <- c(
  "AlpConF","MedConF","MedSclF","SMedMarF","TempDecF",
  "AlpScrGr","AtlScrub","MedScrub","AridScr","Grass",
  "Crop","Reforestation",
  "HaloVeg","HydroRipVeg","SaltMar"
)

# 🔠 Diccionario nombres completos
nombres_vegetacion <- c(
  AlpScrGr="Alpine Scrublands and Grasslands",
  AlpConF="Alpine Coniferous Forest",
  MedConF="Mediterranean Coniferous Forests",
  MedScrub="Mediterranean Scrublands",
  TempDecF="Temperate Broadleaf Deciduous Forest",
  SMedMarF="Semi-Mediterranean Marcescent Forests",
  MedSclF="Mediterranean Sclerophyllous Forest",
  AtlScrub="Atlantic Scrublands",
  AridScr="Hyperxerophilous Garrigues and Scrublands",
  Grass="Mediterranean and Atlantic Grasslands",
  Reforestation="Reforestation",
  Crop="Crops",
  HaloVeg="Halophilous Vegetation",
  HydroRipVeg="Hydrophilous and Riparian Vegetation",
  SaltMar="Salt Marshes"
)

# 🏷️ Títulos de leyenda por variable
titulo_leyenda <- list(
  tmax = "Tmax scale (months)",
  tmin = "Tmin scale (months)",
  pr   = "Precipitation scale (months)",
  spei = "SPEI scale (months)"
)

# 🚀 Función para generar figura por variable
generar_figura <- function(variable) {
  
  if (!variable %in% names(paletas_variables)) {
    stop(paste("❌ No hay paleta definida para la variable:", variable))
  }
  
  colores <- paletas_variables[[variable]]
  titulo <- titulo_leyenda[[variable]]
  
  # 📁 CSV de entrada (negativos)
  csv_path <- sprintf(paste0(ROOT, "/articulo_2/resultados_v2/snv_respuesta_%s_%s_absoluto_negativo.csv"), variable, indice)
  df <- read_csv(csv_path, show_col_types = FALSE)
  
  # 🧼 Preparar datos
  df_plot <- df %>%
    filter(Id_SNVeg %in% tipos_interes) %>%
    mutate(
      escala = factor(escala, levels = names(colores)),
      mes = factor(mes, levels = 1:12, labels = meses),
      Nombre_SNVeg = factor(nombres_vegetacion[Id_SNVeg],
                            levels = nombres_vegetacion[tipos_interes])
    ) %>%
    group_by(Nombre_SNVeg, mes, escala) %>%
    summarise(
      n = sum(n, na.rm = TRUE),
      n_total = sum(n_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(prop = n / n_total * 100)
  
  # 📊 Gráficos por cobertura
  plots <- df_plot %>%
    split(.$Nombre_SNVeg) %>%
    {
      lapply(seq_along(.), function(i) {
        subdf <- .[[i]]
        es_primera_col <- ((i - 1) %% 3 == 0)
        
        ggplot(subdf, aes(x = mes, y = prop, fill = escala)) +
          geom_bar(stat = "identity", position = "stack") +
          scale_fill_manual(values = colores, name = titulo) +
          labs(
            x = NULL,
            y = if (es_primera_col) "% Significant Pixels" else NULL,
            title = unique(subdf$Nombre_SNVeg)
          ) +
          ylim(0, 100) +
          theme_minimal(base_size = 13) +
          theme(
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            plot.title = element_text(size = 13, face = "bold", hjust = 0),
            axis.title.y = element_text(size = if (es_primera_col) 10 else 0)
          )
      })
    }
  
  # 🖼️ Composición final
  figura <- wrap_plots(plots, ncol = 3) +
    plot_layout(guides = "collect") & 
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
  
  # 💾 Guardar figura en carpeta de negativos
  out_path_tif <- sprintf("%s/snv_response_%s_ponderado_global_negativo.tif", out_base, variable)
  out_path_pdf <- sprintf("%s/snv_response_%s_ponderado_global_negativo.pdf", out_base, variable)
  out_path_svg <- sprintf("%s/snv_response_%s_ponderado_global_negativo.svg", out_base, variable)
  
  ggsave(out_path_tif, figura, width = 14, height = 14, dpi = 300)
  ggsave(out_path_pdf, figura, width = 14, height = 14)
  ggsave(out_path_svg, figura, width = 14, height = 14)
}

# ▶️ Ejecutar para las variables deseadas
walk(c("pr", "tmax", "tmin"), generar_figura)

