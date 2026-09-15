# ==============================================================================
# FIGURA 9: CORRELACIONES POSITIVAS SIGNIFICATIVAS POR TIPO DE VEGETACIÓN
# (red completa de parques)
# ==============================================================================
# Un panel por tipo de vegetación con, para cada mes y variable (Pr, SPEI,
# Tmax, Tmin), el porcentaje de combinaciones píxel x escala con correlación
# positiva significativa, agregado sobre todos los parques:
#     100 * sum(n_sig) / sum(n_total)   sobre parques y escalas.
# Como un píxel puede ser significativo en varias escalas, los valores son
# menores que en las figuras basadas en máximos, donde cada píxel cuenta una
# sola vez.
#
# Entrada: resultados_v2/snv_respuesta_clima_mensual_positivo_<indice>_completo.csv
#          (combinar_csvs.R)
# Salida : figuras/SNV_v2/fdr/global_pos/snv_response_clima_ponderado.{tif,pdf,svg}
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(svglite)

indice     <- "ndvi"    # "ndvi" o "kndvi"

csv_in <- sprintf(paste0(ROOT, "/articulo_2/resultados_v2/snv_respuesta_clima_mensual_positivo_%s_completo.csv"), indice)
if (!file.exists(csv_in)) stop("No existe ", csv_in, "\nEjecuta antes combinar_csvs.R con indice <- \"", indice, "\".")

out_dir <- paste0(ROOT, "/articulo_2/figuras/SNV_v2/fdr/global_pos")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("entrada:", csv_in, "\nsalida :", out_dir, "\n")

# 📁 Cargar datos
df <- read_csv(csv_in, show_col_types = FALSE)

# 🎯 Tipos de vegetación de interés
tipos_interes <- c(
  "AlpConF", "MedConF", "MedSclF", "SMedMarF", "TempDecF",  # Bosques
  "AlpScrGr", "AtlScrub", "MedScrub", "AridScr", "Grass",   # Matorrales y herbáceas
  "Crop", "Reforestation",                                   # Antrópicas
  "HaloVeg", "HydroRipVeg", "SaltMar"                        # Halófitas e hidrófitas
)

# 🔠 Diccionario nombres completos
nombres_vegetacion <- c(
  AlpScrGr = "Alpine Scrublands and Grasslands",
  AlpConF = "Alpine Coniferous Forest",
  MedConF = "Mediterranean Coniferous Forests",
  MedScrub = "Mediterranean Scrublands",
  TempDecF = "Temperate Broadleaf Deciduous Forest",
  SMedMarF = "Semi-Mediterranean Marcescent Forests",
  MedSclF = "Mediterranean Sclerophyllous Forest",
  AtlScrub = "Atlantic Scrublands",
  AridScr = "Hyperxerophilous Garrigues and Scrublands",
  Grass = "Mediterranean and Atlantic Grasslands",
  Reforestation = "Reforestation",
  Crop = "Crops",
  HaloVeg = "Halophilous Vegetation",
  HydroRipVeg = "Hydrophilous and Riparian Vegetation",
  SaltMar = "Salt Marshes"
)

# 📆 Etiquetas de mes
meses <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# 🧼 Filtrar y calcular porcentaje ponderado
df_ponderado <- df %>%
  filter(Id_SNVeg %in% tipos_interes) %>%
  mutate(Nombre_SNVeg = factor(nombres_vegetacion[Id_SNVeg], levels = nombres_vegetacion[tipos_interes])) %>%
  group_by(Nombre_SNVeg, mes, variable) %>%
  summarise(
    prop_sig = sum(n_sig, na.rm = TRUE) / sum(n_total, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(mes_label = factor(meses[mes], levels = meses))

# 🎨 Colores por variable
colores_var <- c("pr" = "#003f5c", "spei" = "#ffa600", "tmax" = "#ef5675", "tmin" = "#7a5195")

# ylim() descarta en silencio las barras que superan el límite del eje:
# se comprueba antes de dibujar
Y_MAX <- 80
pico <- max(df_ponderado$prop_sig, na.rm = TRUE)
cat(sprintf("maximo del panel: %.1f%%  (limite del eje: %d)\n", pico, Y_MAX))
if (pico > Y_MAX) {
  stop(sprintf("Hay barras de %.1f%%, por encima del limite del eje (%d): se dibujarian cortadas.\n  Sube Y_MAX y el ylim() de abajo.", pico, Y_MAX))
}

# 📊 Crear lista de gráficos individuales
plots <- df_ponderado %>%
  split(.$Nombre_SNVeg) %>%
  { 
    nombres <- names(.)
    lapply(seq_along(.), function(i) {
      subdf <- .[[i]]
      nombre <- nombres[i]
      ggplot(subdf, aes(x = mes_label, y = prop_sig, fill = variable)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_fill_manual(
          values = colores_var,
          name = "Climate Variable",
          labels = c("pr" = "Pr", "spei" = "SPEI", "tmin" = "Tmin", "tmax" = "Tmax")
        ) +
        labs(
          x = NULL,
          y = if ((i - 1) %% 3 == 0) "% Significant Pixels" else NULL,
          title = nombre
        ) +
        ylim(0, Y_MAX) +
        theme_minimal(base_size = 13) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(size = 13, face = "bold", hjust = 0),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.8, "cm"),
          legend.title.align = 0.5,
          legend.margin = margin(t = 3)
        )
    })
  }

# 🖼️ Composición final con leyenda recogida abajo
figura <- wrap_plots(plots, ncol = 3) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# 💾 Mostrar o guardar
print(figura)
ggsave(file.path(out_dir, "snv_response_clima_ponderado.tif"), figura, width = 14, height = 12, dpi = 300)
ggsave(file.path(out_dir, "snv_response_clima_ponderado.pdf"), figura, width = 14, height = 12)
ggsave(file.path(out_dir, "snv_response_clima_ponderado.svg"), figura, width = 14, height = 12)
cat("✅ Figura 9 guardada en:\n  ", out_dir, "\n")
