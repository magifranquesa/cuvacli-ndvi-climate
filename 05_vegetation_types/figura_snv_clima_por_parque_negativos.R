# ==============================================================================
# FIGURAS S69-S79 (Y PANELES DE LA FIGURA 10): CORRELACIONES NEGATIVAS
# SIGNIFICATIVAS POR TIPO DE VEGETACIÓN, UNA FIGURA POR PARQUE
# ==============================================================================
# Para cada parque, un panel por tipo de vegetación presente (coberturas_por_
# parque) con, para cada mes y variable (Pr, Tmax, Tmin), el porcentaje de
# píxeles del tipo con correlación negativa significativa. Se dibuja el valor
# de cada escala de acumulación superpuesto (position_dodge sobre la variable),
# de modo que la barra visible corresponde a la escala con mayor fracción
# significativa en ese mes.
#
# Entrada: resultados_v2/snv_respuesta_clima_mensual_negativo_<indice>_completo.csv
#          (combinar_csvs_neg.R)
# Salida : figuras/SNV_v2/fdr/por_parque_neg/snv_response_negativo_<parque>.{tif,pdf,svg}
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(svglite)

# =============================
# ⚙️ Parámetros
# =============================
indice <- "ndvi"    # "ndvi" o "kndvi"

csv_in <- sprintf(paste0(ROOT, "/articulo_2/resultados_v2/snv_respuesta_clima_mensual_negativo_%s_completo.csv"), indice)
if (!file.exists(csv_in)) stop("No existe ", csv_in, "\nEjecuta antes combinar_csvs_neg.R con indice <- \"", indice, "\".")

out_dir <- paste0(ROOT, "/articulo_2/figuras/SNV_v2/fdr/por_parque_neg")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("entrada:", csv_in, "\nsalida :", out_dir, "\n")
variables <- c("tmax", "tmin", "pr")

# Tabla de parques
parques_info <- data.frame(
  parque = c("aiguestortes","cabaneros","cabrera","daimiel","donana",
             "guadarrama","cies","monfrague","ordesa","picoseuropa","sierranevada"),
  nombre_parque = c("Aigüestortes","Cabañeros","Cabrera","Daimiel","Doñana",
                    "Guadarrama","Islas Atlánticas","Monfragüe",
                    "Ordesa y Monte Perdido","Picos de Europa","Sierra Nevada"),
  ncol = c(3,3,3,3,3,3,3,3,3,3,3),
  nrow = c(1,2,1,2,2,2,1,2,2,1,2),
  stringsAsFactors = FALSE
)

# 📁 Cargar datos NEGATIVOS
df <- read_csv(csv_in, show_col_types = FALSE)

# 🧾 Diccionario nombres completos
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
  Shrub = "Shrubs",
  Dehesa = "Dehesa",
  Grass = "Mediterranean and Atlantic Grasslands",
  Reforestation = "Reforestation",
  Crop = "Crops",
  HaloVeg = "Halophilous Vegetation",
  HydroRipVeg = "Hydrophilous and Riparian Vegetation",
  SaltMar = "Salt Marshes",
  DuneSand = "Sand Dunes",
  BareArea = "Bare area",
  WaterSurf = "WaterSurf"
)

# 🧾 Coberturas por parque
coberturas_por_parque <- list(
  aiguestortes   = c("AlpScrGr", "AlpConF"),
  cabaneros      = c("MedConF","MedScrub", "MedSclF","Grass","SMedMarF","Dehesa"),
  cabrera        = c("MedConF", "AridScr"),
  daimiel        = c("HydroRipVeg","Grass","Dehesa","HaloVeg","Crop", "WaterSurf"),
  donana         = c("MedConF", "MedScrub", "SaltMar", "HydroRipVeg", "DuneSand","WaterSurf"),
  guadarrama     = c("AlpScrGr", "MedScrub", "AlpConF", "Reforestation"),
  cies           = c("AtlScrub", "Grass", "Reforestation"),
  monfrague      = c("MedScrub", "MedSclF", "WaterSurf", "Reforestation", "Dehesa","Shrub"),
  ordesa         = c("AlpScrGr", "AlpConF", "TempDecF", "Grass"),
  picoseuropa    = c("AlpScrGr", "TempDecF", "Grass"),
  sierranevada   = c("AlpScrGr", "MedScrub", "AlpConF", "Grass")
)

# 📆 Etiquetas de mes
meses <- c("Jan","Feb","Mar","Apr","May","Jun",
           "Jul","Aug","Sep","Oct","Nov","Dec")

# 🎨 Colores por variable (sin SPEI)
colores_var <- c("pr" = "#003f5c", "tmax" = "#ef5675", "tmin" = "#7a5195")

# =============================
# 🔁 Loop por parques
# =============================
for (i in 1:nrow(parques_info)) {
  parque <- parques_info$parque[i]
  nombre_parque <- parques_info$nombre_parque[i]
  ncol <- parques_info$ncol[i]
  nrow <- parques_info$nrow[i]
  
  coberturas <- coberturas_por_parque[[parque]]
  
  df_sub <- df %>%
    filter(parque == !!parque, Id_SNVeg %in% coberturas) %>%
    mutate(
      Nombre_SNVeg = factor(nombres_vegetacion[Id_SNVeg],
                            levels = nombres_vegetacion[coberturas]),
      mes_label = factor(meses[mes], levels = meses)
    )
  
  # Crear subgráficos
  plots <- df_sub %>%
    split(.$Nombre_SNVeg) %>%
    lapply(function(subdf) {
      j <- which(levels(df_sub$Nombre_SNVeg) == unique(subdf$Nombre_SNVeg))
      es_primera_col <- ((j - 1) %% ncol == 0)
      
      ggplot(subdf, aes(x = mes_label, y = prop_sig, fill = variable)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_fill_manual(
          values = colores_var,
          name = "Climate Variable",
          labels = c("pr" = "Pr", "tmin" = "Tmin", "tmax" = "Tmax")
        ) +
        labs(
          x = NULL,
          y = if (es_primera_col) "% Significant Pixels" else NULL,
          title = unique(subdf$Nombre_SNVeg)
        ) +
        ylim(0, 100) +
        theme_minimal(base_size = 13) +
        theme(
          axis.title.y = element_text(size = 10),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(size = 12, hjust = 0)
        )
    })
  
  # Ensamblar figura
  figura <- wrap_plots(plots, ncol = ncol) +
    plot_annotation(
      title = nombre_parque,
      theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0))
    ) +
    plot_layout(guides = "collect") & 
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
  
  # 📐 Tamaño dinámico
  width  <- ncol * 4
  if (nrow < 2) {
    height <- nrow * 2.8
  } else {
    height <- nrow * 2.3
  }
  
  # 💾 Guardar figura
  ggsave(
    file.path(out_dir, sprintf("snv_response_negativo_%s.tif", parque)),
    figura, width = width, height = height, dpi = 300
  )
  ggsave(
    file.path(out_dir, sprintf("snv_response_negativo_%s.pdf", parque)),
    figura, width = width, height = height
  )
  ggsave(
    file.path(out_dir, sprintf("snv_response_negativo_%s.svg", parque)),
    figura, width = width, height = height
  )
  
  cat("✅ Figura guardada para", nombre_parque, "\n")
}
