# ==============================================================================
# FIGURAS S24-S67: CORRELACIONES POSITIVAS MÁXIMAS POR TIPO DE VEGETACIÓN,
# APILADAS POR ESCALA DE ACUMULACIÓN, UNA FIGURA POR PARQUE Y VARIABLE
# ==============================================================================
# Para cada parque y variable (SPEI, Pr, Tmax, Tmin), un panel por tipo de
# vegetación presente (coberturas_por_parque) con, para cada mes, el
# porcentaje de píxeles del tipo con correlación positiva significativa,
# apilado por la escala en que se alcanza la correlación máxima. Cada píxel
# cuenta una sola vez. Lee directamente los máximos con umbral de años y FDR
# (03_fdr/); los píxeles descartados por el umbral (significance NA) no
# cuentan ni en el numerador ni en el denominador.
#
# Los tipos sin ningún píxel significativo en todo el año se dibujan como
# panel vacío (es un resultado) y se listan por consola.
#
# Entradas: correlaciones_v2/maximos_cor_positive_fdr/<parque>_max_pos_correlation_<var>_<indice>[_selscales].nc
#           mapas_vegetacion/sistemasNaturales_vegetacion/SNV_parques/SNV_<parque>.shp
# Salida  : figuras/figuras_mensuales_apilados_absolutos_SNV_v2/fdr/<indice>/por_parque_pos/
#             snv_response_<var>_scales_absoluto_<parque>.{tif,pdf,svg}
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(svglite)

# =============================
# ⚙️ Parámetros
# =============================

indice <- "ndvi"    # "ndvi" o "kndvi"

DIR_MAX <- paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_positive_fdr")
variables <- c("tmax","tmin","pr","spei")

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

# 🗓️ Meses
meses <- c("Jan","Feb","Mar","Apr","May","Jun",
           "Jul","Aug","Sep","Oct","Nov","Dec")
n_meses <- length(meses)

# 🎨 Paletas por variable
paletas_variables <- list(
  spei = c("1"="#ffffb2","3"="#fecc5c","6"="#fd8d3c","9"="#f03b20","12"="#bd0026"),
  pr   = c("1"="#ffffcc","3"="#a1dab4","6"="#41b6c4","9"="#2c7fb8","12"="#253494"),
  tmax = c("1"="#ffffb2","3"="#fecc5c","6"="#fd8d3c","9"="#f03b20","12"="#bd0026"),
  tmin = c("1"="#edf8fb","3"="#b3cde3","6"="#8c96c6","9"="#8856a7","12"="#810f7c")
)

titulo_leyenda <- list(
  spei = "SPEI scale (months)",
  pr   = "Precipitation scale (months)",
  tmax = "Tmax scale (months)",
  tmin = "Tmin scale (months)"
)

# 🧾 Diccionario nombres vegetación
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
  Shrub="Shrubs",
  Dehesa="Dehesa",
  Grass="Mediterranean and Atlantic Grasslands",
  Reforestation="Reforestation",
  Crop="Crops",
  HaloVeg="Halophilous Vegetation",
  HydroRipVeg="Hydrophilous and Riparian Vegetation",
  SaltMar="Salt Marshes",
  DuneSand="Sand Dunes",
  BareArea="Bare area",
  WaterSurf="WaterSurf"
)

# 🧾 Coberturas válidas por parque
coberturas_por_parque <- list(
  aiguestortes   = c("AlpScrGr","AlpConF"),
  cabaneros      = c("MedConF","MedScrub","MedSclF","Grass","SMedMarF","Dehesa"),
  cabrera        = c("MedConF","AridScr"),
  daimiel        = c("HydroRipVeg","Grass","Dehesa","HaloVeg","Crop","WaterSurf"),
  donana         = c("MedConF","MedScrub","SaltMar","HydroRipVeg","DuneSand","WaterSurf"),
  guadarrama     = c("AlpScrGr","MedScrub","AlpConF","Reforestation"),
  cies           = c("AtlScrub","Grass","Reforestation"),
  monfrague      = c("MedScrub","MedSclF","WaterSurf","Reforestation","Dehesa","Shrub"),
  ordesa         = c("AlpScrGr","AlpConF","TempDecF","Grass"),
  picoseuropa    = c("AlpScrGr","TempDecF","Grass"),
  sierranevada   = c("AlpScrGr","MedScrub","AlpConF","Grass")
)

# =============================
# 🚀 Bucle principal
# =============================
out_dir <- sprintf(paste0(ROOT, "/articulo_2/figuras/figuras_mensuales_apilados_absolutos_SNV_v2/fdr/%s/por_parque_pos"), indice)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (variable in variables) {
  for (i in 1:nrow(parques_info)) {
    parque <- parques_info$parque[i]
    nombre_parque <- parques_info$nombre_parque[i]
    ncol <- parques_info$ncol[i]
    nrow <- parques_info$nrow[i]
    
    cat("\n📍 Procesando:", nombre_parque, "–", variable, "\n")
    
    # Archivo NetCDF
    if (variable == "spei") {
      nc_path <- sprintf("%s/%s_max_pos_correlation_%s_%s.nc",
                         DIR_MAX, parque, variable, indice)
    } else {
      nc_path <- sprintf("%s/%s_max_pos_correlation_%s_%s_selscales.nc",
                         DIR_MAX, parque, variable, indice)
    }
    
    snv_path <- sprintf(paste0(ROOT, "/mapas_vegetacion/sistemasNaturales_vegetacion/SNV_parques/SNV_%s.shp"), parque)
    
    # Saltar si falta el archivo
    if (!file.exists(nc_path) | !file.exists(snv_path)) {
      warning("❌ Falta archivo para ", parque, " - ", variable)
      next
    }
    
    # Leer NetCDF
    r_esc <- rast(nc_path, subds="scale")
    r_sig <- rast(nc_path, subds="significance")
    
    # Leer coberturas
    snv <- vect(snv_path)
    snv_proj <- project(snv, crs(r_esc))
    snv_proj$Id_SNVeg <- as.factor(snv_proj$Id_SNVeg)
    snv_proj$ID <- 1:nrow(snv_proj)
    
    # Cálculo
    resultados <- list()
    for (m in 1:n_meses) {
      escala_m <- r_esc[[m]]
      sig_m <- r_sig[[m]]
      names(escala_m) <- "scale"

      # Píxeles válidos: los descartados por el umbral de años tienen
      # significance NA aunque conserven scale
      escala_m_valida <- ifel(is.na(sig_m), NA, escala_m)
      names(escala_m_valida) <- "scale"
      
      df_total <- terra::extract(escala_m_valida, snv_proj) %>%
        left_join(as.data.frame(snv_proj)[,c("ID","Id_SNVeg")], by="ID") %>%
        filter(!is.na(scale),!is.na(Id_SNVeg)) %>%
        count(Id_SNVeg) %>% rename(n_total=n)
      
      # Píxeles significativos (ifel es seguro con NA en significance)
      escala_m <- ifel(is.na(sig_m) | sig_m != 2, NA, escala_m)
      names(escala_m) <- "scale"
      df_extract <- terra::extract(escala_m, snv_proj) %>%
        left_join(as.data.frame(snv_proj)[,c("ID","Id_SNVeg")], by="ID") %>%
        filter(!is.na(scale),!is.na(Id_SNVeg)) %>%
        mutate(escala=as.character(scale), mes=meses[m])
      
      # Rejilla completa tipo x escala del mes, con n = 0 donde no hay píxeles
      # significativos. Los niveles de escala se toman de la paleta para que
      # coincidan con los de scale_fill_manual.
      resumen <- df_total %>%
        crossing(escala = names(paletas_variables[[variable]])) %>%
        left_join(count(df_extract, Id_SNVeg, escala), by = c("Id_SNVeg", "escala")) %>%
        mutate(n = ifelse(is.na(n), 0L, n),
               prop = 100 * n / n_total,
               mes = meses[m]) %>%
        ungroup()
      resultados[[length(resultados)+1]] <- resumen
    }
    
    # DataFrame final
    coberturas_validas <- coberturas_por_parque[[parque]]
    df_final <- bind_rows(resultados) %>%
      filter(Id_SNVeg %in% coberturas_validas) %>%
      mutate(
        escala=factor(escala, levels=names(paletas_variables[[variable]])),
        mes=factor(mes, levels=meses),
        Nombre_SNVeg=factor(nombres_vegetacion[as.character(Id_SNVeg)],
                            levels=nombres_vegetacion[coberturas_validas])
      )
    
    # Crear plots
    # Rejilla completa tipo x mes x escala con ceros, para que todos los
    # paneles tengan el mismo eje X y la misma leyenda. split() sobre el factor
    # crea un grupo por nivel, incluidos los tipos sin ningún píxel
    # significativo, que se dibujan como panel vacío y se listan por consola.
    df_final <- complete(df_final, Nombre_SNVeg, mes, escala, fill = list(prop = 0))
    df_split <- split(df_final, df_final$Nombre_SNVeg)
    vacias <- names(which(tapply(df_final$prop, df_final$Nombre_SNVeg,
                                 function(v) all(is.na(v) | v == 0))))
    if (length(vacias))
      message("   sin ningun pixel significativo (panel vacio): ",
              paste(vacias, collapse = ", "))
    plots <- lapply(seq_along(df_split), function(i) {
        subdf <- df_split[[i]]
        es_primera_col <- ((i - 1) %% ncol == 0)
        
        ggplot(subdf, aes(x=mes, y=prop, fill=escala)) +
          geom_bar(stat="identity", position="stack") +
          scale_fill_manual(values=paletas_variables[[variable]],
                            name=titulo_leyenda[[variable]]) +
          labs(
            x=NULL,
            y=if (es_primera_col) "% Significant Pixels" else NULL,
            title=names(df_split)[i]
          ) +
          ylim(0,100) +
          theme_minimal(base_size=13) +
          theme(
            axis.text.x=element_text(angle=0,hjust=0.5),
            plot.title=element_text(size=12,hjust=0),
            axis.title.y=element_text(size=if (es_primera_col) 10 else 0)
          )
      })
    
    figura <- wrap_plots(plots, ncol=ncol) +
      plot_annotation(
        title=nombre_parque,
        theme=theme(plot.title=element_text(size=13, face="bold", hjust=0))
      ) +
      plot_layout(guides="collect") & 
      theme(legend.position="bottom",
            legend.text=element_text(size=10),
            legend.title=element_text(size=10))
    
    # 📐 Ajuste tamaño automático
    width <- ncol * 4
    if (nrow < 2) {
      height <- nrow * 2.8
    } else {
      height <- nrow * 2.3
    }
    
    # Guardar figura
    ggsave(sprintf("%s/snv_response_%s_scales_absoluto_%s.tif", out_dir, variable, parque),
           figura, width=width, height=height, dpi=300)
    ggsave(sprintf("%s/snv_response_%s_scales_absoluto_%s.pdf", out_dir, variable, parque),
           figura, width=width, height=height)
    ggsave(sprintf("%s/snv_response_%s_scales_absoluto_%s.svg", out_dir, variable, parque),
           figura, width=width, height=height)
    
  }
}
