# ============================================================
# FIGURA 7: SUPERFICIE CON CORRELACIÓN POSITIVA SIGNIFICATIVA NDVI-SPEI POR PARQUE
# ============================================================
# Un panel por parque con la superficie que muestra correlación positiva
# significativa con el SPEI en cada mes, apilada por la escala del SPEI
# (1, 3, 6, 9, 12 meses) en que se alcanza la correlación máxima.
#
# La superficie se expresa como fracción de los píxeles válidos del parque
# (con n >= 20 años en al menos un mes). Las etiquetas sobre las barras dan el
# total mensual.
#
# Cambia "indice" a "kndvi" para generar la versión con kNDVI.
#
# Entrada: correlaciones_v2/maximos_cor_positive_fdr/<parque>_max_pos_correlation_spei_<indice>.nc
# Salida : figuras/figuras_mensuales_v2/fdr/cor_spei_pos/
# ============================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(svglite)
library(scales)  # percent_format

# =============================
# ⚙️ Parámetros
# =============================
parques <- c("aiguestortes", "cabaneros", "cabrera", "daimiel", "donana",
             "guadarrama", "cies", "monfrague", "ordesa", "picoseuropa",
             "sierranevada", "sierranieves")

variables <- c("spei")  # solo SPEI
indice <- "ndvi"

# Máximos con umbral de años válidos y FDR (03_fdr/)
base_dir_pos <- paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_positive_fdr/")
out_dir_comb <- paste0(ROOT, "/articulo_2/figuras/figuras_mensuales_v2/fdr/cor_spei_pos/")
dir.create(out_dir_comb, showWarnings = FALSE, recursive = TRUE)

# 🏞️ Nombres de parques legibles
nombre_parques <- c(
  "aiguestortes"  = "Aigüestortes",     "cabaneros"     = "Cabañeros",
  "cabrera"       = "Cabrera",          "daimiel"       = "Tablas de Daimiel",
  "donana"        = "Doñana",           "guadarrama"    = "Sierra de Guadarrama",
  "cies"          = "Islas Atlánticas", "monfrague"     = "Monfragüe",
  "ordesa"        = "Ordesa y Monte Perdido", "picoseuropa"   = "Picos de Europa",
  "sierranevada"  = "Sierra Nevada",    "sierranieves"  = "Sierra de las Nieves"
)

# 🔠 Etiquetas
etiquetas_variables <- c(spei = "SPEI")
etiquetas_indices   <- c(ndvi = "NDVI", kndvi = "kNDVI")

# 🎨 Paleta para SPEI = misma que Tmax
paletas_variables <- list(
  spei = c("1" = "#ffffb2", "3" = "#fecc5c", "6" = "#fd8d3c", "9" = "#f03b20", "12" = "#bd0026")
)

niveles_meses  <- factor(month.abb, levels = month.abb)
niveles_escala <- factor(c(1,3,6,9,12), levels = c(1,3,6,9,12))

# =============================
# 🧩 Función de panel (solo barras, altura = % superficie)
# =============================
panel_superficie_spei_pos <- function(parque, colores_escala, titulo_y,
                                      mostrar_titulo = TRUE, mostrar_y = TRUE) {
  # patrón: {parque}_max_pos_correlation_spei_ndvi.nc
  archivo_nc <- file.path(base_dir_pos, paste0(parque, "_max_pos_correlation_spei_", tolower(indice), ".nc"))
  
  if (!file.exists(archivo_nc)) {
    warning("Archivo no encontrado: ", archivo_nc)
    df_empty <- data.frame(mes = niveles_meses, escala = niveles_escala[1], pct = 0)
    return(
      ggplot(df_empty, aes(x = mes, y = pct, fill = escala)) +
        geom_col() +
        scale_y_continuous(labels = if (mostrar_y) percent_format(accuracy = 1) else NULL,
                           limits = c(0,1)) +
        scale_fill_manual(values = colores_escala, drop = FALSE) +
        labs(title = if (mostrar_titulo) nombre_parques[parque] else NULL,
             x = NULL, y = if (mostrar_y) titulo_y else NULL) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(face = "bold"),
          axis.text.y = if (!mostrar_y) element_blank() else element_text(size = 7),
          axis.ticks.y = if (!mostrar_y) element_blank() else element_line(),
          legend.position = "none"
        )
    )
  }
  
  nc  <- nc_open(archivo_nc)
  cor <- ncvar_get(nc, "max_correlation")
  sig <- ncvar_get(nc, "significance")
  esc <- ncvar_get(nc, "scale")
  nc_close(nc)
  
  dims <- dim(cor)  # lon x lat x mes
  # Denominador: píxeles válidos del parque (significance no NA en algún mes).
  # Los descartados por n < 20 conservan max_correlation pero tienen
  # significance = NA, y no cuentan ni en el numerador ni en el denominador.
  total_pixeles <- sum(!apply(sig, c(1,2), function(x) all(is.na(x))))
  
  # Solo significativas positivas, con máscara a prueba de NA (una asignación
  # con índice NA dejaría la celda intacta)
  mask_sig <- !is.na(sig) & sig == 2
  cor[!mask_sig] <- NA
  esc[!mask_sig] <- NA
  
  df <- expand.grid(lon = 1:dims[1], lat = 1:dims[2], mes = 1:12)
  df$cor    <- as.vector(cor)
  df$escala <- as.vector(esc)
  
  df <- df %>%
    filter(!is.na(cor)) %>%
    mutate(
      mes    = factor(month.abb[mes], levels = month.abb),
      escala = factor(escala, levels = levels(niveles_escala))
    )
  
  df_escala <- df %>%
    count(mes, escala, name = "n") %>%
    complete(mes = niveles_meses, escala = niveles_escala, fill = list(n = 0)) %>%
    mutate(pct = n / total_pixeles)
  
  df_total_mes <- df_escala %>%
    group_by(mes) %>% summarise(pct_mes = sum(pct), .groups = "drop")
  
  ggplot(df_escala, aes(x = mes, y = pct, fill = escala)) +
    geom_col() +
    scale_fill_manual(values = colores_escala, drop = FALSE) +
    scale_y_continuous(labels = if (mostrar_y) percent_format(accuracy = 1) else NULL,
                       limits = c(0, 1)) +
    labs(
      title = if (mostrar_titulo) nombre_parques[parque] else NULL,
      x = NULL,
      y = if (mostrar_y) titulo_y else NULL,
      fill = "Scale (months)"
    ) +
    geom_text(
      data = df_total_mes,
      aes(x = mes, y = pmin(pct_mes, 0.98), label = paste0(round(100*pct_mes), "%")),
      inherit.aes = FALSE, size = 3, vjust = -0.2, colour = "black"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title  = element_text(face = "bold"),
      axis.text.x = element_text(size = 7),
      axis.text.y = if (!mostrar_y) element_blank() else element_text(size = 7),
      axis.ticks.y = if (!mostrar_y) element_blank() else element_line(),
      legend.title  = element_text(size = 10),
      legend.text   = element_text(size = 9),
      legend.key.size = unit(0.6, "cm"),
      legend.title.align = 0.5,
      legend.margin = margin(t = 3)
    )
}

# =============================
# 🧱 Construcción de figura única (ncol = 4 por fila)
# =============================
ncol_paneles <- 4  # cambia a 3 si lo prefieres
for (var in variables) {
  colores_escala <- paletas_variables[[var]]
  var_label      <- etiquetas_variables[[var]]
  indice_label   <- etiquetas_indices[[indice]]
  titulo_y       <- paste0(var_label, " ~ ", indice_label)  # "SPEI ~ NDVI"
  
  paneles <- map(seq_along(parques), function(i) {
    # Mostrar eje Y solo en la primera columna de cada fila
    mostrar_y <- ((i - 1) %% ncol_paneles) == 0
    mostrar_titulo <- TRUE  # título en todos (puedes ponerlo solo en los que quieras)
    panel_superficie_spei_pos(
      parque = parques[i],
      colores_escala = colores_escala,
      titulo_y = titulo_y,
      mostrar_titulo = mostrar_titulo,
      mostrar_y = mostrar_y
    )
  })
  
  figura_comb <- wrap_plots(paneles, ncol = ncol_paneles) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal")
  
  # Medidas: más compacto al tener 4 por fila
  altura_por_fila <- 0.9  # ajusta si lo quieres más o menos alto
  num_filas <- ceiling(length(parques) / ncol_paneles)
  width  <- 14
  height <- altura_por_fila * num_filas
  
  ggsave(file.path(out_dir_comb, paste0("figura_mensual_pos_", var, "_", indice, "_SPEI.tif")),
         figura_comb, width = width, height = 6, dpi = 300)
  ggsave(file.path(out_dir_comb, paste0("figura_mensual_pos_", var, "_", indice, "_SPEI.pdf")),
         figura_comb, width = width, height = 6)
  ggsave(file.path(out_dir_comb, paste0("figura_mensual_pos_", var, "_", indice, "_SPEI.svg")),
         figura_comb, width = width, height = 6)
}
