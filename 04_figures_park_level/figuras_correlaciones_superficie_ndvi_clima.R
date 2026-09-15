# ============================================================
# FIGURAS 4-6: SUPERFICIE CON CORRELACIÓN SIGNIFICATIVA NDVI-CLIMA POR PARQUE
# ============================================================
# Una figura por variable (Tmax, Tmin, Pr) con todos los parques: para cada
# parque, dos paneles (correlaciones positivas a la izquierda, negativas a la
# derecha) con la superficie del parque que muestra correlación significativa
# en cada mes, apilada por la escala de acumulación (1, 3, 6, 9, 12 meses) en
# que se alcanza la correlación máxima.
#
# La superficie se expresa como fracción de los píxeles válidos del parque
# (con n >= 20 años en al menos un mes). Las etiquetas sobre las barras dan el
# total mensual.
#
# Cambia "indice" a "kndvi" para generar la versión con kNDVI.
#
# Entrada: correlaciones_v2/maximos_cor_{positive,negative}_fdr/
#            <parque>_max_{pos,neg}_correlation_<var>_<indice>_selscales.nc
# Salida : figuras/figuras_mensuales_v2/fdr/cor_combined/
#            figura_mensual_pos_neg_<var>_<indice>_v2.{tif,pdf,svg}
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

variables <- c("tmax", "tmin", "pr")
indice <- "ndvi"

# Máximos con umbral de años válidos y FDR (03_fdr/)
base_dir_pos <- paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_positive_fdr/")
base_dir_neg <- paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_negative_fdr/")
out_dir_comb <- paste0(ROOT, "/articulo_2/figuras/figuras_mensuales_v2/fdr/cor_combined/")
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

# 🔠 Etiquetas de variables e índice
etiquetas_variables <- c(tmax = "Tmax", tmin = "Tmin", pr = "Pr")
etiquetas_indices   <- c(ndvi = "NDVI", kndvi = "kNDVI")

# 🎨 Paletas por variable
paletas_variables <- list(
  tmax = c("1" = "#ffffb2", "3" = "#fecc5c", "6" = "#fd8d3c", "9" = "#f03b20", "12" = "#bd0026"),
  tmin = c("1" = "#edf8fb", "3" = "#b3cde3", "6" = "#8c96c6", "9" = "#8856a7", "12" = "#810f7c"),
  pr   = c("1" = "#ffffcc", "3" = "#a1dab4", "6" = "#41b6c4", "9" = "#2c7fb8", "12" = "#253494")
)

niveles_meses  <- factor(month.abb, levels = month.abb)
niveles_escala <- factor(c(1,3,6,9,12), levels = c(1,3,6,9,12))

# =============================
# 🧩 Función de panel
# =============================
panel_superficie <- function(parque, var, signo = c("pos","neg"),
                             base_dir, colores_escala,
                             titulo_y, mostrar_titulo = TRUE, mostrar_y = TRUE) {
  signo <- match.arg(signo)
  sufijo <- if (signo == "pos") "max_pos_correlation" else "max_neg_correlation"
  archivo_nc <- file.path(base_dir, paste0(parque, "_", sufijo, "_", var, "_", tolower(indice), "_selscales.nc"))
  
  if (!file.exists(archivo_nc)) {
    warning("Archivo no encontrado: ", archivo_nc)
    df_empty <- data.frame(mes = niveles_meses, escala = niveles_escala[1], pct = 0)
    return(
      ggplot(df_empty, aes(x = mes, y = pct, fill = escala)) +
        geom_col() +
        scale_y_continuous(labels = if (mostrar_y) scales::percent_format(accuracy = 1) else NULL,
                           limits = c(0,1)) +
        scale_fill_manual(values = colores_escala, drop = FALSE) +
        labs(title = if (mostrar_titulo) nombre_parques[parque] else NULL,
             x = NULL, y = if (mostrar_y) titulo_y else NULL) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(face = "bold"),
          axis.text.y = if (!mostrar_y) element_blank() else element_text(size = 7),
          axis.ticks.y = if (!mostrar_y) element_blank() else element_line(),
          legend.position = "none",
          panel.border = element_rect(colour = if (signo == "pos") "#2E7D32" else "#C62828",
                                      fill = NA, linewidth = 0.6)
        )
    )
  }
  
  nc  <- nc_open(archivo_nc)
  cor <- ncvar_get(nc, "max_correlation")
  sig <- ncvar_get(nc, "significance")
  esc <- ncvar_get(nc, "scale")
  nc_close(nc)
  
  dims <- dim(cor)
  # Denominador: píxeles válidos del parque (significance no NA en algún mes).
  # Los descartados por n < 20 conservan max_correlation pero tienen
  # significance = NA, y no cuentan ni en el numerador ni en el denominador.
  total_pixeles <- sum(!apply(sig, c(1,2), function(x) all(is.na(x))))
  
  # Máscara de significancia a prueba de NA: donde significance es NA la
  # condición sig == cod_sig también lo es, y una asignación con índice NA
  # dejaría la celda intacta.
  cod_sig  <- if (signo == "pos") 2 else -2
  mask_sig <- !is.na(sig) & sig == cod_sig
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
  
  borde_col <- if (signo == "pos") "#2E7D32" else "#C62828"
  
  p <- ggplot(df_escala, aes(x = mes, y = pct, fill = escala)) +
    geom_col() +
    scale_fill_manual(values = colores_escala, drop = FALSE) +
    scale_y_continuous(labels = if (mostrar_y) scales::percent_format(accuracy = 1) else NULL,
                       limits = c(0, 1)) +
    labs(
      title = if (mostrar_titulo) nombre_parques[parque] else NULL,
      x = NULL,
      y = if (mostrar_y) titulo_y else NULL,
      fill = paste0(var_label, " scale (months)")
    )+
    geom_text(
      data = df_total_mes,
      aes(x = mes, y = pmin(pct_mes, 0.98), 
          label = paste0(round(100*pct_mes), "%")),
      inherit.aes = FALSE, size = 3, vjust = -0.2,
      colour = "black"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title  = element_text(face = "bold"),
      axis.text.x = element_text(size = 7),
      axis.text.y = if (!mostrar_y) element_blank() else element_text(size = 7),
      axis.ticks.y = if (!mostrar_y) element_blank() else element_line(),
      legend.position = "none",
      panel.border = element_rect(colour = borde_col, fill = NA, linewidth = 0.7)
    )
  
  return(p)
}

# =============================
# 🧱 Construcción de figuras
# =============================
# =============================
# 🧱 Construcción de figura única (4 columnas: Pos | Neg | Pos | Neg …)
# =============================
for (var in variables) {
  colores_escala <- paletas_variables[[var]]
  var_label      <- etiquetas_variables[[var]]
  indice_label   <- etiquetas_indices[[indice]]
  titulo_y       <- paste0(var_label, " ~ ", indice_label)
  
  # Lista intercalada de paneles
  paneles <- flatten(
    map(seq_along(parques), function(i) {
      parque <- parques[i]
      list(
        panel_superficie(parque, var, "pos", base_dir_pos, colores_escala, titulo_y,
                         mostrar_titulo = TRUE, mostrar_y = TRUE),
        panel_superficie(parque, var, "neg", base_dir_neg, colores_escala, titulo_y,
                         mostrar_titulo = FALSE, mostrar_y = FALSE)
      )
    })
  )
  
  # Títulos de columnas
  pos_title <- ggplot() + theme_void() +
    ggtitle("Positive correlations") +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11)) +
    guides(fill = "none")
  
  neg_title <- ggplot() + theme_void() +
    ggtitle("Negative correlations") +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11)) +
    guides(fill = "none")
  
  fila_titulos <- wrap_plots(pos_title, neg_title, pos_title, neg_title, ncol = 4)
  
  # Paneles
  figura_paneles <- wrap_plots(paneles, ncol = 4)
  
  # Combinar títulos + paneles (ajustar altura fila de títulos)
  figura_comb <- fila_titulos / figura_paneles +
    plot_layout(heights = c(0.01, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal")
  
  # Guardar
  altura_por_parque <- 0.9
  width <- 14
  height <- altura_por_parque * length(parques) + 1
  
  ggsave(file.path(out_dir_comb, paste0("figura_mensual_pos_neg_", var, "_", indice, "_v2.tif")),
         figura_comb, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir_comb, paste0("figura_mensual_pos_neg_", var, "_", indice, "_v2.pdf")),
         figura_comb, width = width, height = height)
  ggsave(file.path(out_dir_comb, paste0("figura_mensual_pos_neg_", var, "_", indice, "_v2.svg")),
         figura_comb, width = width, height = height)
}