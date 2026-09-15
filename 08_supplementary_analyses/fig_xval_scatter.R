# ==============================================================================
# FIGURA: DISPERSIÓN PÍXEL A PÍXEL DEL COEFICIENTE NDVI FRENTE A kNDVI
# ==============================================================================
# Mismas doce combinaciones parque-variable y mismo orden que fig_xval_kndvi.R.
# Se representa la correlación máxima sobre las cinco escalas de acumulación,
# con todos los píxeles válidos en los dos índices, sean significativos o no.
#
# Entrada: correlaciones_v2/maximos_cor_positive_fdr/ (ndvi y kndvi)
# Salida : figuras/complementarias/fig_xval_scatter.{png,tif}
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

suppressMessages({library(ncdf4); library(ggplot2); library(dplyr); library(patchwork)})

B   <- paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_positive_fdr/")
OUT <- paste0(ROOT, "/articulo_2/figuras/complementarias/")
N_SUB <- 8000   # pixeles por mes y combinacion (12 meses -> ~96.000 puntos/panel)
set.seed(42)

combos <- data.frame(
  parque = c("monfrague","monfrague","cabrera","cabrera","cabaneros","cabaneros",
             "donana","donana","daimiel","daimiel","sierranevada","sierranevada"),
  var    = c("spei","pr","spei","pr","spei","pr","spei","pr","spei","pr","spei","pr"),
  stringsAsFactors = FALSE)
bonito  <- c(monfrague="Monfragüe", cabrera="Cabrera", cabaneros="Cabañeros",
             donana="Doñana", daimiel="Daimiel", sierranevada="Sierra Nevada")
vbonito <- c(spei="SPEI", pr="Pr")

abre <- function(p, v, ix) {
  cand <- c(sprintf("%s%s_max_pos_correlation_%s_%s_selscales.nc", B, p, v, ix),
            sprintf("%s%s_max_pos_correlation_%s_%s.nc",           B, p, v, ix))
  nc_open(cand[file.exists(cand)][1])
}
mes <- function(nc, v, m) as.vector(ncvar_get(nc, v, start = c(1,1,m), count = c(-1,-1,1)))

paneles <- list()
for (k in seq_len(nrow(combos))) {
  p <- combos$parque[k]; v <- combos$var[k]
  a <- abre(p, v, "ndvi"); b <- abre(p, v, "kndvi")
  X <- Y <- numeric(0)
  for (m in 1:12) {
    r1 <- mes(a,"max_correlation",m); s1 <- mes(a,"significance",m)
    r2 <- mes(b,"max_correlation",m); s2 <- mes(b,"significance",m)
    ok <- which(!is.na(s1) & !is.na(s2))
    if (!length(ok)) next
    id <- if (length(ok) > N_SUB) sample(ok, N_SUB) else ok
    X <- c(X, r1[id]); Y <- c(Y, r2[id])
  }
  nc_close(a); nc_close(b)
  rr <- cor(X, Y)
  bajo <- 100 * mean(Y < X)
  d <- data.frame(x = X, y = Y)
  paneles[[k]] <- ggplot(d, aes(x, y)) +
    geom_hex(bins = 55) +
    scale_fill_gradient(low = "#dbe7f0", high = "#123f63", trans = "log10", guide = "none") +
    geom_abline(slope = 1, intercept = 0, colour = "#c0392b", linewidth = 0.45) +
    coord_fixed(xlim = c(-0.4, 1), ylim = c(-0.4, 1), expand = FALSE) +
    scale_x_continuous(breaks = c(-0.25, 0.25, 0.75)) +
    scale_y_continuous(breaks = c(-0.25, 0.25, 0.75)) +
    labs(title = sprintf("%s — %s", bonito[[p]], vbonito[[v]]),
         subtitle = sprintf("r = %.2f   %.0f %% below 1:1", rr, bajo),
         x = if (k > 8) "r (NDVI)" else NULL,
         y = if (k %% 4 == 1) "r (kNDVI)" else NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 7.2, colour = "grey30"),
          axis.title = element_text(size = 8, colour = "grey30"),
          panel.grid.minor = element_blank())
  rm(X, Y, d); gc(verbose = FALSE)
  cat(sprintf("  %-14s %-5s  r = %.3f   %.0f %% por debajo de 1:1\n", p, v, rr, bajo))
}

fig <- wrap_plots(paneles, ncol = 4) +
  plot_annotation(
    title = "NDVI vs kNDVI: pixel-wise maximum correlation coefficients",
    subtitle = "Same twelve park–variable combinations as in the profile figure. Red line is 1:1. Darker hexagons contain more pixels (log scale).",
    caption = "All pixels valid in both indices, significant or not; maximum over the five accumulation timescales. Points below the 1:1 line are pixels where kNDVI gives a weaker correlation than NDVI.",
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 8.6, colour = "grey30"),
                  plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0)))

ggsave(paste0(OUT, "fig_xval_scatter.png"), fig, width = 10.5, height = 9.4, dpi = 200, bg = "white")
ggsave(paste0(OUT, "fig_xval_scatter.tif"), fig, width = 10.5, height = 9.4, dpi = 300, compression = "lzw", bg = "white")
cat("\nguardado en", OUT, "\n")
