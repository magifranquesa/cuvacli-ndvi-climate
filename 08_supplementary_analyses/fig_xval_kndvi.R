# ==============================================================================
# FIGURA: PERFILES MENSUALES NDVI FRENTE A kNDVI
# ==============================================================================
# Doce paneles, uno por cada combinación parque-variable en la que el NDVI da
# en media anual más del 20 % de superficie con correlación positiva
# significativa. En cada panel, el perfil mensual de superficie significativa
# con NDVI y con kNDVI, y los estadísticos de acuerdo de esa combinación
# (correlación de perfiles, correlación píxel a píxel, kappa).
#
# Entrada: correlaciones_v2/maximos_cor_positive_fdr/ (ndvi y kndvi)
# Salida : figuras/complementarias/fig_xval_kndvi.{png,tif}
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

suppressMessages({library(ncdf4); library(ggplot2); library(dplyr); library(patchwork)})

B   <- paste0(ROOT, "/articulo_2/correlaciones_v2/maximos_cor_positive_fdr/")
OUT <- paste0(ROOT, "/articulo_2/figuras/complementarias/")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

combos <- data.frame(
  parque = c("monfrague","monfrague","cabrera","cabrera","cabaneros","cabaneros",
             "donana","donana","daimiel","daimiel","sierranevada","sierranevada"),
  var    = c("spei","pr","spei","pr","spei","pr","spei","pr","spei","pr","spei","pr"),
  stringsAsFactors = FALSE)
bonito <- c(monfrague="Monfragüe", cabrera="Cabrera", cabaneros="Cabañeros",
            donana="Doñana", daimiel="Daimiel", sierranevada="Sierra Nevada")
vbonito <- c(spei="SPEI", pr="Pr")

abre <- function(p, v, ix) {
  cand <- c(sprintf("%s%s_max_pos_correlation_%s_%s_selscales.nc", B, p, v, ix),
            sprintf("%s%s_max_pos_correlation_%s_%s.nc",           B, p, v, ix))
  nc_open(cand[file.exists(cand)][1])
}
mes <- function(nc, v, m) as.vector(ncvar_get(nc, v, start = c(1,1,m), count = c(-1,-1,1)))

paneles <- list(); resumen <- NULL
for (k in seq_len(nrow(combos))) {
  p <- combos$parque[k]; v <- combos$var[k]
  a <- abre(p, v, "ndvi"); b <- abre(p, v, "kndvi")
  pn <- pk <- numeric(12); R1 <- R2 <- S1 <- S2 <- c()
  for (m in 1:12) {
    r1 <- mes(a,"max_correlation",m); s1 <- mes(a,"significance",m)
    r2 <- mes(b,"max_correlation",m); s2 <- mes(b,"significance",m)
    ok <- !is.na(s1) & !is.na(s2); d <- sum(ok)
    pn[m] <- 100*sum(s1[ok] == 2)/d; pk[m] <- 100*sum(s2[ok] == 2)/d
    id <- sample.int(d, min(d, 12000))            # submuestra para el pixel-r
    R1 <- c(R1, r1[ok][id]); R2 <- c(R2, r2[ok][id])
    S1 <- c(S1, s1[ok] == 2); S2 <- c(S2, s2[ok] == 2)
  }
  nc_close(a); nc_close(b)
  po <- mean(S1 == S2); q1 <- mean(S1); q2 <- mean(S2)
  pe <- q1*q2 + (1-q1)*(1-q2); kap <- (po - pe)/(1 - pe)
  rp <- cor(pn, pk); rx <- cor(R1, R2)
  resumen <- rbind(resumen, data.frame(parque = p, var = v, perfil = rp, pixel = rx, kappa = kap))

  d <- rbind(data.frame(mes = 1:12, pct = pn, idx = "NDVI"),
             data.frame(mes = 1:12, pct = pk, idx = "kNDVI"))
  d$idx <- factor(d$idx, levels = c("NDVI","kNDVI"))
  paneles[[k]] <- ggplot(d, aes(mes, pct, colour = idx, linetype = idx)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
    scale_colour_manual(values = c(NDVI = "#1b6ca8", kNDVI = "#e07b39"), name = NULL) +
    scale_linetype_manual(values = c(NDVI = "solid", kNDVI = "22"), name = NULL) +
    scale_x_continuous(breaks = c(1,4,7,10), labels = c("Jan","Apr","Jul","Oct")) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = sprintf("%s — %s", bonito[[p]], vbonito[[v]]),
         subtitle = sprintf("profile r = %.2f   pixel r = %.2f   \u03BA = %.2f", rp, rx, kap),
         x = NULL, y = if (k %% 4 == 1) "% of park area" else NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 7.2, colour = "grey30"),
          axis.title.y = element_text(size = 8, colour = "grey30"),
          legend.position = "none", panel.grid.minor = element_blank())
  rm(R1, R2, S1, S2); gc(verbose = FALSE)
}

leyenda <- paneles[[1]] + theme(legend.position = "bottom",
                                legend.text = element_text(size = 9))
fig <- wrap_plots(paneles, ncol = 4) +
  plot_annotation(
    title = "NDVI vs kNDVI: monthly percentage of park area with significant positive correlation",
    subtitle = "The twelve park–variable combinations exceeding 20 % of significant area on annual average",
    caption = "Agreement statistics computed within each combination: profile r over the 12 monthly values; pixel r between the correlation coefficients; \u03BA on the significance decision.",
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey30"),
                  plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(paste0(OUT, "fig_xval_kndvi.png"), fig, width = 11, height = 7.6, dpi = 200, bg = "white")
ggsave(paste0(OUT, "fig_xval_kndvi.tif"), fig, width = 11, height = 7.6, dpi = 300, compression = "lzw", bg = "white")
cat("\n=== resumen de las 12 combinaciones ===\n")
print(within(resumen, {perfil <- round(perfil,3); pixel <- round(pixel,3); kappa <- round(kappa,3)}), row.names = FALSE)
cat("\nguardado en", OUT, "\n")
