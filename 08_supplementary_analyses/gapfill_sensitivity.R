# ==============================================================================
# SENSIBILIDAD AL RELLENO FOCAL DE LA REJILLA CLIMÁTICA
# ==============================================================================
# En los dos parques con mayor superficie rellenada (Cabrera e Islas
# Atlánticas), compara la superficie con correlación positiva significativa
# mensual calculada con todos los píxeles frente a la calculada solo con los
# píxeles cuyo clima es nativo (sin relleno focal). Requiere la rejilla
# climática de 1,1 km para reconstruir la máscara de celdas nativas.
#
# Entrada: correlaciones_v2/maximos_cor_positive_fdr/, CLIMA/agregados_lag/pr_scale_1.nc
# Salida : resultados_v2/analisis_complementarios/gapsens.rds y tabla por consola
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")
OUT <- file.path(ROOT, "articulo_2/resultados_v2/analisis_complementarios"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressMessages({library(terra); library(ncdf4)})
terraOptions(progress = 0); options(width = 140)
path <- paste0(ROOT, '/CLIMA/agregados_lag/')
vip  <- paste0(ROOT, '/MOSAICOS/mosaicos_all/')
shp  <- paste0(ROOT, '/limites_red/desc_Red_PN_LIM_Enero2023/Limites_PN.shp')
FDR  <- paste0(ROOT, '/articulo_2/correlaciones_v2/maximos_cor_positive_fdr/')
pv_all <- project(vect(shp), "EPSG:4326")

res <- NULL
for (parque in c("cabrera","cies")) {
  vi_r <- rast(paste0(vip, parque, "_NDVI_filled_clean.nc"))[[1]]
  pv   <- project(pv_all[pv_all$name_pn == parque, ], crs(vi_r))
  vi_r <- mask(vi_r, pv, touches = FALSE)
  cs   <- rast(paste0(path, "pr_scale_1.nc")); crs(cs) <- "EPSG:23030"
  s    <- cs[[277:756]]
  nat_r <- resample(mask(project(s, crs(vi_r)), pv, touches = TRUE), vi_r, "bilinear")[[1]]
  nat_r <- !is.na(nat_r)                       # TRUE donde el clima es NATIVO

  for (v in c("pr","spei","tmax","tmin")) {
    cand <- c(sprintf("%s%s_max_pos_correlation_%s_ndvi_selscales.nc", FDR, parque, v),
              sprintf("%s%s_max_pos_correlation_%s_ndvi.nc",           FDR, parque, v))
    f <- cand[file.exists(cand)][1]; nc <- nc_open(f)
    lon <- nc$dim$lon$vals; lat <- nc$dim$lat$vals
    pts <- as.matrix(expand.grid(lon = lon, lat = lat))          # orden [lon, lat] = el del array
    nat <- terra::extract(nat_r, pts)[,1]; nat[is.na(nat)] <- FALSE
    todo <- nativo <- numeric(12)
    for (m in 1:12) {
      sg <- as.vector(ncvar_get(nc, "significance", start=c(1,1,m), count=c(-1,-1,1)))
      ok <- !is.na(sg)
      todo[m]   <- 100*sum(sg[ok] == 2)/sum(ok)
      kn <- ok & nat
      nativo[m] <- if (sum(kn)) 100*sum(sg[kn] == 2)/sum(kn) else NA
    }
    nc_close(nc)
    res <- rbind(res, data.frame(parque, var=v, mes=1:12, todo, nativo))
  }
}
saveRDS(res, file.path(OUT, "gapsens.rds"))
cat("=== SENSIBILIDAD AL GAP-FILLING, datos v2 con FDR ===\n")
cat("    % de superficie significativa: todos los pixeles vs solo clima nativo\n\n")
cat(sprintf("%-10s %-6s %10s %10s %9s %9s %9s\n","parque","var","media tod","media nat","dif media","dif mediana","r perfil"))
for (p in unique(res$parque)) for (v in c("pr","spei","tmax","tmin")) {
  x <- res[res$parque==p & res$var==v, ]
  d <- x$nativo - x$todo
  cat(sprintf("%-10s %-6s %10.2f %10.2f %+9.2f %+9.2f %9.4f\n", p, v,
      mean(x$todo), mean(x$nativo,na.rm=TRUE), mean(d,na.rm=TRUE), median(d,na.rm=TRUE),
      suppressWarnings(cor(x$todo, x$nativo, use="complete.obs"))))
}
d <- res$nativo - res$todo
cat(sprintf("\nGLOBAL: diferencia mediana absoluta %.2f pp | percentil 95 %.2f pp | maxima %.2f pp\n",
    median(abs(d), na.rm=TRUE), quantile(abs(d), .95, na.rm=TRUE), max(abs(d), na.rm=TRUE)))
