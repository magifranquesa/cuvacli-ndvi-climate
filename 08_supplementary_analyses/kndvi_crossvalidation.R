# ==============================================================================
# VALIDACIÓN CRUZADA NDVI - kNDVI
# ==============================================================================
# Compara los máximos con FDR de los dos índices para las 96 combinaciones
# parque x variable x signo, sobre los píxeles válidos en ambos: correlación
# entre los perfiles mensuales de superficie significativa, correlación píxel a
# píxel de los coeficientes, acuerdo en signo, kappa de la decisión de
# significancia y acuerdo en la escala seleccionada.
#
# Entrada: correlaciones_v2/maximos_cor_{positive,negative}_fdr/ (ndvi y kndvi)
# Salida : resultados_v2/analisis_complementarios/xval.rds y tabla por consola
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")
OUT <- file.path(ROOT, "articulo_2/resultados_v2/analisis_complementarios"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressMessages(library(ncdf4))
options(width = 190)
B  <- paste0(ROOT, "/articulo_2/correlaciones_v2/")
PQ <- c("aiguestortes","cabaneros","cabrera","cies","daimiel","donana",
        "guadarrama","monfrague","ordesa","picoseuropa","sierranevada","sierranieves")
VA <- c("tmax","tmin","pr","spei")

abre <- function(signo, p, va, ix) {
  d <- sprintf("%smaximos_cor_%s_fdr/", B, if (signo=="pos") "positive" else "negative")
  s <- if (signo=="pos") "max_pos_correlation" else "max_neg_correlation"
  cand <- c(sprintf("%s%s_%s_%s_%s_selscales.nc", d,p,s,va,ix),
            sprintf("%s%s_%s_%s_%s.nc",           d,p,s,va,ix))
  f <- cand[file.exists(cand)][1]
  if (is.na(f)) return(NULL)
  nc_open(f)
}
mes <- function(nc, v, m) as.vector(ncvar_get(nc, v, start=c(1,1,m), count=c(-1,-1,1)))

res <- NULL
for (signo in c("pos","neg")) {
 cod <- if (signo=="pos") 2 else -2
 for (p in PQ) for (va in VA) {
  a <- abre(signo,p,va,"ndvi"); b <- abre(signo,p,va,"kndvi")
  if (is.null(a) || is.null(b)) { cat("falta:",signo,p,va,"\n"); next }
  pn <- pk <- numeric(12)
  R1 <- R2 <- S1 <- S2 <- E1 <- E2 <- c(); den <- 0
  for (m in 1:12) {
    r1 <- mes(a,"max_correlation",m); s1 <- mes(a,"significance",m); e1 <- mes(a,"scale",m)
    r2 <- mes(b,"max_correlation",m); s2 <- mes(b,"significance",m); e2 <- mes(b,"scale",m)
    ok <- !is.na(s1) & !is.na(s2)                 # pixel valido en LOS DOS indices
    if (m == 1) den <- sum(ok)
    d <- sum(ok); if (d == 0) next
    pn[m] <- 100*sum(s1[ok]==cod)/d; pk[m] <- 100*sum(s2[ok]==cod)/d
    R1 <- c(R1, r1[ok]); R2 <- c(R2, r2[ok])
    S1 <- c(S1, s1[ok]==cod); S2 <- c(S2, s2[ok]==cod)
    E1 <- c(E1, e1[ok]); E2 <- c(E2, e2[ok])
  }
  nc_close(a); nc_close(b)
  if (!length(R1)) next
  # kappa de Cohen sobre la decision de significacion
  po <- mean(S1 == S2); p1 <- mean(S1); p2 <- mean(S2)
  pe <- p1*p2 + (1-p1)*(1-p2); k <- if (pe < 1) (po-pe)/(1-pe) else NA
  amb <- S1 & S2
  res <- rbind(res, data.frame(signo, parque=p, variable=va,
    perfil_r = suppressWarnings(cor(pn, pk)),
    dif_media_pp = mean(abs(pn - pk)),
    med_ndvi = mean(pn), med_kndvi = mean(pk),
    px_r  = suppressWarnings(cor(R1, R2)),
    signo_ac = 100*mean(sign(R1) == sign(R2)),
    sig_ac = 100*po, kappa = k,
    escala_ac = if (sum(amb)) 100*mean(E1[amb] == E2[amb]) else NA_real_))
  rm(R1,R2,S1,S2,E1,E2); gc(verbose=FALSE)
 }
}
saveRDS(res, file.path(OUT, "xval.rds"))
cat("\n\n################ VALIDACION CRUZADA NDVI vs kNDVI ################\n")
for (sg in c("pos","neg")) {
  x <- res[res$signo==sg,]
  cat(sprintf("\n===== correlaciones %s =====\n", toupper(sg)))
  cat(sprintf("%-14s %-5s %7s %7s %7s %8s %8s %8s %7s %8s\n","parque","var",
      "perfilR","difpp","medNDVI","medkNDVI","pixelR","signo%","sig%","kappa"))
  for (i in 1:nrow(x)) with(x[i,], cat(sprintf("%-14s %-5s %7.3f %7.2f %7.1f %8.1f %8.3f %8.1f %7.1f %8.3f\n",
      parque, variable, perfil_r, dif_media_pp, med_ndvi, med_kndvi, px_r, signo_ac, sig_ac, kappa)))
  cat(sprintf("\n  MEDIANA  perfilR %.3f | dif %.2f pp | pixelR %.3f | signo %.1f%% | sig %.1f%% | kappa %.3f | escala %.1f%%\n",
      median(x$perfil_r,na.rm=TRUE), median(x$dif_media_pp), median(x$px_r,na.rm=TRUE),
      median(x$signo_ac), median(x$sig_ac), median(x$kappa,na.rm=TRUE), median(x$escala_ac,na.rm=TRUE)))
  cat(sprintf("  RANGO    perfilR %.3f-%.3f | pixelR %.3f-%.3f | kappa %.3f-%.3f\n",
      min(x$perfil_r,na.rm=TRUE), max(x$perfil_r,na.rm=TRUE),
      min(x$px_r,na.rm=TRUE), max(x$px_r,na.rm=TRUE),
      min(x$kappa,na.rm=TRUE), max(x$kappa,na.rm=TRUE)))
}
