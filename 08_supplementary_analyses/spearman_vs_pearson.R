# ==============================================================================
# PEARSON FRENTE A SPEARMAN
# ==============================================================================
# Compara los dos coeficientes en la submuestra reproducible de píxeles en la
# que 01_correlations/ calcula Spearman (5 %), para los 12 parques, 4 variables,
# 5 escalas y 12 meses. A cada coeficiente se le aplica el umbral de años y el
# FDR por campo parque-mes-escala, y se mide: la correlación entre ambos
# coeficientes, el acuerdo en signo, el acuerdo en la decisión de significancia
# (kappa de Cohen) y la superficie significativa con cada uno.
#
# Entrada: correlaciones_v2/<parque>/correlation_ndvi_<var>_scale_<s>_<parque>.nc
# Salida : resultados_v2/analisis_complementarios/pvs.rds y resumen por consola
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")
OUT <- file.path(ROOT, "articulo_2/resultados_v2/analisis_complementarios"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressMessages(library(ncdf4))
options(width = 140)
B  <- paste0(ROOT, "/articulo_2/correlaciones_v2/")
PQ <- c("aiguestortes","cabaneros","cabrera","cies","daimiel","donana",
        "guadarrama","monfrague","ordesa","picoseuropa","sierranevada","sierranieves")
VA <- c("tmax","tmin","pr","spei"); ES <- c(1,3,6,9,12); Q <- 0.05; NMIN <- 20

res <- NULL
for (p in PQ) for (v in VA) for (e in ES) {
  f <- sprintf("%s%s/correlation_ndvi_%s_scale_%d_%s.nc", B, p, v, e, p)
  if (!file.exists(f)) next
  nc <- nc_open(f)
  rp <- rs <- pp <- ps <- nn <- vector("list", 12)
  for (m in 1:12) {
    st <- c(1,1,m); ct <- c(-1,-1,1)
    s <- as.vector(ncvar_get(nc, "correlation_spearman", start=st, count=ct))
    k <- which(!is.na(s))                       # la submuestra del 5 %
    if (!length(k)) next
    rs[[m]] <- s[k]
    rp[[m]] <- as.vector(ncvar_get(nc, "correlation",       start=st, count=ct))[k]
    pp[[m]] <- as.vector(ncvar_get(nc, "p_value",           start=st, count=ct))[k]
    ps[[m]] <- as.vector(ncvar_get(nc, "p_value_spearman",  start=st, count=ct))[k]
    nn[[m]] <- as.vector(ncvar_get(nc, "n_obs",             start=st, count=ct))[k]
  }
  nc_close(nc)
  for (m in 1:12) {
    if (is.null(rs[[m]])) next
    ok <- !is.na(rp[[m]]) & !is.na(rs[[m]]) & !is.na(nn[[m]]) & nn[[m]] >= NMIN
    if (sum(ok) < 50) next
    a <- rp[[m]][ok]; b <- rs[[m]][ok]
    qa <- p.adjust(pp[[m]][ok], "BH"); qb <- p.adjust(ps[[m]][ok], "BH")
    sa <- qa <= Q; sb <- qb <= Q
    po <- mean(sa == sb); f1 <- mean(sa); f2 <- mean(sb)
    pe <- f1*f2 + (1-f1)*(1-f2)
    res <- rbind(res, data.frame(parque=p, var=v, escala=e, mes=m, n=sum(ok),
      r_coef = cor(a, b), signo = mean(sign(a) == sign(b)),
      acuerdo = po, kappa = if (pe < 1) (po-pe)/(1-pe) else NA,
      sig_p = 100*f1, sig_s = 100*f2))
  }
  rm(rp, rs, pp, ps, nn); gc(verbose=FALSE)
}
saveRDS(res, file.path(OUT, "pvs.rds"))
cat(sprintf("combinaciones parque-variable-escala-mes evaluadas: %d\n", nrow(res)))
cat(sprintf("pixeles-mes en total: %s\n\n", format(sum(res$n), big.mark=".")))
m <- function(x) median(x, na.rm=TRUE)
cat("=== GLOBAL (mediana sobre todas las combinaciones) ===\n")
cat(sprintf("  correlacion entre los dos coeficientes : %.3f   (rango %.2f-%.2f)\n",
    m(res$r_coef), quantile(res$r_coef,.01,na.rm=TRUE), max(res$r_coef,na.rm=TRUE)))
cat(sprintf("  acuerdo en el SIGNO                    : %.1f %%\n", 100*m(res$signo)))
cat(sprintf("  acuerdo en la DECISION de significacion: %.1f %%  (kappa %.3f)\n",
    100*m(res$acuerdo), m(res$kappa)))
cat(sprintf("  superficie significativa: Pearson %.1f %%  Spearman %.1f %%\n",
    mean(res$sig_p), mean(res$sig_s)))
