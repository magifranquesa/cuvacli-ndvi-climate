# ==============================================================================
# INVERSIÓN ESTACIONAL DE LA RESPUESTA A Tmax POR PARQUE
# ==============================================================================
# Para cada parque, máximo mensual de superficie con correlación positiva con
# Tmax en invierno (dic-mar) frente a máximo de correlación negativa en verano
# (jun-ago), y el menor de los dos como medida de la coexistencia de ambos
# comportamientos.
#
# Entrada: correlaciones_v2/maximos_cor_{positive,negative}_fdr/ (Tmax, NDVI)
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

suppressMessages(library(ncdf4)); options(width=130)
B <- paste0(ROOT, "/articulo_2/correlaciones_v2/")
PQ <- c("aiguestortes","cabaneros","cabrera","cies","daimiel","donana",
        "guadarrama","monfrague","ordesa","picoseuropa","sierranevada","sierranieves")
perf <- function(sg,p){
  d <- sprintf("%smaximos_cor_%s_fdr/", B, if(sg=="pos")"positive" else "negative")
  s0 <- if(sg=="pos")"max_pos_correlation" else "max_neg_correlation"
  nc <- nc_open(sprintf("%s%s_%s_tmax_ndvi_selscales.nc", d, p, s0))
  s <- ncvar_get(nc,"significance"); nc_close(nc)
  den <- sum(apply(!is.na(s), c(1,2), any)); cod <- if(sg=="pos")2 else -2
  r <- sapply(1:12, function(m) 100*sum(s[,,m]==cod, na.rm=TRUE)/den); rm(s); gc(verbose=FALSE); r }
cat("INVERSION ESTACIONAL DE Tmax: positiva en invierno, negativa en verano\n\n")
cat(sprintf("%-14s %22s %22s %9s\n","parque","max invierno (dic-mar)","max verano (jun-ago)","el menor"))
r <- NULL
for (p in PQ) {
  a <- perf("pos",p); b <- perf("neg",p)
  inv <- c(12,1,2,3); ver <- 6:8
  mi <- inv[which.max(a[inv])]; mv <- ver[which.max(b[ver])]
  r <- rbind(r, data.frame(p, w=a[mi], wm=month.abb[mi], s=b[mv], sm=month.abb[mv], min=min(a[mi],b[mv])))
}
for (i in order(-r$min)) cat(sprintf("%-14s %14.1f%% (%s) %14.1f%% (%s) %8.1f%%\n",
    r$p[i], r$w[i], r$wm[i], r$s[i], r$sm[i], r$min[i]))
