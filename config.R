# ============================================================================
# config.R — single place for paths and analysis parameters
#
# Every script in this repository starts with
#     if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")
# so you can either (a) set the environment variable CUVACLI_ROOT, (b) source this
# file before running a script interactively, or (c) edit the default below.
#
# Expected layout under ROOT (see README.md, "Data"):
#   MOSAICOS/mosaicos_all/<park>_NDVI_filled_clean.nc       monthly NDVI / kNDVI, 30 m (Zenodo)
#   CLIMA/agregados_lag/<var>_scale_<s>.nc                   1.1 km climate, accumulated (not public)
#   CLIMA/indices/SPEI_monthly/spei<s>.nc                    SPEI (Spanish Drought Monitor)
#   limites_red/desc_Red_PN_LIM_Enero2023/Limites_PN.shp     park boundaries (OAPN)
#   articulo_2/correlaciones_v2/                             outputs of phases 1–3 (Zenodo)
#   articulo_2/resultados_v2/                                CSV tables (copy in results/)
#   articulo_2/figuras/*_v2/                                 figures
# ============================================================================
ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

# Analysis parameters used throughout (the scripts define them locally with the
# same values; they are repeated here for reference)
INDICE     <- "ndvi"   # "ndvi" for the manuscript; "kndvi" for the cross-validation
N_MIN      <- 20       # minimum valid years per pixel
Q_FDR      <- 0.05     # Benjamini–Hochberg false-discovery rate, per park–month field
ESCALAS    <- c(1, 3, 6, 9, 12)   # accumulation scales, months
ANIOS      <- 1984:2023
