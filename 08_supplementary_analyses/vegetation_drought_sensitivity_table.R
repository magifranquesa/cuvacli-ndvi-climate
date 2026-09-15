# ==============================================================================
# SENSIBILIDAD A LA SEQUÍA POR TIPO DE VEGETACIÓN (red completa)
# ==============================================================================
# Para cada tipo de vegetación: mes y valor del máximo de superficie con
# correlación positiva significativa NDVI-SPEI, escala de acumulación dominante
# entre los píxeles significativos y proporción de respuestas a 9-12 meses.
#
# Entrada: resultados_v2/snv_respuesta_spei_ndvi_absoluto_positivo.csv
# ==============================================================================

if (!exists("ROOT")) ROOT <- Sys.getenv("CUVACLI_ROOT", "E:/IPE/PROYECTOS/CUVACLI")

suppressMessages({library(dplyr); library(readr)})
options(width = 140)
R2 <- paste0(ROOT, "/articulo_2/resultados_v2/")
d <- read_csv(paste0(R2,"snv_respuesta_spei_ndvi_absoluto_positivo.csv"), show_col_types = FALSE)
nom <- c(AridScr="Hyperxerophilous garrigues & scrub", MedScrub="Mediterranean scrublands",
         HaloVeg="Halophilous vegetation", MedSclF="Mediterranean sclerophyllous forest",
         MedConF="Mediterranean coniferous forest", Crop="Crops", SMedMarF="Semi-Med. marcescent forest",
         HydroRipVeg="Hydrophilous & riparian veg.", Grass="Grasslands", SaltMar="Salt marshes",
         AlpScrGr="Alpine scrub & grassland", Reforestation="Reforestation",
         AlpConF="Alpine coniferous forest", AtlScrub="Atlantic scrublands",
         TempDecF="Temperate broadleaf deciduous forest")

# % de superficie de cada tipo con correlacion significativa, por mes (red completa)
mes <- d %>% filter(Id_SNVeg %in% names(nom)) %>%
  group_by(Id_SNVeg, mes) %>%
  summarise(sig = sum(n), tot = sum(n_total[escala == 1]), .groups = "drop") %>%
  mutate(p = 100 * sig / tot)
pico <- mes %>% group_by(Id_SNVeg) %>% slice_max(p, n = 1, with_ties = FALSE) %>%
  select(Id_SNVeg, mes_pico = mes, p_pico = p)

# distribucion de escalas entre los pixeles significativos
esc <- d %>% filter(Id_SNVeg %in% names(nom)) %>%
  group_by(Id_SNVeg, escala) %>% summarise(n = sum(n), .groups = "drop") %>%
  group_by(Id_SNVeg) %>% mutate(pc = 100 * n / sum(n))
dom <- esc %>% slice_max(pc, n = 1, with_ties = FALSE) %>%
  select(Id_SNVeg, esc_dom = escala, esc_pc = pc)
largo <- esc %>% filter(escala >= 9) %>% group_by(Id_SNVeg) %>%
  summarise(pc_912 = sum(pc), .groups = "drop")

r <- pico %>% left_join(dom, "Id_SNVeg") %>% left_join(largo, "Id_SNVeg") %>% arrange(-p_pico)
cat("SENSIBILIDAD A LA SEQUIA POR TIPO DE VEGETACION (red completa, NDVI-SPEI positivo, v2+FDR)\n\n")
cat(sprintf("%-38s %9s %6s %9s %9s\n","tipo de vegetacion","pico %","mes","escala dom","9-12 m %"))
for (i in 1:nrow(r)) cat(sprintf("%-38s %8.1f%% %6s %6d m %7.0f%% %8.0f%%\n",
    nom[[r$Id_SNVeg[i]]], r$p_pico[i], month.abb[r$mes_pico[i]], r$esc_dom[i], r$esc_pc[i], r$pc_912[i]))
