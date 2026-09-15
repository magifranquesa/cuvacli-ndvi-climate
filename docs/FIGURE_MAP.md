# Figure → script map

Output folders are under `<ROOT>/articulo_2/figuras/`. "Manual" means the final layout was
composed by hand (ArcGIS Pro or image editor) from the listed outputs.

## Manuscript

| Figure | Content | Script | Output |
|---|---|---|---|
| 1 | Study area, park types | — (map, manual) | — |
| 2 | Monthly climate trends 1984–2023 by park | not in this repository (trend analysis: Theil–Sen + Mann–Kendall with `zyp`) | — |
| 3 | Monthly SPEI trends | idem | — |
| 4 | NDVI–Pr, positive and negative, by park, stacked by scale | `04_figures_park_level/figuras_correlaciones_superficie_ndvi_clima.R` | `figuras_mensuales_v2/fdr/cor_combined/` |
| 5 | NDVI–Tmax | idem | idem |
| 6 | NDVI–Tmin | idem | idem |
| 7 | NDVI–SPEI, positive | `04_figures_park_level/figuras_correlaciones_superficie_ndvi_spei.R` | `figuras_mensuales_v2/fdr/cor_spei_pos/` |
| 8 | Maps: strongest NDVI–SPEI correlation, month, scale (mainland + Balearics) | `07_maps/max_cor2tif.R` → ArcGIS Pro (manual) | `mapas_v2/fdr/rasters_max_cor_positive/spei/` |
| 9 | Positive correlations by vegetation type, network level | `05_vegetation_types/figura_snv_clima.R` | `SNV_v2/fdr/global_pos/snv_response_clima_ponderado.*` |
| 10 | Negative correlations, three park × type examples: (a) Cabrera coniferous, (b) Guadarrama alpine scrub/grassland, (c) Doñana salt marsh | panels from `05_vegetation_types/figura_snv_clima_por_parque_negativos.R` (S71, S74, S73), composed manually | `SNV_v2/fdr/por_parque_neg/` |

## Supplement

| Figures | Content | Script | Output |
|---|---|---|---|
| S1–S8 | Maps by variable and sign; S7–S8 Islas Atlánticas | `07_maps/max_cor2tif.R` → ArcGIS Pro (manual) | `mapas_v2/fdr/rasters_max_cor_*/` |
| S9–S19 | Positive correlations by vegetation type, per park (11 parks) | `05_vegetation_types/figura_snv_clima_por_parque.R` | `SNV_v2/fdr/por_parque_pos/` |
| S20–S23 | Positive, network level, stacked by scale (SPEI, Pr, Tmax, Tmin) | `06_figures_stacked_scales/snv_response_clima_apilado_absoluto_positive_global.R` | `figuras_mensuales_apilados_absolutos_SNV_v2/fdr/ndvi/global_pos/` |
| S24–S67 | Positive, per park, stacked by scale (11 parks × 4 variables; SPEI S24–S34, Pr S35–S45, Tmax S46–S56, Tmin S57–S67) | `06_figures_stacked_scales/snv_response_apilado_absoluto_positivo_por_parque.R` | `…/por_parque_pos/` |
| S68 | Negative correlations by vegetation type, network level | `05_vegetation_types/figura_snv_clima_negativos.R` | `SNV_v2/fdr/global_neg/` |
| S69–S79 | Negative, per park | `05_vegetation_types/figura_snv_clima_por_parque_negativos.R` | `SNV_v2/fdr/por_parque_neg/` |
| S80–S82 | Negative, network level, stacked by scale (Pr, Tmax, Tmin) | `06_figures_stacked_scales/snv_response_clima_apilado_absoluto_negative_global.R` | `…/global_neg/` |
| S83–S115 | Negative, per park, stacked by scale (Pr S83–S93, Tmax S94–S104, Tmin S105–S115) | `06_figures_stacked_scales/snv_response_apilado_absoluto_negativo_por_parque.R` | `…/por_parque_neg/` |

Park order within each block is alphabetical by folder name: aiguestortes, cabaneros, cabrera,
daimiel, donana, guadarrama, cies (Islas Atlánticas), monfrague, ordesa, picoseuropa,
sierranevada. Sierra de las Nieves has no SNV coverage and is absent from S9–S115.

## Complementary figures

| Figure | Content | Script |
|---|---|---|
| NDVI vs kNDVI profiles | Monthly profiles for the 12 park × variable combinations with ≥ 20 % NDVI area, with per-panel statistics | `08_supplementary_analyses/fig_xval_kndvi.R` |
| NDVI vs kNDVI scatter | Pixel-level scatter of the two coefficients, same 12 combinations | `08_supplementary_analyses/fig_xval_scatter.R` |

## Two denominators, stated in the captions

- **Figs. 9, S68** (network, by type): `sum(n_sig) / sum(n_total)` over the five scales, i.e.
  the fraction of pixel × scale combinations that are significant. Values are lower than in the
  maximum-based figures because a pixel rarely is significant at all five scales.
- **Figs. 10, S9–S19, S69–S79** (per park, by type): the scripts plot the per-scale fraction
  with `position_dodge()`, so the five scales overlap and the visible bar is the scale with the
  largest significant fraction for that month and variable.
- **Figs. 4–7, S20–S67, S80–S115** (maxima): each pixel counts once, at the scale of its
  strongest significant correlation; bars are % of park (or vegetation-type) area.
