# Workflow: execution order, inputs, outputs, timings

All paths are relative to `ROOT` (see `config.R`). Times are for a desktop workstation
(R 4.4, single thread) and are indicative. Scripts within a phase are independent of each
other unless stated. Parameters common to the whole chain:

| parameter | value | where |
|---|---|---|
| `indice` | `"ndvi"` (`"kndvi"` for the cross-validation) | phases 1–7 |
| `N_MIN` | 20 valid years per pixel | phases 3, 5 |
| `Q_FDR` | 0.05 (Benjamini–Hochberg) | phases 3, 5 |

---

## Phase 1 — correlations per accumulation scale  `01_correlations/`  (needs the climate grid)

| script | what it does | output |
|---|---|---|
| `correlation_vi_nc.r` | For each park, index (NDVI, kNDVI), variable (Pr, Tmax, Tmin) and scale (1, 3, 6, 9, 12 months): detrends both series per pixel and calendar month (OLS against calendar year, after masking missing years), Pearson r, two-sided p (df = n − 3), number of valid years; Spearman ρ on a reproducible 5 % pixel subsample (`set.seed(42)`). | `articulo_2/correlaciones_v2/<park>/correlation_<index>_<var>_scale_<s>_<park>.nc` |
| `correlation_spei_vi.r` | Same for SPEI. | idem, `var = spei` |

480 files (12 parks × 2 indices × 4 variables × 5 scales). ~1 day in total; Sierra Nevada is the slowest park.
Climate cells missing in Islas Atlánticas, Cabrera, Ordesa and Doñana are filled by an iterative 3 × 3 focal mean before bilinear resampling to the 30 m grid (`parques_gapfill`).

## Phase 2 — strongest response across scales  `02_maxima/`

| script | output |
|---|---|
| `cor_maximas_pos_variables_seleccion_escalas.R` | `correlaciones_v2/maximos_cor_positive/<park>_max_pos_correlation_<var>_<index>_selscales.nc` |
| `cor_maximas_pos_spei.R` | `…/maximos_cor_positive/<park>_max_pos_correlation_spei_<index>.nc` |
| `cor_maximas_neg_variables_seleccion_escalas.R` | `…/maximos_cor_negative/…_max_neg_…_selscales.nc` |
| `cor_maximas_neg_spei.R` | `…/maximos_cor_negative/…_max_neg_correlation_spei_<index>.nc` |

Variables: `max_correlation`, `scale`, `p_value` (uncorrected, at the selected scale), `significance` (raw code, p ≤ 0.05), `n_obs`. 96 + 96 files. ~2 h.
Note the file-name asymmetry: `_selscales` for Pr/Tmax/Tmin, none for SPEI; the SPEI chain is kept separate throughout.

## Phase 3 — minimum years and FDR  `03_fdr/`

| script | output |
|---|---|
| `fdr_postproceso_pos.R` | `correlaciones_v2/maximos_cor_positive_fdr/` (same file names) |
| `fdr_postproceso_neg.R` | `correlaciones_v2/maximos_cor_negative_fdr/` |

Sets `significance` to 2 / −2 only where `n_obs ≥ 20` and the BH-adjusted p-value within the park–month field is ≤ 0.05; NA where the pixel is below the years threshold. Minutes.

**These `_fdr` folders are what every figure reads.**

## Phase 4 — park-level figures  `04_figures_park_level/`  (Figs. 4–7)

| script | figures | output |
|---|---|---|
| `figuras_correlaciones_superficie_ndvi_clima.R` | Figs. 4, 5, 6 (Pr, Tmax, Tmin; positive and negative) | `articulo_2/figuras/figuras_mensuales_v2/fdr/cor_combined/` |
| `figuras_correlaciones_superficie_ndvi_spei.R` | Fig. 7 (SPEI, positive) | `…/figuras_mensuales_v2/fdr/cor_spei_pos/` |

Bars: % of park area (pixels with n ≥ 20 in at least one month) significant in that month, stacked by the scale of the maximum. Minutes.

## Phase 5 — vegetation types  `05_vegetation_types/`  (Figs. 9–10, S9–S19, S68–S79)

This chain does **not** use the maxima: it reads the per-scale files of phase 1 and applies n ≥ 20 and FDR itself, so that a pixel can count at more than one scale.

| step | script | notes |
|---|---|---|
| 5.1 | `snv_respuesta_clima_mensual_pos.R` | positive, Pr/Tmax/Tmin — ~1 h |
| 5.2 | `snv_respuesta_clima_spei_mensual.R` | positive, SPEI |
| 5.3 | `snv_respuesta_clima_mensual_neg.R` | negative, Pr/Tmax/Tmin |
| 5.4 | `combinar_csvs.R` | needs 5.1 + 5.2 → `resultados_v2/snv_respuesta_clima_mensual_positivo_ndvi_completo.csv` |
| 5.5 | `combinar_csvs_neg.R` | needs 5.3 → `…_negativo_ndvi_completo.csv` |
| 5.6 | `figura_snv_clima.R` | **Fig. 9** (network, positive; denominator = pixels × 5 scales) |
| 5.7 | `figura_snv_clima_negativos.R` | S68 (network, negative) |
| 5.8 | `figura_snv_clima_por_parque.R` | S9–S19 (per park, positive; bar = scale with the largest significant fraction) |
| 5.9 | `figura_snv_clima_por_parque_negativos.R` | S69–S79 (per park, negative); **Fig. 10** is composed from three of these panels |

Outputs in `resultados_v2/` (copied to `results/` here) and `articulo_2/figuras/SNV_v2/fdr/{global_pos,global_neg,por_parque_pos,por_parque_neg}/`.
Sierra de las Nieves has no SNV shapefile and is absent from this phase.

## Phase 6 — vegetation types stacked by scale  `06_figures_stacked_scales/`  (S20–S67, S80–S115)

Reads the FDR maxima of phase 3 (each pixel counts once, at the scale of its maximum).

| step | script | figures |
|---|---|---|
| 6.1–6.3 | `generar_csv_pos_clima.R`, `generar_csv_pos_spei.R`, `generar_csv_neg_clima.R` | → `resultados_v2/snv_respuesta_<var>_ndvi_absoluto_<sign>.csv` (~30 min) |
| 6.4 | `snv_response_clima_apilado_absoluto_positive_global.R` | S20–S23 (needs 6.1–6.2) |
| 6.5 | `snv_response_clima_apilado_absoluto_negative_global.R` | S80–S82 (needs 6.3) |
| 6.6 | `snv_response_apilado_absoluto_positivo_por_parque.R` | S24–S67 (reads the maxima directly; ~20 min) |
| 6.7 | `snv_response_apilado_absoluto_negativo_por_parque.R` | S83–S115 (idem; the SPEI figures it also produces are not in the supplement) |

Output: `articulo_2/figuras/figuras_mensuales_apilados_absolutos_SNV_v2/fdr/ndvi/`.

## Phase 7 — map layers  `07_maps/`  (Fig. 8, S1–S8)

`max_cor2tif.R` exports, for every park, variable and sign, three GeoTIFFs (strongest correlation, its month, its scale) from the FDR maxima:
`articulo_2/figuras/mapas_v2/fdr/rasters_max_cor_{positive,negative}/<var>/<park>_<var>_<pos|neg>_<index>_{cor,mes,esc}_max.tif`.
The map layouts themselves are composed manually in ArcGIS Pro from these layers.

## Complementary analyses  `08_supplementary_analyses/`

| script | reads | prints |
|---|---|---|
| `kndvi_crossvalidation.R` + `kndvi_crossvalidation_summary.R` | FDR maxima, NDVI and kNDVI | profile r, pixel r, κ, scale agreement for 96 park × variable × sign combinations; the 12 with ≥ 20 % area |
| `fig_xval_kndvi.R` | idem | figure: NDVI vs kNDVI monthly profiles for those 12 combinations |
| `fig_xval_scatter.R` | idem | figure: pixel-level scatter of the two coefficients |
| `spearman_vs_pearson.R` | per-scale files (5 % subsample) | r between coefficients, sign and decision agreement, significant area |
| `gapfill_sensitivity.R` | FDR maxima + climate grid | significant area with all pixels vs native-climate pixels, Cabrera and Cíes |
| `spei_trends_summary.R` | `resultados/clima/tendencias_mensuales_spei_por_parque_1984_2023.csv` | significant SPEI trends per park |
| `vegetation_drought_sensitivity_table.R` | `results/snv_respuesta_spei_ndvi_absoluto_positivo.csv` | peak month, dominant scale by vegetation type |
| `tmax_seasonal_reversal.R` | FDR maxima | winter-positive vs summer-negative Tmax by park |

Small outputs (`.rds`) go to `articulo_2/resultados_v2/analisis_complementarios/`; figures to `articulo_2/figuras/complementarias/`.
