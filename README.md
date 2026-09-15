# Vegetation sensitivity to interannual climate variability in Spain's National Parks

Code for the analyses and figures of

> Franquesa, M., Adell-Michavila, M., Vicente-Serrano, S.M. *Vegetation sensitivity to interannual climate variability in Spain's National Parks: seasonal patterns and accumulation timescales.* International Journal of Geoheritage and Parks (in revision).

Starting from the monthly NDVI and kNDVI series of Franquesa et al. (2025), the workflow computes pixel-wise correlations between vegetation activity and climate at 30 m for the twelve national parks of mainland Spain and the Balearic and Atlantic islands over 1984–2023; selects the strongest response across accumulation timescales; controls the false discovery rate; aggregates the results by park and by vegetation type; and produces every figure of the manuscript and supplement, plus a set of complementary analyses.

Data are on Zenodo: [doi:10.5281/zenodo.18469776](https://doi.org/10.5281/zenodo.18469776) (version 2.0 matches this code).

## Layout

| Folder | Content | Manuscript |
|---|---|---|
| `01_correlations/` | Detrended interannual Pearson correlation between NDVI (or kNDVI) and each climate variable accumulated over 1, 3, 6, 9 and 12 months, for each pixel and calendar month; Spearman on a 5 % pixel subsample. | Sect. 2.3.2 |
| `02_maxima/` | For each pixel and month, the strongest positive (or negative) correlation across the five scales, with the scale at which it occurs. | Sect. 2.3.2 |
| `03_fdr/` | Minimum-years threshold (n ≥ 20) and Benjamini–Hochberg FDR control (q = 0.05) applied independently to each park–month field. | Sect. 2.3.2 |
| `04_figures_park_level/` | Monthly extent of significant correlations by park, stacked by accumulation scale. | Figs. 4–7 |
| `05_vegetation_types/` | Aggregation by Natural Vegetation Systems (SNV) at network and park level, positive and negative correlations. | Figs. 9–10, S9–S19, S68–S79 |
| `06_figures_stacked_scales/` | Vegetation-type figures stacked by accumulation scale. | S20–S67, S80–S115 |
| `07_maps/` | GeoTIFF export (correlation, month, scale) used to compose the maps in ArcGIS Pro. | Fig. 8, S1–S8 |
| `08_supplementary_analyses/` | Complementary analyses: NDVI–kNDVI cross-validation, Pearson vs Spearman, sensitivity to climate gap-filling, SPEI trends, drought sensitivity by vegetation type, seasonal reversal of the Tmax response. | Sects. 2.2, 2.3, 4.4 |
| `results/` | CSV tables behind the vegetation-type figures (network and park level). | Figs. 9–10, S9–S115 |
| `docs/` | Execution order, inputs, outputs and timings; figure-to-script map. | |

`docs/WORKFLOW.md` gives the execution order, inputs, outputs and approximate run times of every script. `docs/FIGURE_MAP.md` says which script produces each figure.

## Running the code

**Language.** R (4.4) with `terra`, `ncdf4`, `dplyr`, `tidyr`, `purrr`, `readr`, `ggplot2`, `patchwork`, `svglite`, `abind` and `sf`. Comments inside the scripts are in Spanish; every folder is described in English here and in `docs/`.

**Paths.** All paths hang from a single root. Set the environment variable `CUVACLI_ROOT` (or edit `config.R`) to the folder that holds the data in the layout described in `config.R`. Nothing else needs editing.

**What can be run from the public data.** Phases 1–3 (`01_`–`03_`) need the 1.1 km gridded climate dataset, which cannot be redistributed; their outputs — the per-scale correlation files and the FDR-screened maxima — are on Zenodo. Everything from `04_` onwards, and most of `08_`, runs from the Zenodo files alone:

```
Zenodo <park>.zip  →  <ROOT>/articulo_2/correlaciones_v2/<park>/                       (correlations/)
                      <ROOT>/articulo_2/correlaciones_v2/maximos_cor_positive_fdr/     (maximum_positive_correlations/)
                      <ROOT>/articulo_2/correlaciones_v2/maximos_cor_negative_fdr/     (maximum_negative_correlations/)
                      <ROOT>/MOSAICOS/mosaicos_all/<park>_NDVI_filled_clean.nc         (ndvi/<park>_NDVI.nc, renamed)
```

Inside the Zenodo zips the maxima files for Pr, Tmax and Tmin are named without the `_selscales` suffix that the scripts expect; either rename them or use the two-candidate lookup already present in `08_supplementary_analyses/*.R`.

**NDVI series.** The 30 m monthly NDVI and kNDVI series (Landsat 5/7/8, 1984–2023) were generated for Franquesa et al. (2025) and are not produced by this repository; they are included in the Zenodo dataset. kNDVI uses a fixed σ = 0.15, so it is not a monotonic transform of NDVI.

**Other inputs.** Park boundaries (OAPN, `Limites_PN.shp`) and the Natural Vegetation Systems maps (`SNV_<park>.shp`, MITECO) are public and are downloaded separately; the SPEI grids come from the Spanish Drought Monitor (https://monitordesequia.csic.es). Sierra de las Nieves has no SNV coverage and therefore does not appear in the vegetation-type figures.

## Method in brief

- Both series are detrended by ordinary least squares against calendar year after masking missing values; significance uses n − 3 degrees of freedom (`01_correlations`).
- Pixels with fewer than 20 valid years are discarded (`03_fdr`).
- The false discovery rate is controlled per park–month field (Benjamini & Hochberg, 1995; Wilks, 2016); no correction is applied for the selection of the best scale, since scales are nested and the maximum is reported as a descriptive quantity (`03_fdr`).
- Significance masks are NA-safe throughout: where `significance` is `NA` (pixel below the years threshold) the pixel counts neither as significant nor in the denominator.
- The vegetation-type figures use two different denominators, stated in each caption: Figs. 9 and S68 divide by pixels × scales; the park-level Figs. 10 and S9–S19/S69–S79 show, for each month and variable, the scale with the largest significant fraction; the maxima-based figures (4–7, S20–S67, S80–S115) count each pixel once.

## Citation

If you use this code, please cite the paper above and the dataset (doi:10.5281/zenodo.18469776). A `CITATION.cff` file is included.

## License

Code: MIT (see `LICENSE`). Data on Zenodo: CC BY 4.0.
