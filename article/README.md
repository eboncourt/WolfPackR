# Scripts of the article

R scripts reproducing the case study of Boncourt, *WolfPackR: An R package for identifying wolf packs based on genetic and spatial data* (Ecological Solutions and Evidence, in revision), with WolfPackR 1.1.0 and the Carpathian data shipped with the package. Run them from the root of the repository; outputs go to `article/output/`.

- `figures_article.R`: all analyses, Table 2, Tables S1–S2, Figs 1–3, and the numbers cited in the text (`numbers_in_text.txt`); about three minutes.
- `figure1_workflow.R`: Fig. 1 (workflow) on its own, handy to adjust its layout.
- `comparison_MCP_v094.R`: the earlier approach based on minimum convex polygons, run with WolfPackR 0.9.4 installed in a separate library.

Needs `sf`, `ggplot2` and `patchwork` (plus `remotes` and `callr` for the MCP comparison).
