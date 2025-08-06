# ggalign: Bridging the Grammar of Graphics and Biological Multilayered Complexity

This repository provides all the code and data used to generate the figures in
the manuscript` ggalign: Bridging the Grammar of Graphics and Biological Multilayered Complexity`. Each figure is supported by a dedicated Quarto (.qmd)
file and an associated data folder to ensure full transparency and
reproducibility.


## Script Structure
All scripts have been organized and renamed to match their corresponding figure
numbers in the manuscript. Each `.qmd` file contains the full analysis pipeline
used to generate that specific figure, including preprocessing, visualization,
and optional statistical analysis.

For example:

`Figure 3.qmd` contains the script used to generate `Figure 3`, along with the
necessary data processing steps.

`Figure 6.qmd` contains the script used to generate `Figure 6`

These files are written in Quarto format to allow seamless rendering into HTML
or PDF, and can also be opened as plain R scripts for interactive use.


## Data Organization
Each figure has a corresponding folder (e.g., `Figure3/`, `Figure6/`) that contains
the raw and processed data used in the analysis:

Raw datasets (when permissible by license and size constraints) are included
directly.

For large or publicly hosted datasets (e.g., from GEO or SRA), download
instructions are provided inside the respective `.qmd` file.

Intermediate and processed results are saved to support step-wise execution and
partial reruns.

>Note: Some files may be too large to upload directly to the repository. In such cases, clear instructions and download links are provided to retrieve the required data manually.
