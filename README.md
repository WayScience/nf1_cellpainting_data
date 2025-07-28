[![DOI](https://zenodo.org/badge/DOI/latest.svg)](https://doi.org/10.5281/zenodo.16539382)

# NF1 Schwann cell Cell Painting data

In this repository, we generate image analysis and image-based profiling pipelines to extract and format single-cell morphological profiles.

We train a machine learning model to predict NF1 genotype, evaluate, and generate figures within a separate repository called: [nf1_schwann_cell_morphology_signature](https://github.com/WayScience/nf1_schwann_cell_morphology_signature).
Please visit the above repository for further information on the generation, validation, and figures from this model for the manuscript.

## Goal

It is important to study Schwann cells from NF1 patients because NF1 causes patients to develop neurofibromas, which are peripheral nerve tumors forming bumps underneath the skin that appear due to the decrease of Ras-GAP neurofibromin production.
This decrease in production occurs when the NF1 gene is mutated (NF1 +/-).

**The goal of this project is to predict NF1 genotype from Schwann cell morphology.**
We apply single cell image analysis to Cell Painting images and use representation learning to extract morphology features.
We will apply machine learning to the morphology features to discover a biomarker of NF1 genotype.
Once we discover a biomarker from these cells, we hope that our method can be used for drug discovery to treat this rare disease.

## Data

The data we use is a modified Cell Painting assay on [Schwann cells](https://www.ncbi.nlm.nih.gov/books/NBK544316/) from patients with [Neurofibromatosis type 1 (NF1)](https://medlineplus.gov/genetics/condition/neurofibromatosis-type-1/).
The images are publicly available on figshare, under the [NF1 Schwann Cell Genotype Cell Painting Assay project](https://figshare.com/projects/NF1_Schwann_Cell_Genotype_Cell_Painting_Assay/161620).

The data is as follows:

| Plate                | DOI                                          | Description                                                                                                                                                                                           |
| -------------------- | -------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Plate 1              | <https://doi.org/10.6084/m9.figshare.22233292> | Preliminary plate of 8 wells with image sets of three Cell Painting channels for wildtype and null cells.                                                                                             |
| Plate 2              | <https://doi.org/10.6084/m9.figshare.22233700> | Preliminary plate of 32 wells with image sets of three Cell Painting channels for wildtype and null cells.                                                                                            |
| Plates 3 and 3 prime | <https://doi.org/10.6084/m9.figshare.22592890> | Plates utilized for modelling. Each contain 48 wells, with Plate 3 treated with 10% FBS and prime treated with 5% FBS. These plate contain all three *NF1* genotypes, with varying seeding densities. |
| Plate 4              | <https://doi.org/10.6084/m9.figshare.23671056> | Plate containing 60 wells with null and wildtype cells either not treated or treated with siRNAs. We do not include this plate for modelling or evaluation. The seeding density is 1000 cells.        |
| Plate 5              | <https://doi.org/10.6084/m9.figshare.26759914> | Plate containing 48 wells with all three *NF1* genotypes used for modelling. The seeding density is 1000 cells.                                                                                       |
| Plate 6              | <https://doi.org/10.6084/m9.figshare.28797659>                                         | Plate containing 60 wells with all three *NF1* genotypes across two cell lines (one from `iNFixion` and one from `MGH`) used to assess generalizability. The seeding density is 1000 cells.           |

There are two versions of the Cell Painting assay in this repository:

In this modified Cell Painting, there are three channels for plates 1 and 2:

- `DAPI` (Nuclei)
- `GFP` (Endoplasmic Reticulum)
- `RFP` (Actin)

![Modified_Cell_Painting.png](example_figures/Modified_Cell_Painting.png)

In this modified Cell Painting, there are four channels for all the rest of the plates:

- `DAPI` (Nuclei)
- `GFP` (Endoplasmic Reticulum)
- `CY5` (Mitochondria)
- `RFP` (Actin)

![Modified_CellPainting_Plate3.png](example_figures/Modified_CellPainting_Plate3.png)

For more information on plate maps and plate map figures, please go to the [metadata folder](./0.download_data/metadata/) in the first module.
All larger files, including `SQLite` outputs from CellProfiler and `parquet` processed data file from pycytominer, will need to be downloaded using git LFS after the repo is cloned.

## Repository Structure

| Module                                                | Purpose                                           | Description                                                                                                                                                    |
| :---------------------------------------------------- | :------------------------------------------------ | :------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [0.download_data](./0.download_data/)                 | Download NF1 data                                 | We download images from each plate of the NF1 dataset for analysis from Figshare                                                                               |
| [1.cellprofiler_ic](./1.cellprofiler_ic/)             | Apply CellProfiler illumination correction (IC)   | We use a CellProfiler pipeline to calculate and apply IC the images and save them                                                                              |
| [2.cellprofiler_analysis](./2.cellprofiler_analysis/) | Perform CellProfiler analysis on corrected images | We use a CellProfiler pipeline to segment single cells and extract features into a SQLite file                                                                 |
| [3.processing_features](./3.processing_features/)     | Process CellProfiler SQLite files                 | We use CytoTable to convert extracted features from SQLite files to parquet files. We then use pycytominer to annotate, normalize, and feature select profiles |
| [4.analyze_data](./4.analyze_data/)                   | Perform various analysis of morphology data       | Using different statistical methods, like linear modeling, we analyze the data to assess the difference in morphology between genotypes                        |

## Main environment

For all modules, we use one main environment for the repository, which includes all packages needed including installing CellProfiler v4.2.4.

To create the environment, run the below code block:

```bash
# Run this command in terminal to create the conda environment
conda env create -f nf1_cellpainting_env.yml
```

**Make sure that the conda environment is activated before running notebooks or scripts:**

```bash
conda activate nf1_cellpainting_data
```

## Licensing

- Code: BSD 3-Clause License (see [LICENSE](./LICENSE))
- Data: CC BY 4.0 License (see [LICENSE_DATA](./LICENSE_DATA))
