GenomicTrackRepository
==============================

Store scripts used to make Bushman Lab genomic tracks   

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data/figures that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── environment        <- Conda environment installation requirements.
    │
    ├── exploratory        <- Exploratory analyses that may or may not be formalized into a notebook.
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Rmarkdown, markdown, and Jupyter notebooks. Naming convention is a number (for ordering) 
    │                         (Ordering is absolute and indicates order each notebook could be run to replicate analysis)
    │                         and `-` delimited description, e.g.
    │                         `1.0_initial_data_exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analyses as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Finished graphics and figures to be used in reports.
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── Rlib           <- Default Helper Rlibraries
    │   │   └── utils.R
    │   │
    │   ├── ProjectPaths   <- Project path management for python.
    │   │   └── ProjectPaths.py
    │
    ├── tests              <- Tests for more formal src packages.
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------

# General idea

notebooks: formal analyses. These can be run to replicate all data and analyses generated. Numbers indicate the order each notebook can be run. For example 0-my_analysis must be run before 1-next_analysis because the numbering indicates that 1-next_analysis may use data generated be 0-my_analysis. 

src: libraries and scripts to generate data.

# Organizing analyses

All Folders are subdivied by the notebook name.

**Notebook analysis**

`notebook/0-first_analysis.rmd`

**Interim data**

`data/intermin/0-first_analysis`

**models**

`models/0-first_analysis`

**Helper R libraries(s)**

`Rlib/0_first_analysis.R `

**Helper python libraries (x at start because filename can't start with a number)**

`lib/x0_first_analysis.py`











# GenomicTrackRepository
