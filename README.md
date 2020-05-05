# STUtility

## Installation

STUtility requires R version >=3.6 and some packages within the Bioconductor suit needs to be installed prior to installing STUtility from github

If you don’t have bioconductor installed:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.11") #OBS! Require R >4.0.0
```

The following packages needs to be manually installed via BiocManager:

```
BiocManager::install("SingleCellExperiment")

BiocManager::install("multtest")

BiocManager::install("spdep")
```

To install STUtility from github (currently only option), you need to have devtools installed:

```
install.packages("devtools")
```

Additionaly, the package `NNLM` currently needs to be installed via github:

```
install_github('linxihui/NNLM')
```
Finally, install STUtility:

```
devtools::install_github(
    "jbergenstrahle/STUtility"
)
```

## How to use the package

Please see  https://ludvigla.github.io/STUtility_web_site/ for online documentation and tutorials in how to use the package
