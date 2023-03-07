# STUtility

## STUtility becomes `semla`

STUtility will no longer be maintained and instead you can use our new R package 
`semla` which is available at https://github.com/ludvigla/semla.

## Installation

To install STUtility from github (currently only option), you need to have devtools installed:

```
install.packages("remotes")
```

Finally, install STUtility:

```
remotes::install_github(
    "jbergenstrahle/STUtility"
)
```

### Install STUtility in an anaconda environment (Mac OS)

From the terminal

```
conda create -n R4.0
conda activate R4.0
conda install -c conda-forge r-essentials r-base r-devtools r-spdep r-hdf5r
conda install -c bioconda r-fftwtools
```

From R

```
devtools::install_github(
    "jbergenstrahle/STUtility"
)
```

## How to use the package

Please see  https://ludvigla.github.io/STUtility_web_site/ for online documentation and tutorials in how to use the package
