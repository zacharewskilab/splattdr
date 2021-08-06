# splattDR
This package implements a new Bayesian test for detecting differential gene expression over multiple dose groups in single cell gene expression studies. 

## Installation
Be sure to install dependencies before installing splattDR
```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("splatter")
```

Then install splattDR directly from Github
```{r}
library("devtools")
devtools::install_github("zacharewskilab/splattdr")
```

## Simulating dose-response data
In order to develop new testing methods for differential expression analysis we need data for which the ground truth is known.
Here we simulate dose-response data based on known dose-response models and using parameters derived from real data.

_dose-response models derived from [BMD Express](https://bmdexpress-2.readthedocs.io/en/feature-readthedocs/) and the [US EPA BMDS Documentation](https://www.epa.gov/bmds/benchmark-dose-software-bmds-version-27-user-manual)_

* Hill
* Exponential 2 - 5
* Power
* Linear
* Polynomial 2-4

