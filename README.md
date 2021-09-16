# splattdr
The splattdr package is a wrapper for [`Splatter`][Splatter] allowing the simulation of dose-response single-cell RNA sequencing data.

## Installation
Splattdr depends on the [`Splatter`][Splatter] package (tested with version 1.14.1) which can be installed as follows
```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("splatter")
```

Then install splattdr directly from Github
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

## Citing splattdr
Please cite ["Nault, R., Saha, S., Bhattacharya, S., Dodson, J., Sinha, S., Maiti, T. and Zacharewski, T. (2021). Benchmarking of a Bayesian single cell RNAseq differential gene expression test for dose-response study designs. bioRxiv; doi.org/10.1101/2021.09.08.459475"][paper]

```
@article{,
   author = {Nault, Rance and Saha, Satabdi and Bhattacharya, Sudin and Dodson, Jack and Sinha, Samiran and Maiti, Tapabrata and Zacharewski, Tim},
   title = {Benchmarking of a Bayesian single cell RNAseq differential gene expression test for dose-response study designs},
   journal = {bioRxiv},
   pages = {2021.09.08.459475},
   DOI = {10.1101/2021.09.08.459475},
   url = {http://biorxiv.org/content/early/2021/09/10/2021.09.08.459475.abstract},
   year = {2021},
   type = {Journal Article}
}
```

[Splatter]: http://bioconductor.org/packages/release/bioc/html/splatter.html
[paper]: https://www.biorxiv.org/content/10.1101/2021.09.08.459475v1.full