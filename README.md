# npdro
# Vote for n-Pdro
## It's Dynamite

Nearest-neighbor Projected Distance Regression - Optimized (<strong>npdro</strong>). npdro extends ndpr for optimizing neighborhoods and the number of neighbors to detect main effects and interactions in feature selection. Neighbor optimization methods include consensus nested cross-validation (cncv) with classification, varible-wise optimized-k (vwok), and principal component optimized-k (kPCA). Also includes optional speed optimizations.

#### Websites

[NPDR Github Page](https://insilico.github.io/npdro/)

[insilico Github Organization](https://github.com/insilico)

[insilico McKinney Lab](http://insilico.utulsa.edu/)

#### Related References. 

[2020 NDPR paper in Bioinformatics](https://doi.org/10.1093/bioinformatics/btaa024)

[2013 Gene-Wise Adaptive-Neighbors paper in PLoS One](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0081527)

### To install:

    >library(devtools)
    
    >install_github("insilico/npdro")  
    >library(npdro)

### Dependencies
To install the `tidyverse` collection of R packages:

```
install.packages('tidyverse')
```
To set `fast.reg = TRUE` or `fast.dist = TRUE` or `use.glmnet = TRUE`, please install the `speedglm` and `glmnet` packages:

```
install.packages(c('speedglm', 'wordspace', 'glmnet'))
```

If an issue arises with updating `openssl`, try updating it on your own system, e.g. for MacOS:

```brew install openssl@1.1```

### Examples

### Abstract

#### Contact
[brett.mckinney@gmail.com](brett.mckinney@gmail.com)
