# TANDEM
A two-stage regression method that can be used when various input data types are correlated, for example gene expression and methylation in drug response prediction. In the first stage it uses the upstream features (such as methylation) to predict the response variable (such as drug response), and in the second stage it uses the downstream features (such as gene expression) to predict the residuals of the first stage. In our manuscript [(Aben et al., 2016)](https://doi.org/10.1093/bioinformatics/btw449), we show that using TANDEM prevents the model from being dominated by gene expression and that the features selected by TANDEM are more interpretable.

The R package is available on [CRAN](https://cran.r-project.org/package=TANDEM/index.html). It can be installed within R using:
```r
install.packages("TANDEM")
```
A quick start can be found [here](https://cran.r-project.org/package=TANDEM/vignettes/my-vignette.html). The reference manual can be found [here](https://cran.r-project.org/package=TANDEM/TANDEM.pdf).
