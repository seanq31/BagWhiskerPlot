# BagWhiskerPlot

Implementation for the paper 'The Bag-and-Whisker Plot: A New Bagplot for Bivariate Data'.

This package exposes a single high-level function:

- `bag_whisker(x, y, ...)` — compute and (by default) plot a bag-and-whisker plot. `x` and `y` are numeric vectors or `x` may be a two-column matrix/data.frame.

## Usage example

```r
library(devtools)
install_github('seanq31/BagWhiskerPlot')

library(BagWhiskerPlot)

set.seed(1)
n <- 200
clean <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, 0.6, 0.6, 1), 2))

contam_frac1 <- 0.10
k1 <- floor(n * contam_frac1)
contam1 <- cbind(rnorm(k1, 6, 1), rnorm(k1, -6, 1))
dat_clean <- clean
dat_cont1 <- clean
dat_cont1[1:k1, ] <- contam1

bag_whisker(dat_clean, type1 = 'FDR', q = 0.1, main = "Clean")
bag_whisker(dat_cont1, type1 = 'FDR', q = 0.1, main = paste0("Contaminated ", contam_frac1 * 100, "%"))
```

Code for the figures in the paper can be found under the `vignette` folder.