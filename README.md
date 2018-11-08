cellAggregator <img src="man/figures/hex.png" align="right"  height="250" width="250"/>
  ======================================================

  CellAggregator

Overview
--------

  **cellAggregator** is a network method that simulates cell-cell aggregation assays
  
--------
  
<img src="man/figures/animation.gif" align="right"  height="250" width="1000"/>

--------

Installation
--------

```r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("shazanfar/cellAggregator")
library(cellAggregator)
```

Usage
-----

  The following example simulates the cell-cell aggregation of three populations expressing the same combinations of three separate proteins

```r
numCells = c(25,25,25)
numProtsPerCell = rbind(c(5,5,5),c(5,5,5), c(5,5,5))
bindingAffinity = cbind(c(1,0,0), c(0,1,0), c(0,0,1))

cellAggregationResult = cellAggregator(
  numCells,
  numProtsPerCell,
  bindingAffinity,
  timesteps = 100,
  burnIn = 0.75,
  verbose = TRUE,
  includeListOfGraphs = TRUE,
  plot = TRUE
)

cellAggregationResult$t_index
unlist(cellAggregationResult$listoftIndex)
cellAggregationBarplot(cellAggregationResult)


# another scenario, 2 populations and 3 proteins
numCells = c(25,25)
numProtsPerCell = rbind(c(5,0,0),c(5,5,5))
bindingAffinity = cbind(c(1,0,0), c(0,1,0), c(0,0,1))

cellAggregationResult = cellAggregator(
  numCells,
  numProtsPerCell,
  bindingAffinity,
  timesteps = 100,
  burnIn = 0.75,
  verbose = TRUE,
  includeListOfGraphs = TRUE,
  plot = TRUE
)

cellAggregationResult$t_index
unlist(cellAggregationResult$listoftIndex)
cellAggregationBarplot(cellAggregationResult)
```

## Author

* **Shila Ghazanfar**  - [@shazanfar](https://twitter.com/shazanfar)

