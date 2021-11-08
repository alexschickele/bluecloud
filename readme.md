# Projecting functional plankton gene clusters

A. Schickele, P. Debeljack, L. Guidi, S. Ayata, J.O. Irisson

Sorbonne Université, CNRS, Laboratoire d'Océanographie de Villefranche (LOV), Villefranche-sur-Mer, FRANCE.

Corresponding author/maintainer: alexandre.schickele\@imev-mer.fr

## Introduction:

This document is oriented for expert users with knowledge in R programming, environmental genomics and ecology. The corresponding code is described in the following sections.

[...]

## Loading functions:

The code is structured in four modeling steps each of them loading a configuration file. The latter contains the list of packages, functions and general parameters necessary for each step. Below is a preview of some parameters that may be modified by an expert user:

```{r}

# --- Model specific parameters
NBOOST <- 3000 # maximum number of boosting rounds
N_FOLD <- 5 # number of k-fold cross validation runs
MAX_CLUSTER <- 20 # maximum number of clusters for parallel computing

NBOOTSTRAP <- 20 # number of bootstrap rounds for script 03
```

Once the general parameters are set, we load the functions corresponding to each step as following:

```{r}
source("01_query_data.R")
source("02a_model_param.R")
source("O2b_model_eval.R")
source("03a_bootstrap_predict.R")
```

[...]

## Building the data:

[...]

## Building the model:

[...]

## Evaluating the model:

[...]

## Building the spatial projections:

[...]

## Outputs:

[...]
