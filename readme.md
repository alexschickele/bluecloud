# Demonstrator 2 - Notebook 2:

# Projecting plankton functional gene clusters

------------------------------------------------------------------------

A. Schickele, P. Debeljack, L. Guidi, S. Ayata, J.O. Irisson

**Corresponding author/maintainer:** `alexandre.schickele@imev-mer.fr`

# 1. Introduction:

The genetic potential of plankton is influenced by the environmental context. In this notebook, we provide a series of tools to explore the relationship between the abundance of plankton genes and the environmental context in order to project the biogeography of key metabolic pathways and as yet unknown plankton genes.

This notebook is divided in two services, designed for expert and non-expert users, respectively:

-   service 1: R modeling framework. This service describes the entire modeling pipeline and codes. It is essentially designed for expert users with a strong genomic and modeling background.

-   service 2: R Shiny application. This services provides a ready-to-use tool to explore the biogeography of plankton genetic diversity related to key metabolic pathways. It is designed for general users and easily accessible at the following link: `https://`

The code is accessible upon registration at: `https://`

------------------------------------------------------------------------

# 2. Data sources

## 2.1. Environmental data

First, we retrieved a large set of annual climatologies from World Ocean Atlas and CMEMS. These data encompasses the 2005 to 2017 period, at 1° x 1° resolution and on a geographical domain ranging from -180 to +180°E and -90 to +90°. The environmental variables included are:

-   **AOU**: Apparent Oxygen Utilization

-   **nitrate**: Nitrate concentration

-   **o2sat**: Percent Oxygen Saturation

-   **oxygen**: Dissolved oxygen concentration

-   **phosphate**: Phosphate concentration

-   **salinity**: Sea surface salinity

-   **silicate**: Silicate concentration

-   **temperature**: Sea surface temperature

-   **MLD**: Mixed Layer Depth, calculated using the method from Boyer de Montégut (2004)

-   **CHL**: Chlorophyll A concentration

-   **ZE**: depth of the lower limit of the Euphotic Zone, calculated using the method from Lavigne et al. (2012)

-   **bathymetry**: sea bottom bathymetry

-   **distcoast**: calculated distance to the nearest coast

For all environmental variables, we calculated the mean (`mean`), standard deviation (`sd`), median (`med`), mean average distance (`mad`), monthly minimum (`min`) and monthly maximum (`max`) over the 2005 - 2017 period.

## 2.2. Genomic data

Second, we use MetaGenomic data retrieved from the Marine Atlas of Tara Ocean Unigenes (MATOU) described in the Notebook 1. These data correspond to the number of reads (i.e. proxy of abundance) per genes and Tara Ocean station, for the surface and 0.5 to 3 µm organisms (i.e. pico-plankton). To be able to explore the genomic diversity of plankton, the data are annotated as following:

-   If possible, each gene is functionally annotated using the KEGG (Kyoto Encyclopedia of Genes and Genomes) framework: <https://www.genome.jp/kegg/pathway.html>. This correspond to a two-level annotation, i.e. the KEGG functional annotation, contained in KEGG metabolic pathways.

-   Genes are grouped by similarity of sequences into clusters, containing annotated genes (i.e. KEGG + pathways) and as yet unknown genes.

In other words, we are able to model the biogeography of gene clusters related to a given metabolic pathways, and then the genes related to one of the clusters.

[...]

------------------------------------------------------------------------

# 3. **Service 1:** R modeling framework for expert users

## 3.1. Framework structure

### 3.1.1. Folder

The files are structured in different folders as following:

-   the root contains this document that serves as a master script for service 1 and the `app.R` file to launch the Shiny application from R.

-   the `/code` folder contains all scripts corresponding to each step of the R modeling framework, as well as the `ui.R`, `server.R` and `global.R` files necessary to run the Shiny application from R.

-   the `/function` folder corresponds to custom function necessary during the R modeling framework

-   the `/data` folder corresponds to temporary data file written during the modeling process

### 3.1.2. Functions

The code is structured in four modeling steps each of them loading a configuration file: `00a_config.R`.

The latter contains the list of packages, functions and general parameters necessary for each step. Below is a preview of some parameters that may be modified by an expert user:

```{r}

# --- Model specific parameters
NBOOST <- 3000 # maximum number of boosting rounds
N_FOLD <- 5 # number of k-fold cross validation runs
MAX_CLUSTER <- 20 # maximum number of clusters for parallel computing

NBOOTSTRAP <- 20 # number of bootstrap rounds for script 03
```

Once the general parameters are set, we can load the functions corresponding to each step as following:

```{r}

source("./code/01_query_data.R")
source("./code/02a_model_param.R")
source("./code/02b_model_eval.R")
source("./code/03a_bootstrap_predict.R")
```

[...]

## 3.2. Selecting the genomic and environmental data

### 3.2.1. General principle

Now that we described the data accessible by the notebook, we can select which part to select as input for the model using the `query_data()` function.

This first step sends a query to an SQL database that contains the genomic data described in section 2.2. as well as the environmental data corresponding to each TARA Ocean station. The `query_data()` function extracts:

-   the genomic data corresponding to a user-defined metabolic pathway to model, its underlying gene clusters

-   the environmental data at the necessary TARA Ocean stations corresponding to the selected variables

The `query_data()` function is called as following:

```{r}

query <- query_data(config_file = "./code/00a_config.R",
                    KEGG_p = "00195",
                    CLUSTER_SELEC = list(MIN_STATIONS = 80, MIN_GENES = 5, MAX_GENES = 25),
                    ENV_METRIC = c("mean","sd","dist","bathy"))
```

### 3.2.2. Inputs

The function take as user-defined input parameters:

-   `config_file`: the location of the config file `00a_config.R`

-   `KEGG_p`: a string corresponding to the ID of the KEGG pathway to select. The ID are available at <https://www.genome.jp/kegg/pathway.html>

-   `CLUSTER_SELECT`: a list of integers defining the size of the clusters to consider within a metabolic pathways: minimum number of stations where it is present, minimum and maximum number of genes per clusters (e.g. discarding clusters of two-genes or ubiquitous clusters)

-   `ENV_METRIC`: a vector of strings containing the environmental variables to include in the features. The strings correspond to the variable names or format described in section 2.1.

### 3.2.3. Outputs

The function returns two outputs, written as `X.feather` and `Y.feather` tables in the `/data` folder and as a list of dataframes in the working environment (i.e. needed for the Shiny application; Service 2). The two output have the following format:

-   **X**: *n-environmental variables x n-stations* table of environmental values, later referred as "**features**"

-   **Y**: *n-clusters x n-stations* table of reads, later referred as "**targets**". To alleviate spurious model results, the reads are normalized per stations and re-scaled globally between 0 and 1. The latter will be referred as "**normalized abundances**".

In addition, the size of the outputs (i.e. number of stations, targets and features) are printed in the console at the end of the query.

[...]

## 3.3. Training the model

### 3.3.1. General principle

The machine learning algorithm used in this notebook is based on a Python library (MBTR; see <https://mbtr.readthedocs.io/en/latest/>) and performs multivariate gradient boosting (i.e. several target modeled at once). The underlying Python implementation is invisible to the user.

In line with best practices in machine learning, the `X.feather` and `Y.feather` tables retrieved from section 3.2. are split between a training set and an test set. The latter is performed using *n*-fold cross-validation splits: i.e. the data are split in *n*-equal size groups, *n*-models are trained on *n-1* splits, holding the last split for model evaluation only. The test split is different for each of the *n*-models.

In order to select the best model, we tested different set of hyperparameters, leading to one algorithm per cross-validation fold and set of hyperparameters. The fit of each algorithm is measure by a loss function at each boosting round. **The best model is selected according to the hyperparameters and number of boosting round that produced the minimum loss averaged between all cross-validation folds.**

The model training and hyperparameters selection is called using the `model_run()` function:

```{r}

run <- model_run(config_file = "./code/00a_config.R",
                 HYPERPARAMETERS = data.frame(LEARNING_RATE = c(1e-2, 1e-2, 1e-2, 1e-2),
                                              N_Q = c(5, 10, 20, 50),
                                              MEAN_LEAF = c(30, 40, 50, 60)),
                 relative = TRUE,
                 verbose = TRUE)
```

### **3.3.2. Inputs**

This function uses the following user-defined input parameters:

-   `config_file`: the location of the config file `00a_config.R`

-   `HYPERPARAMETERS`: numeric dataframe of hyperparameters to be tested. With `LEARNING_RATE` numeric vector of learning rates to test, `N_Q` integer vector of quantiles to use for tree splitting, `MEAN_LEAF` : minimum number of data per leaf (i.e. after splitting). For more detailed informations, see <https://mbtr.readthedocs.io/en/latest/mbtr.html#module-mbtr.mbtr>.

-   `relative`: if `TRUE`, the targets are calculated as relative abundance between each clusters. This option is better suited to explore the genomic composition of plankton.

-   `verbose`: if `TRUE`, a plot of the loss per boosting round and the best set of hyperparameters is displayed.

Supplementary model input parameters can be found in the `00a_config.R` file as well as in the `mbtr_fit()` function within the function described here. **However, we do not recommend modifying these parameters unless you are familiar with the MBTR Python library and aware of the computing requirements.** Please note that depending on the number of boosting rounds, MBTR may take a few minutes to run (c.a. 10min for 1000 boosting rounds). A progress bar is displayed in the console to inform the user.

### **3.3.3. Outputs**

The function returns two outputs, written as `m.pickle` and `HYPERPARAMETERS.feather` in the `/data folder`. The outputs are described hereafter:

-   `m.pickle`: output file from the MBTR library corresponding to the selected model.

-   `HYPERPARAMETERS.feather`: dataframe of the same format as `HYPERPARAMETERS` input. However, only the line corresponding to the best set of hyperparameters is kept. A supplementary columns informing on the optimal boosting rounds has also been added.

[...]

## 3.4. Evaluating the model

### 3.4.1. General principle

As described in the previous section, the performance of our select model is evaluated by using the test split corresponding to each cross-validation fold. In this step, we perform three evaluation metrics: i.e. r-squared (`R2`), mean squared error (`mse`) and root mean squared error (`rmse`). For each cross-validation fold, these metrics compare the prediction of our model on `X_te.feather` (referred as "y_hat" )with the true values written in `Y_te.feather`.

To evaluate the performance of our selected model, we call the `model_eval()` function:

```{r}

eval <- model_eval(config_file = "./code/00a_config.R",
                   var_importance = FALSE)
```

### 3.4.2. Inputs

This function takes the following parameters as input:

-   `config_file`: the location of the config file `00a_config.R`

-   `var_importance`: if `TRUE`, the importance of each environmental variable to describe the distribution of genomic data is computed. This is done by counting the number of times each variable is used to split the trees. Please note that this features is time and computing intensive and should not be used by default for prototyping.

The number of user-defined inputs is relatively limited as this function directly reads the `m.pickle`, `X_te.feather` and `Y_te.feather` on the disk.

### 3.4.3. Outputs

This function provides a series of graphical outputs and evaluation metrics listed below:

-   Variable importance plot: this plot provides the importance of each environmental variable in explaining the relationship between environmental patterns and abundance of gene clusters.

-   Prediction vs test set plot: for each cross-validation fold, this function provides a plot comparing the true test values (`Y_te.feather`) of the targets and the model predictions on this test set (`y_hat`).

-   Evaluation metrics: as described in section 3.4.1., the average `R2`, `mse` and `rmse` and the corresponding standard deviation between cross-validation folds are printed in the console.

## 3.5. Spatial projections

### 3.5.1. General principle

Now that we selected the best hyperparameter and evaluated the corresponding MBTR model, we are still missing the spatial projections of the modeled response of gene clusters to environmental patterns. In line with best practices in machine learning and biogeography, we need to project both the average environmental response and an estimation of its uncertainty.

Therefore, this last step performs *m*-bootstrap round projections. At each round, the entire data (i.e. training and test sets) is randomly re-sampled with replacement. The *m*-projections are then represented by means of their average and standard deviation.

The spatial projections are calculated by the `model_proj()` function:

```{r}

proj <- model_proj(config_file = "./code/00a_config.R",
                   ENV_METRIC = c("mean","sd","med","mad","dist","bathy"))
```

Before describing the inputs and outputs of `model_proj()`, we need to understand how the underlying `colmat()` function works. In order to plot both average projection and uncertainty on the same plot, we required saturation-dependent color palettes, that are defined by:

-   the color (i.e. chroma) , that represents the different abundance levels (e.g. red for low abundance levels, blue for high abundance levels)

-   the color saturation, that represents the different uncertainty levels (i.e. lower color saturation reflects high uncertainty between bootstrap runs).

The `colmat()` function therefore builds a *n*-color level x *m*-saturation level matrix, with a unique index corresponding to each color x saturation combination. The projections retrieved from `model_proj()` do not represent the cluster abundance or uncertainty but the combination of both through this color x saturation index.

### 3.5.2. Inputs

The `model_proj()` function reads the `X.feather`, `Y.feather` and `HYPERPARAMETERS.feather` files written in section 3.2.3 and 3.3.3. respectively. In addition, it takes as input the following user-defined parameters:

-   `config_file`: the location of the config file `00a_config.R`

-   `ENV_METRIC`: a vector of strings containing the environmental variables to include in the features. The strings correspond to the variable names or format described in section 2.1. **To avoid spurious projections, please note that this parameter needs to contain the same variable names as in section 3.2.2. (i.e. the variables on which the model has been trained).**

### 3.5.3. Outputs

The `model_proj()` function provides a series of outputs stored in a list and necessary for graphical outputs. The output list has the following structure:

-   `col_matrix`: a *n*-color level x *m*-saturation level matrix containing the corresponding color references

-   `cutx`: a numeric vector corresponding to the values of the different uncertainty levels

-   `cuty`: a numeric vector corresponding to the values of the different abundance levels

<!-- -->

-   `proj`: a raster stack with each layers containing the abundance x uncertainty projections for one target (i.e. gene cluster)

-   `col`: a list of color vectors, each corresponding to the abundance x uncertainty index of one `proj$layer`. This output could be omitted and re-calculated from `col_matrix`, with unnecessary computing cost, however.

### 3.5.4. Graphical outputs

To avoid unnecessary computing, the graphical outputs of this notebook are provided by separate functions using the outputs described in section 3.5.3., namely `legend_proj()`, `map_proj()` and `cor_proj()`.

The `legend_proj()` function plot the saturation-dependent color palette that is used in the projection maps. It takes the following fixed inputs:

```{r}

legend_proj(col_matrix = proj$col_matrix, cutx = proj$cutx, cuty = proj$cuty)
```

The `map_proj()` function plots the targets called by the user as following:

```{r}

# Display plot for all targets:
map_proj(proj = proj$proj, col = proj$col, targetID = seq(1:nlayers(proj)))

# Display plot a specific target (e.g. target n°1):
map_proj(proj = proj$proj, col = proj$col, targetID = 1)
```

This function takes the `proj` and `col` objects from section 3.5.3. and a user defined `targetID`:

-   `targetID`: numeric vector of target numbers to plot. Accept one or several numeric values ranging from 1 to the maximum number of targets

Finally, the `cor_proj()` function provide a spatial correlation plot of all targets, based on the pearson correlation coefficient. It takes the `proj` object as input.

------------------------------------------------------------------------

# 4. **Service 2:** R Shiny application for non-expert users
