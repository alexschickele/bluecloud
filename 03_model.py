
from gbdtmo import load_lib
LIB = load_lib("/home/aschickele/.local/lib/python3.8/site-packages/GBDTMO/build/gbdtmo.so")


"""
This script uses X.feather and Y.feather as input for a multivariate
boosted regression tree model.
The MBTR model is run from python as the recitulate package does not work
for the syntax of the fitting function
See MBTR_installation.Rmd for more information on the custom loss function that
was implemented
"""

input.wd = "~/workspace/bluecloud descriptor/data"

# --- Loading Python modules
import numpy as np
import mbtr.utils as ut
from mbtr.mbtr import MBT
from mbtr.utils import set_figure
from pandas import read_feather

# --- Load data
X = read_feather(input.wd+"/X.feather").to_numpy()
Y = read_feather(input.wd+"/Y.feather").to_numpy()
N = len(X)

# --- Plot raw data ?

# --- Train and test datasets
ID.tr <- sample(seq(1,N), 0.7*N)
ID.te <- seq(1,N)[-ID.tr]

X.tr <- X[ID.tr,]
X.te <- X[ID.te,]

Y.tr <- Y[ID.tr,]
Y.te <- Y[ID.te,]

# --- Model fit
m = MBT(loss_type = 'mse_custom', n_boosts=30,  min_leaf=100, lambda_weights=1e-3).fit(X, Y, do_plot=True)
