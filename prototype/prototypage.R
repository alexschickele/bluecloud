
#' Testing some custom loss function on the y and y_hat

Y <- matrix(data = c(0,0.3,0.7,0.3,0,0.7,0,0.7,0.3),
            nrow = 3,
            ncol = 3)

Y_HAT_good <- Y*10

Y_HAT_bad <- matrix(data = rep((1/3),10),
                     nrow = 3,
                     ncol = 3)

Y_HAT_bad <- matrix(data = c(0,0.5,0.5,0.4,0,0.6,0,0.7,0.3),
                    nrow = 3,
                    ncol = 3)

sqrt(mean((Y - Y_HAT_bad)^2))
sqrt(mean((Y - Y_HAT_good)^2))

sqrt(mean((Y - Y_HAT_bad/(sum(Y_HAT_bad)/sum(Y)))^2))
sqrt(mean((Y - Y_HAT_good/(sum(Y_HAT_good)/sum(Y)))^2)) # works 

# NOW how to calculate gradient and hessian of this ? :D
g <- Y - Y_HAT_bad/(sum(Y_HAT_bad)/sum(Y))
h <- Y/Y




# PROTOTYPE DE FONCTION

git clone https://github.com/zzd1992/GBDTMO.git

from gbdtmo import load_lib, GBDTMulti, GBDTSingle
import numpy as np
from pandas import read_feather

LIB = load_lib("/home/aschickele/workspace/custom package/GBDTMO/build/gbdtmo.so")


inp_dim, out_dim = 10, 5
params = {"max_depth": 5, "lr": 0.1, 'loss': b"mse", "num_threads": 1}
booster = GBDTMulti(LIB, out_dim=out_dim, params=params)

inp_dim, out_dim = 5, 10
x_train, y_train = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)
x_valid, y_valid = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)


x_train, y_train = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)
x_valid, y_valid = np.random.rand(10000, inp_dim), np.random.rand(10000, out_dim)
booster.set_data((x_train, y_train), (x_valid, y_valid))






