
# Testing fold validation for mtbr function

import numpy as np
from pandas import read_feather
path="~/workspace/bluecloud descriptor"
x = read_feather(path+"/data/X.feather").to_numpy()
y = read_feather(path+"/data/Y.feather").to_numpy()
n_val=5
x_lr=None

x=np.linspace(1,100, 100, dtype = int)
y=np.linspace(1,100, 100, dtype = int)
id=np.linspace(1,20, 20, dtype = int)

id=None
if (id is not None):
  print('hi')
  
x_tr, x_val = [x[:-id],x[-id:]]

randidx = np.random.permutation(x.shape[0])
n_val = int(0.2*x.shape[0])
tr_idx, val_idx = [randidx[:-n_val],randidx[-n_val:]]
x_tr, x_val = [x[tr_idx],x[val_idx]]
y_tr, y_val = [y[tr_idx], y[val_idx]]
x_lr_tr, x_lr_val = [x_lr[tr_idx], x_lr[val_idx]]

x = read_feather(path+"/data/X.feather").to_numpy()
id = np.linspace(1,20, 20, dtype = int)-1
x_val, x_tr = [x[id,], np.delete(x,id,0)]



def gbdtmo_fit(path)

  """
  Function to fit a MBTR model according to the train and test datasets generated
  by the 03_model.R script.
  A default custom loss function is implemented in the form of :
  
  Loss = mean((Y - Y_HAT/(sum(Y_HAT)/sum(Y)))^2)
  g <- Y - Y_HAT/(sum(Y_HAT)/sum(Y))
  h <- Y/Y
  
  parameter g and h are given as default in the function but can be changed
  if needed
  """
  
  from gbdtmo import load_lib, GBDTMulti
  import numpy as np
  from pandas import read_feather
  
  # Parameters
  path = "~/workspace/bluecloud descriptor"
  LIB = load_lib("/home/aschickele/.local/lib/python3.8/site-packages/GBDTMO/build/gbdtmo.so")

  #  Code
  X_tr = read_feather(path+"/data/X_tr.feather").to_numpy()
  Y_tr = read_feather(path+"/data/Y_tr.feather").to_numpy()
  
  X_te = read_feather(path+"/data/X_tr.feather").to_numpy()
  Y_te = read_feather(path+"/data/Y_tr.feather").to_numpy()
  
  p = {"max_depth": 5, "lr": 0.1, 'loss': b"mse"}
  m = GBDTMulti(LIB, out_dim=3, params = p)
  m.set_data((X_tr, Y_tr), (X_te, Y_te))
  
  def MSE_prop(y, y_hat)
    g = y - y_hat/(sum(y_hat)/sum(y))
    h = np.ones_like(y)
    return g, h
  
  g, h = MSE(m.preds_train.copy(), m.label.copy())
  m._set_gh(g, h)
  m.boost()
  
  return m
