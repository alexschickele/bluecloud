
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
