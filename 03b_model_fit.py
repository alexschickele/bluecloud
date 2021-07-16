def mbtr_fit(path, loss_type='mse_custom', n_boosts=200, learning_rate=0.1, min_leaf=100, lambda_weights=0.001, alphas=None, n_harmonics=None, x_lr=None):

  """
  Function to fit a MBTR model according to the train and test datasets generated
  by the 03_model.R script.
  A custom loss function still needs to be implemented in the losses.py script
  of the mbtr module.
  """

  import numpy as np
  import mbtr.utils as ut
  from mbtr.mbtr import MBT
  from mbtr.utils import set_figure
  from pandas import read_feather
  
  X_tr = read_feather(path+"/data/X_tr.feather").to_numpy()
  Y_tr = read_feather(path+"/data/Y_tr.feather").to_numpy()
  
  m = MBT(loss_type = loss_type, n_boosts=n_boosts, learning_rate=learning_rate, 
  min_leaf=min_leaf,lambda_weights=lambda_weights, verbose = 1, alphas=alphas, 
  n_harmonics=n_harmonics, early_stopping_rounds=3).fit(X_tr, Y_tr, x_lr=x_lr, do_plot=True)
  
  X_te = read_feather(path+"/data/X_te.feather").to_numpy()
  
  y_hat = m.predict(X_te)
  
  return m, y_hat


  
  
