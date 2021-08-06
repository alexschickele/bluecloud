def mbtr_fit(path, loss_type: str='mse', n_boosts: int=1000, n_q: int=10, learning_rate: float=0.1, 
val_ratio: float=0.2, val_path: str='~/workspace/bluecloud descriptor',
min_leaf: int=20, lambda_weights: float=0.0001, lambda_leaves: float=0.0001):

  """
  Fit a MBTR model under mean square error loss function. The inputs corresponds
  to the MBT().fit() function of the MBTR module except:
    
  :input x_tr: training set generated by the 03_model_param.R script
  :input y_tr: training set generated by the 03_model_param.R script
  
  (optionnal)
  :input x_val: external validation set generated by the 03_model_param.R script
  :input x_val: external validation set generated by the 03_model_param.R script
  
  :param val_path: path to the validation inputs x_val and y_val as generated by the
                   03_model_param.R script. This is necessary in the MBTR 1.1 release
                   as there is a bug concerning the size of x_tr that needs to be equal
                   to x. Therefore we give x = x_tr to MBT().fit() and then load external
                   x_val and y_val datasets.
  
  This function returns a list containing:
  [[1]] the MBT model object
  [[2]] a list of the loss computed on the validation datasets for each boosting round
  
  See the MBT().fit() and https://mbtr.readthedocs.io/en/latest/mbtr.html#module-mbtr.mbtr
  for more details on parameters.
  """

  import numpy as np
  import mbtr.utils as ut
  from mbtr.mbtr import MBT
  from mbtr.utils import set_figure
  from pandas import read_feather
  
  X_tr = read_feather(path+"/data/X_tr.feather").to_numpy()
  Y_tr = read_feather(path+"/data/Y_tr.feather").to_numpy()
  
  m = MBT(loss_type = loss_type,
          n_boosts=n_boosts, 
          n_q=n_q,
          learning_rate=learning_rate, 
          min_leaf=min_leaf,
          lambda_weights=lambda_weights, lambda_leaves=lambda_leaves,
          verbose = 0, 
          val_ratio = val_ratio,
          val_path= val_path,
          early_stopping_rounds=n_boosts).fit(X_tr, Y_tr, do_plot=True)
  
  return m

def mbtr_predict(model, X_pred):

  """
  Predict the y_hat target dataframe according to a x_pred feature matrix and
  a previously fitted MBTR model.
  
  :param model: a MBT model object
  :param X_pred: a Nfeature * Nobs matrix on which the model will predict the
                 target values y_hat
                 
  See the MBT().fit() and https://mbtr.readthedocs.io/en/latest/mbtr.html#module-mbtr.mbtr
  for more details on parameters.
  """

  import numpy as np
  import mbtr.utils as ut
  from mbtr.mbtr import MBT
  from mbtr.utils import set_figure
  
  X_pred = X_pred.to_numpy()
  y_hat = model.predict(X_pred)
  
  return y_hat
  
  
