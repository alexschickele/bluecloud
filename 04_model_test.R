# =========================== 3. TEST FINAL MODEL ==============================

m <- mbtr_test(path=input.wd, loss_type='mse', min_leaf=10, 
               learning_rate=best_hp$LEARNING_RATE,
               lambda_weights=best_hp$LAMBDA_WEIGHTS,
               lambda_leaves=best_hp$LAMBDA_LEAVES)

# --- Does the lines sum to 1 ?
range(apply(m[[2]], 1, sum))

# --- Plot

plot(Y_te[,1], type = 'l', ylim = c(0,max(Y)), main = "Y_te (dotted) VS Y_hat (full)", col = "red", lty = "dotted", lwd=2)
lines(Y_te[,2], col = "blue", lty = "dotted", lwd=2)
lines(Y_te[,3], col = "green", lty = "dotted", lwd=2)

lines(m[[2]][,1], col = "red")
lines(m[[2]][,2], col = "blue")
lines(m[[2]][,3], col = "green")