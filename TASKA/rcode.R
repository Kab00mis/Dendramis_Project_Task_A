rm(list=ls())
par(mfrow=c(1,3))
# =================================================================
# A
# =================================================================
# Q1: Simulation of DGPs
# =================================================================
T <- 1000 
sigma2 <- 1
c <- 0.1

# 1. MA(1)
theta1 <- 0.7
e <- rnorm(T, mean = 0, sd = sqrt(sigma2))
yy_ma1 <- numeric(T)

for (i in 2:T) {
  yy_ma1[i] <- c + e[i] + theta1 * e[i-1]
}
plot(yy_ma1, type = 'l', main = "Simulated MA(1)", ylab = "y(t)")

# 2. AR(2)
phi1 <- 0.5
phi2 <- 0.1
yy_ar2 <- numeric(T)
for (i in 3:T) {
  yy_ar2[i] <- c + phi1*yy_ar2[i-1] + phi2*yy_ar2[i-2] + e[i]
}
plot(yy_ar2, type='l', main="Simulated AR(2)", ylab="y(t)")

# 3. ARMA(2,1)
yy_arma <- numeric(T)
for (i in 3:T) {
  yy_arma[i] <- 0.1 + 0.5*yy_arma[i-1] + 0.1*yy_arma[i-2] + e[i] + 0.7*e[i-1]
}
plot(yy_arma, type='l', main="Simulated ARMA(2,1)", ylab="y(t)")

# =================================================================
# Q2a: Autocorrelation Function (ACF)
# =================================================================
autocorr_1 <- function(y, l) {
  T_len = length(y)
  yy_dev = y - mean(y)
  v = yy_dev[(l + 1):T_len]
  v_l = yy_dev[1:(T_len - l)]
  g_l = sum(v * v_l) / T_len
  rho_l = g_l / var(y)
  return(rho_l)
}

plot_my_acf <- function(data, title) {
  lags <- 20
  rho <- numeric(lags)
  for (l in 1:lags) { rho[l] <- autocorr_1(data, l) }
  plot(rho, type='h', main=title, xlab="Lag", ylab="ACF", ylim=c(-1,1))
  abline(h=0)
  abline(h = c(1.96/sqrt(T), -1.96/sqrt(T)), col = "blue", lty = 2)
}

par(mfrow=c(1,3))
plot_my_acf(yy_ma1, "ACF: MA(1)")
plot_my_acf(yy_ar2, "ACF: AR(2)")
plot_my_acf(yy_arma, "ACF: ARMA(2,1)")

plot_my_pacf <- function(data, title) {
  lags <- 20
  partial_rho <- numeric(lags)
  for (k in 1:lags) {
    df <- embed(data, k + 1)
    y <- df[, 1]
    x <- df[, 2:(k + 1)]
    model <- lm(y ~ x)
    partial_rho[k] <- coef(model)[k + 1]
  }
  plot(partial_rho, type='h', main=title, xlab="Lag", ylab="PACF", ylim=c(-1,1))
  abline(h=0) 
  abline(h = c(1.96/sqrt(length(data)), -1.96/sqrt(length(data))), col = "red", lty = 2)
}

par(mfrow=c(1,3))
plot_my_pacf(yy_ma1, "PACF: MA(1)")
plot_my_pacf(yy_ar2, "PACF: AR(2)")
plot_my_pacf(yy_arma, "PACF: ARMA(2,1)")

# =================================================================
# Q2b: AIC / BIC 
# =================================================================
# 1.OLS ESTIMATION
# MA-Simple model

out_ma3 <- lm(yy_ma1~1)
summary(out_ma3)
e_ma <- out_ma3$residuals
acf(e_ma)
e_lag1 <- c(NA, e_ma[1:NROW(e_ma)-1])
e_lag2 <- c(NA, NA, e_ma[1:(NROW(e_ma) - 2)])
e_lag3 <- c(NA, NA, NA, e_ma[1:(NROW(e_ma) - 3)])

#OLS
out_ma3 <- lm(yy_ma1~e_lag1 + e_lag2 + e_lag3)
summary(out_ma3)

#AR
y_lag1 <- c(NA, yy_ar2[1:NROW(yy_ar2)-1])
y_lag2 <- c(NA, NA, yy_ar2[1:(NROW(yy_ar2)-2)])
y_lag3 <- c(NA, NA, NA, yy_ar2[1:(NROW(yy_ar2) - 3)])

#OLS
out_ar3 <- lm(yy_ar2 ~ y_lag1 + y_lag2 + y_lag3, na.action = na.omit)
summary(out_ar3)


#Normal pdf
log_norm_pdf1 = function(x, m, ss) {
  (-1/2) * log(2 * pi * ss) - ((x - m)^2) / (2 * ss)
}

# 1. MA Models
likel_maq = function(theta, y, q) {
  T_len = length(y); const = theta[1]; ma_coefs = theta[2:(q+1)]; sig2 = theta[q+2]
  if(sig2 <= 0) return(1e10)
  u = numeric(T_len); lik = 0
  for (i in (q + 1):T_len) {
    m = const + sum(ma_coefs * u[(i-1):(i-q)])
    u[i] = y[i] - m
    lik = lik + log_norm_pdf1(y[i], m, sig2)
  }
  return(-lik)
}

results_ma <- matrix(0, 3, 2, dimnames=list(c("MA(1)","MA(2)","MA(3)"), c("AIC","BIC")))
thetas_start <-coef(out_ma3)

for (q in 1:3) {
  intercept <- thetas_start[1]
  ma_parts <- thetas_start[2:(q+1)] 
  sigma2 <- var(yy_ma1)
  theta0 = c(intercept, ma_parts, sigma2)
  fit = optim(theta0, likel_maq, y = yy_ma1, q = q, method = "BFGS")
  k = length(theta0)
  results_ma[q, ] = c(-2*(-fit$value) + 2*k, -2*(-fit$value) + k*log(T))
}

# 2. AR Models
likel_arp = function(phi, y, p) {
  T_len = length(y); const = phi[1]; phi_coefs = phi[2:(p+1)]; sig2 = phi[p+2]
  if(sig2 <= 0) return(1e10)
  lik = 0
  for (i in (p + 1):T_len) {
    m = const + sum(phi_coefs * y[(i-1):(i-p)])
    lik = lik + log_norm_pdf1(y[i], m, sig2)
  }
  return(-lik)
}

results_ar <- matrix(0, 3, 2, dimnames=list(c("AR(1)","AR(2)","AR(3)"), c("AIC","BIC")))
phi_start <- coef(out_ar3)

for(p in 1:3){
  intercept <- phi_start[1]
  ar_parts <- phi_start[2:(q+1)] 
  sigma2 <- var(yy_ar2)
  theta0 = c(intercept, ar_parts, sigma2)
  fit <- optim(theta0, function(t) likel_arp(t, yy_ar2, p), method = "BFGS")
  k <- length(theta0)
  results_ar[p, ] = c(-2*(-fit$value) + 2*k, -2*(-fit$value) + k*log(T))
}

# 3. ARMA(2,1)
likel_arma21 = function(theta, y) {
  const=theta[1]; p1=theta[2]; p2=theta[3]; t1=theta[4]; s2=theta[5]
  if(s2 <= 0) return(1e10)
  u = numeric(length(y)); lik = 0
  for (i in 3:length(y)) {
    m = const + p1*y[i-1] + p2*y[i-2] + t1*u[i-1]
    u[i] = y[i] - m
    lik = lik + log_norm_pdf1(y[i], m, s2)
  }
  return(-lik)
}
fit_arma = optim(c(0.1, 0.5, 0.1, 0.7, var(yy_arma)), likel_arma21, y=yy_arma, method="BFGS")
results_arma <- matrix(c(-2*(-fit_arma$value)+10, -2*(-fit_arma$value)+5*log(T)), 1, 2, 
                       dimnames=list("ARMA(2,1)", c("AIC","BIC")))

print(results_ma); print(results_ar); print(results_arma)

# =================================================================
# Q2d: Cross-Validation (End-of-Sample)
# =================================================================
split <- 800
# MA CV
mse_ma <- numeric(3)
for(q in 1:3) {
  f <- arima(yy_ma1[1:split], order=c(0,0,q))
  p <- predict(f, n.ahead=200)$pred
  mse_ma[q] <- mean((yy_ma1[(split+1):T] - p)^2)
}
#AR CV 
mse_ar <- numeric(3)
for(p in 1:3) {
  
  f <- arima(yy_ar2[1:split], order=c(p,0,0))
  preds <- predict(f, n.ahead=200)$pred
  mse_ar[p] <- mean((yy_ar2[(split+1):T] - preds)^2)
}


# ARMA CV
f_arma <- arima(yy_arma[1:split], order=c(2,0,1))
mse_arma <- mean((yy_arma[(split+1):T] - predict(f_arma, n.ahead=200)$pred)^2)

cat("\n--- MSE RESULTS ---\n")
print(mse_ma); print(mse_ar); cat("ARMA(2,1) MSE:", mse_arma, "\n")

# Q3. Estimated coefficients, and their CI.

summary_table <- function(fit_obj, param_names) {
  estimates <- fit_obj$par
  #Standard Errors από Hessian 
  std_errors <- sqrt(diag(solve(fit_obj$hessian)))
  # CI
  lower_ci <- estimates - 1.96 * std_errors
  upper_ci <- estimates + 1.96 * std_errors
  # Final table 
  res <- cbind(Estimate = estimates, 
               Std_Err = std_errors, 
               Lower_CI = lower_ci, 
               Upper_CI = upper_ci)
  rownames(res) <- param_names
  return(res)
}

# MA(1) 
fit_ma1 <- optim(c(0.1, 0.7, var(yy_ma1)), likel_maq, y = yy_ma1, q = 1, method = "BFGS", hessian = TRUE)
print(summary_table(fit_ma1, c("Constant", "MA(1) coeff", "Sigma^2")))

# AR(2) 
fit_ar2 <- optim(c(0.1, 0.1, 0.1, var(yy_ar2)), function(t) likel_arp(t, yy_ar2, 2), method = "BFGS", hessian = TRUE)
print(summary_table(fit_ar2, c("Constant", "Phi 1", "Phi 2", "Sigma^2")))

# ARMA(2,1) 
fit_arma <- optim(c(0.1, 0.1, 0.1, 0.1, var(yy_arma)), likel_arma21, y = yy_arma, method = "BFGS", hessian = TRUE)
print(summary_table(fit_arma, c("Constant", "Phi 1", "Phi 2", "Theta 1", "Sigma^2")))


# Q4
# ΓΙΑ MA
study_q4 = function(T_val, sig2_val) {
  
  # 1. Παραγωγή δεδομένων (DGP 1: MA1)
  e <- rnorm(T_val, mean = 0, sd = sqrt(sig2_val))
  y <- numeric(T_val)
  for (i in 2:T_val) { y[i] <- 0.1 + e[i] + 0.7 * e[i-1] }
  
  # 2. Έλεγχος AIC/BIC για MA(1), MA(2), MA(3)
  aics <- numeric(3)
  bics <- numeric(3)
  for (q in 1:3) {
    th0 = c(0.1, rep(0.1, q), var(y))
    fit = optim(th0, likel_maq, y = y, q = q, method = "BFGS")
    k = length(th0)
    aics[q] = 2*fit$value + 2*k
    bics[q] = 2*fit$value + k*log(T_val)
  }
  
  # Ποιο μοντέλο επιλέχθηκε; (Θέλουμε το 1)
  best_aic = which.min(aics)
  best_bic = which.min(bics)
  
  # 3. Cross-Validation (MSE) για το MA(1)
  split <- round(T_val * 0.8)
  f_cv <- arima(y[1:split], order=c(0,0,1))
  preds <- predict(f_cv, n.ahead = (T_val - split))$pred
  mse_val <- mean((y[(split+1):T_val] - preds)^2)
  
  # Επιστροφή αποτελεσμάτων
  return(c(T = T_val, Sigma2 = sig2_val, AIC_Best = best_aic, BIC_Best = best_bic, MSE = mse_val))
}

# Εκτέλεση για όλους τους συνδυασμούς
T_levels = c(200, 1000, 4000)
S_levels = c(0.1, 1, 5)
res_matrix = matrix(NA, 9, 5)
colnames(res_matrix) = c("T", "Sigma2", "Best_q_AIC", "Best_q_BIC", "MSE")

row = 1
for(t in T_levels) {
  for(s in S_levels) {
    res_matrix[row, ] = study_q4(t, s)
    row = row + 1
  }
}

print(res_matrix)

#ΓΙΑ AR
study_ar_q4 = function(T_val, sig2_val) {
  
  # 1. Παραγωγή δεδομένων (DGP 2: AR2)
  e <- rnorm(T_val, mean = 0, sd = sqrt(sig2_val))
  y <- numeric(T_val)
  for (i in 3:T_val) { y[i] <- 0.1 + 0.5*y[i-1] + 0.1*y[i-2] + e[i] }
  
  # 2. Έλεγχος AIC/BIC για AR(1), AR(2), AR(3)
  aics <- numeric(3)
  for (p in 1:3) {
    th0 = c(mean(y), rep(0.1, p), var(y))
    # Χρησιμοποιούμε τον wrapper για την likel_arp που φτιάξαμε πριν
    fit = optim(th0, function(t) likel_arp(t, y, p), method = "BFGS")
    aics[p] = 2*fit$value + 2*(p+2)
  }
  best_ar_aic = which.min(aics)
  
  # 3. Εκτίμηση παραμέτρων (για να δούμε την ακρίβεια του phi1=0.5)
  fit_final = optim(c(0.1, 0.5, 0.1, sig2_val), function(t) likel_arp(t, y, 2), method = "BFGS")
  phi1_hat = fit_final$par[2]
  
  return(c(T = T_val, Sigma2 = sig2_val, Best_p = best_ar_aic, Phi1_Est = phi1_hat))
}

# Εκτέλεση Loop
res_ar_matrix = matrix(NA, 9, 4)
colnames(res_ar_matrix) = c("T", "Sigma2", "Best_p_AIC", "Phi1_Estimate")

row = 1
for(t in T_levels) {
  for(s in S_levels) {
    res_ar_matrix[row, ] = study_ar_q4(t, s)
    row = row + 1
  }
}

print(res_ar_matrix)


# Q4
# ΓΙΑ MA
study_q4 = function(T_val, sig2_val) {
  
  # 1. Παραγωγή δεδομένων (DGP 1: MA1)
  e <- rnorm(T_val, mean = 0, sd = sqrt(sig2_val))
  y <- numeric(T_val)
  for (i in 2:T_val) { y[i] <- 0.1 + e[i] + 0.7 * e[i-1] }
  
  # 2. Έλεγχος AIC/BIC για MA(1), MA(2), MA(3)
  aics <- numeric(3)
  bics <- numeric(3)
  for (q in 1:3) {
    th0 = c(0.1, rep(0.1, q), var(y))
    fit = optim(th0, likel_maq, y = y, q = q, method = "BFGS")
    k = length(th0)
    aics[q] = 2*fit$value + 2*k
    bics[q] = 2*fit$value + k*log(T_val)
  }
  
  # Ποιο μοντέλο επιλέχθηκε; (Θέλουμε το 1)
  best_aic = which.min(aics)
  best_bic = which.min(bics)
  
  # 3. Cross-Validation (MSE) για το MA(1)
  split <- round(T_val * 0.8)
  f_cv <- arima(y[1:split], order=c(0,0,1))
  preds <- predict(f_cv, n.ahead = (T_val - split))$pred
  mse_val <- mean((y[(split+1):T_val] - preds)^2)
  
  # Επιστροφή αποτελεσμάτων
  return(c(T = T_val, Sigma2 = sig2_val, AIC_Best = best_aic, BIC_Best = best_bic, MSE = mse_val))
}

# Εκτέλεση για όλους τους συνδυασμούς
T_levels = c(200, 1000, 4000)
S_levels = c(0.1, 1, 5)
res_matrix = matrix(NA, 9, 5)
colnames(res_matrix) = c("T", "Sigma2", "Best_q_AIC", "Best_q_BIC", "MSE")

row = 1
for(t in T_levels) {
  for(s in S_levels) {
    res_matrix[row, ] = study_q4(t, s)
    row = row + 1
  }
}

print(res_matrix)

#ΓΙΑ AR
study_ar_q4 = function(T_val, sig2_val) {
  
  # 1. Παραγωγή δεδομένων (DGP 2: AR2)
  e <- rnorm(T_val, mean = 0, sd = sqrt(sig2_val))
  y <- numeric(T_val)
  for (i in 3:T_val) { y[i] <- 0.1 + 0.5*y[i-1] + 0.1*y[i-2] + e[i] }
  
  # 2. Έλεγχος AIC/BIC για AR(1), AR(2), AR(3)
  aics <- numeric(3)
  for (p in 1:3) {
    th0 = c(mean(y), rep(0.1, p), var(y))
    # Χρησιμοποιούμε τον wrapper για την likel_arp που φτιάξαμε πριν
    fit = optim(th0, function(t) likel_arp(t, y, p), method = "BFGS")
    aics[p] = 2*fit$value + 2*(p+2)
  }
  best_ar_aic = which.min(aics)
  
  # 3. Εκτίμηση παραμέτρων (για να δούμε την ακρίβεια του phi1=0.5)
  fit_final = optim(c(0.1, 0.5, 0.1, sig2_val), function(t) likel_arp(t, y, 2), method = "BFGS")
  phi1_hat = fit_final$par[2]
  
  return(c(T = T_val, Sigma2 = sig2_val, Best_p = best_ar_aic, Phi1_Est = phi1_hat))
}

# Εκτέλεση Loop
res_ar_matrix = matrix(NA, 9, 4)
colnames(res_ar_matrix) = c("T", "Sigma2", "Best_p_AIC", "Phi1_Estimate")

row = 1
for(t in T_levels) {
  for(s in S_levels) {
    res_ar_matrix[row, ] = study_ar_q4(t, s)
    row = row + 1
  }
}

print(res_ar_matrix)

#ARMA
# ΓΙΑ ARMA
study_arma_q4 = function(T_val, sig2_val) {
  
  # 1. Παραγωγή δεδομένων (DGP 3: ARMA(2,1))
  # y[i] = 0.1 + 0.5*y[i-1] + 0.2*y[i-2] + e[i] + 0.7*e[i-1]
  e <- rnorm(T_val, mean = 0, sd = sqrt(sig2_val))
  y <- numeric(T_val)
  for (i in 3:T_val) { 
    y[i] <- 0.1 + 0.5*y[i-1] + 0.2*y[i-2] + e[i] + 0.7*e[i-1] 
  }
  
  # 2. Έλεγχος AIC/BIC για 3 συνδυασμούς: (1,1), (2,1), (2,2)
  # Ορίζουμε τους συνδυασμούς (p, q) που θα τεστάρουμε
  orders <- list(c(1,1), c(2,1), c(2,2))
  aics <- numeric(3)
  bics <- numeric(3)
  
  for (i in 1:3) {
    p_curr <- orders[[i]][1]
    q_curr <- orders[[i]][2]
    
    th0 = c(mean(y), rep(0.1, p_curr), rep(0.1, q_curr), var(y))
    
    # Χρήση της likel_arma_pq 
    fit = optim(th0, likel_arma_pq, y = y, p = p_curr, q = q_curr, method = "BFGS")
    
    k = length(th0)
    aics[i] = 2*fit$value + 2*k
    bics[i] = 2*fit$value + k*log(T_val)
  }
  
  # Ποιο μοντέλο επιλέχθηκε; (Αντιστοιχία: 1->(1,1), 2->(2,1), 3->(2,2))
  best_aic_idx = which.min(aics)
  best_bic_idx = which.min(bics)
  
  # 3. Εκτίμηση παραμέτρων (για το αληθές μοντέλο ARMA(2,1))
  # Θέλουμε να δούμε πόσο κοντά πέφτει το AR1 (0.5) και το MA1 (0.7)
  fit_final = optim(c(0.1, 0.5, 0.2, 0.7, sig2_val), 
                    likel_arma_pq, y = y, p = 2, q = 1, method = "BFGS")
  
  ar1_hat = fit_final$par[2]
  ma1_hat = fit_final$par[4]
  
  return(c(T = T_val, Sigma2 = sig2_val, 
           Best_Model_AIC = best_aic_idx, 
           AR1_Est = ar1_hat, MA1_Est = ma1_hat))
}

# 4. Εκτέλεση Loop για όλους τους συνδυασμούς
res_arma_matrix = matrix(NA, 9, 5)
colnames(res_arma_matrix) = c("T", "Sigma2", "Best_Idx_AIC", "AR1_Hat", "MA1_Hat")

row = 1
for(t in T_levels) {
  for(s in S_levels) {
    res_arma_matrix[row, ] = study_arma_q4(t, s)
    row = row + 1
  }
}

print(res_arma_matrix)

# =================================================================
#B
# =================================================================
M = 100 # Αριθμός επαναλήψεων Monte Carlo
theta1_true = 0.7     
c_true = 0.1

# Μετρητές επιτυχίας (Success Counters)
# Θα μετράμε πόσες φορές το κάθε κριτήριο επιλέγει το σωστό μοντέλο (q=1)

wins_aic = 0
wins_bic = 0
wins_cv = 0

#Η ΛΟΥΠΑ MONTE CARLO (100 ΕΠΑΝΑΛΗΨΕΙΣ)
for (m in 1:M) {
  # Παραγωγή νέων δεδομένων MA(1) σε κάθε επανάληψη
  u = rnorm(T, mean = 0, sd = sqrt(sigma2))
  yy_sim = matrix(0, T, 1)
  for (i in 2:T) {
    yy_sim[i] = c_true + u[i] + theta1_true * u[i-1]
  }
  # Προσωρινοί πίνακες για την αποθήκευση των τιμών των κριτηρίων για q=1,2,3
  temp_aic = numeric(3)
  temp_bic = numeric(3)
  temp_mse = numeric(3)
  
  # Διαχωρισμός 80/20 για το Cross-Validation
  n_train = round(0.8 * T)
  y_train = yy_sim[1:n_train]
  y_test = yy_sim[(n_train + 1):T]
  
  
  # Loop για την εκτίμηση των 3 ανταγωνιστικών μοντέλων (q=1, q=2, q=3)
  for (q in 1:3) {
    # Α. Υπολογισμός AIC & BIC (στο πλήρες δείγμα)
    theta0 = c(0.1, rep(0.1, q), var(yy_sim))
    fit = optim(theta0, likel_maq, y = yy_sim, q = q, method = "BFGS") #βρισκει το μαξ θ απο likelihood
    log_lik = -fit$value
    k = length(theta0)
    temp_aic[q] = -2 * log_lik + 2 * k
    temp_bic[q] = -2 * log_lik + k * log(T)
    # Β. Υπολογισμός Cross-Validation (80/20)
    # Εκτίμηση παραμέτρων μόνο στο 80% (Training Set)
    theta0_cv = c(0.1, rep(0.1, q), var(y_train))
    fit_cv = optim(theta0_cv, likel_maq, y = y_train, q = q, method = "BFGS")
    p = fit_cv$par # Οι εκτιμημένες παράμετροι από το training
    # Πρόβλεψη στο Test Set (20%)
    # Υπολογίζουμε τα σφάλματα u για όλη τη σειρά yy_sim χρησιμοποιώντας τις παραμέτρους του training
    u_full = matrix(0, T, 1)
    for (t in (q + 1):T) {
      ma_part = sum(p[2:(q+1)] * u_full[(t-1):(t-q)])
      u_full[t] = yy_sim[t] - p[1] - ma_part
    }
    # Το MSE υπολογίζεται από τα σφάλματα στο κομμάτι που "δεν είδε" το μοντέλο (το τελευταίο 20%)
    errors_test = u_full[(n_train + 1):T]
    temp_mse[q] = mean(errors_test^2)
  }
  
  # Καταγραφή Νικητών
  # Αν το ελάχιστο κριτήριο δείχνει το μοντέλο 1 (q=1), αυξάνουμε τον μετρητή
  if (which.min(temp_aic) == 1) wins_aic = wins_aic + 1
  if (which.min(temp_bic) == 1) wins_bic = wins_bic + 1
  if (which.min(temp_mse) == 1) wins_cv = wins_cv + 1
}


summary_results <- data.frame(
  Criterion = c("AIC", "BIC", "Cross-Validation (80/20)"),
  Success_Rate_Percentage = c(wins_aic, wins_bic, wins_cv)
)


print(summary_results)