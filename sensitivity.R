###############################################################################
#                             SENSITIVITY ANALYSIS                            #
# This script contains the following:                                         #
#     -  scripts needed for Table 3, 4 & 5                                    #
# This script depends on:                                                     #
#     - functions.R                                                           #
############################################################################### 
library(crayon)
library(tidyverse)
library(deSolve)
library(glmnet)


data <- readRDS('data.RDS')
source('functions.R')

mse.int = 1:75
train = 125:350       # days used for training
test = 350:450        # days used for testing
n = min(test)

T1 = 1:26
T2 = 26:51
T3 = 51:76

## 1: Classical SIR
cat(crayon::bold('Classical SIR\n'))
initial_vals = c(S = 17.4e6-data$I[n]-data$R[n], I = data$I[n], R = data$R[n])

cat(crayon::italic('changing beta\n'))
beta.range = c(0.132,0.121,0.1109516,0.099, 0.088)     # roughly +20%, +10%, actual, -10%, -20%
for (beta in beta.range){
  parameters = c(gamma = 0.09405225, beta = beta)
  output <- as.data.frame(ode(y=initial_vals, func = classical.SIR, parms=parameters, times = test)) %>% 
    pivot_longer(2:4, names_to = 'State')
  actual = unlist(data[test, 'I'])
  predict = unlist(subset(output, State == 'I')[, 'value'])
  cat('\t beta =', beta, 'gives MSE of', mse(actual[mse.int], predict[mse.int])/1e6, 'million\n')
}

cat(crayon::italic('\nchanging gamma\n'))
gamma.range = c(0.113,0.103,0.09405225,0.085, 0.075)     # roughly +20%, +10%, actual, -10%, -20%
for (gamma in gamma.range){
  parameters = c(gamma = gamma, beta = 0.1109516)
  output <- as.data.frame(ode(y=initial_vals, func = classical.SIR, parms=parameters, times = test)) %>% 
    pivot_longer(2:4, names_to = 'State')
  actual = unlist(data[test, 'I'])
  predict = unlist(subset(output, State == 'I')[, 'value'])
  cat('\t gamma =', gamma, 'gives MSE of', mse(actual[mse.int], predict[mse.int])/1e6, 'million\n')
}

cat(crayon::italic('\nchanging time-period\n'))
parameters = c(gamma = 0.09405225, beta = 0.1109516)
output <- as.data.frame(ode(y=initial_vals, func = classical.SIR, parms=parameters, times = test)) %>% 
  pivot_longer(2:4, names_to = 'State')
actual = unlist(data[test, 'I'])
predict = unlist(subset(output, State == 'I')[, 'value'])
cat('\tperiod T1 gives MSE of', mse(actual[T1], predict[T1])/1e6, 'million\n')
cat('\tperiod T2 gives MSE of', mse(actual[T2], predict[T2])/1e6, 'million\n')
cat('\tperiod T3 gives MSE of', mse(actual[T3], predict[T3])/1e6, 'million\n')



## 2: Time-Dependent SIR
cat(crayon::bold('\nTime-Dependent SIR\n'))
initial_vals = c(S = 17.4e6-data$I[n]-data$R[n], I = data$I[n], R = data$R[n])

J = K = 3

cat(crayon::italic('\nchanging J\n'))
J.range = 5:2
for (J.val in J.range){
  beta.ridge = doRidge(data$beta.t[train], seq(0,0.01,0.00002), J.val, length(test))
  gamma.ridge = doRidge(data$gamma.t[train], seq(0,0.01,0.00002), K, length(test))
  
  parameters = c('beta' = beta.ridge$param.hat, 'gamma' = gamma.ridge$param.hat)
  output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, parms=parameters, times = 1:length(test[-1]))) %>% 
    pivot_longer(2:4, names_to = 'State')
  
  actual = unlist(data[test[-length(test)], 'I'])
  predict = unlist(subset(output, State == 'I')[, 'value'])
  cat('\t J =', J.val, 'gives MSE of', mse(actual[mse.int], predict[mse.int])/1e6, 'million\n')
}

cat(crayon::italic('\nchanging K\n'))
K.range = 5:2
for (K.val in K.range){
  beta.ridge = doRidge(data$beta.t[train], seq(0,0.01,0.00002), J, length(test))
  gamma.ridge = doRidge(data$gamma.t[train], seq(0,0.01,0.00002), K.val, length(test))
  
  parameters = c('beta' = beta.ridge$param.hat, 'gamma' = gamma.ridge$param.hat)
  output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, parms=parameters, times = 1:length(test[-1]))) %>% 
    pivot_longer(2:4, names_to = 'State')
  
  actual = unlist(data[test[-length(test)], 'I'])
  predict = unlist(subset(output, State == 'I')[, 'value'])
  cat('\t K =', K.val, 'gives MSE of', mse(actual[mse.int], predict[mse.int])/1e6, 'million\n')
}

cat(crayon::italic('\nchanging time-period\n'))
beta.ridge = doRidge(data$beta.t[train], seq(0,0.01,0.00002), 3, length(test))
gamma.ridge = doRidge(data$gamma.t[train], seq(0,0.01,0.00002), 3, length(test))
parameters = c('beta' = beta.ridge$param.hat, 'gamma' = gamma.ridge$param.hat)
output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, parms=parameters, times = 1:length(test[-1]))) %>% 
  pivot_longer(2:4, names_to = 'State')

actual = unlist(data[test[-length(test)], 'I'])
predict = unlist(subset(output, State == 'I')[, 'value'])
cat('\tperiod T1 gives MSE of', mse(actual[T1], predict[T1])/1e6, 'million\n')
cat('\tperiod T2 gives MSE of', mse(actual[T2], predict[T2])/1e6, 'million\n')
cat('\tperiod T3 gives MSE of', mse(actual[T3], predict[T3])/1e6, 'million\n')


## 3: SIAR-T (I) model
cat(crayon::bold('\n SIAR-T (I)\n'))
initial_vals = c(S = 17.4e6-data$I[n]-data$R[n], Is = data$I[n]/2, 
                 A = data$I[n]/2, R = data$R[n])    

cat(crayon::italic('\nchanging beta_a\n'))
range.ba = c(0.144, 0.132,  0.12, 0.108, 0.096)
w.s = 0.3
coefs = coef(lm( log(beta.t) ~ I(-1*time), data = data[train,]))
beta.0 = exp(coefs[1]); omega = coefs[2]
beta.s = beta.0 * exp(-omega * test)
for (beta.a in range.ba){
  parameters = c('beta.s' = beta.s, 'gamma' = gamma.ridge$param.hat, 
                 'beta.a' = beta.a, w.s = w.s, w.a = 1-w.s)
  
  output <- as.data.frame(ode(y=initial_vals, func = asymptomatic.SIR, parms=parameters, times = 1:(length(test[-1])))) %>% 
    mutate(I = Is + A) %>% 
    pivot_longer(2:6, names_to = 'State')
  
  actual = unlist(data[test[-length(test)], 'I'])
  predict = unlist(subset(output, State == 'I')[, 'value'])
  cat('\tbeta_a =', beta.a, 'gives MSE of', mse(actual[mse.int], predict[mse.int])/1e6, 'million\n')
}

cat(crayon::italic('\nchanging w_s\n'))
range.ws = c(0.36, 0.33,  0.30, 0.27, 0.24)
beta.a = 0.12
coefs = coef(lm( log(beta.t) ~ I(-1*time), data = data[train,]))
beta.0 = exp(coefs[1]); omega = coefs[2]
beta.s = beta.0 * exp(-omega * test)
for (w.s in range.ws){
  parameters = c('beta.s' = beta.s, 'gamma' = gamma.ridge$param.hat, 
                 'beta.a' = beta.a, w.s = w.s, w.a = 1-w.s)
  
  output <- as.data.frame(ode(y=initial_vals, func = asymptomatic.SIR, parms=parameters, times = 1:(length(test[-1])))) %>% 
    mutate(I = Is + A) %>% 
    pivot_longer(2:6, names_to = 'State')
  
  actual = unlist(data[test[-length(test)], 'I'])
  predict = unlist(subset(output, State == 'I')[, 'value'])
  cat('\tw_s =', w.s, 'gives MSE of', mse(actual[mse.int], predict[mse.int])/1e6, 'million\n')
}

cat(crayon::italic('\nchanging time-period\n'))
parameters = c('beta.s' = beta.s, 'gamma' = gamma.ridge$param.hat, 
               'beta.a' = 0.12, w.s = 0.30, w.a = 1-0.30)
output <- as.data.frame(ode(y=initial_vals, func = asymptomatic.SIR, parms=parameters, times = 1:(length(test[-1])))) %>% 
  mutate(I = Is + A) %>% 
  pivot_longer(2:6, names_to = 'State')

actual = unlist(data[test[-length(test)], 'I'])
predict = unlist(subset(output, State == 'I')[, 'value'])
cat('\tperiod T1 gives MSE of', mse(actual[T1], predict[T1])/1e6, 'million\n')
cat('\tperiod T2 gives MSE of', mse(actual[T2], predict[T2])/1e6, 'million\n')
cat('\tperiod T3 gives MSE of', mse(actual[T3], predict[T3])/1e6, 'million\n')

## 4: SIAR-T (II) model
cat(crayon::bold('\n SIAR-T (II)\n'))

cat(crayon::italic('\nchanging w_s\n'))
range.ws = c(0.36, 0.33,  0.30, 0.27, 0.24)
for (w.s in range.ws){
  parameters = c('beta.s' = 1.5*beta.ridge$param.hat, 'gamma' = gamma.ridge$param.hat, 
                 'beta.a' = 0.75*beta.ridge$param.hat, w.s = w.s, w.a = 1-w.s)
  initial_vals = c(S = 17.4e6-data$I[n]-data$R[n], Is = data$I[n]*w.s, 
                   A = data$I[n]*(1-w.s), R = data$R[n])    
  
  output <- as.data.frame(ode(y=initial_vals, func = asymptomatic2.SIR, parms=parameters, times = 1:length(test[-1]))) %>% 
    mutate(I = Is + A) %>% 
    pivot_longer(2:6, names_to = 'State')
  
  actual = unlist(data[test[-length(test)], 'I'])
  predict = unlist(subset(output, State == 'I')[, 'value'])
  cat('\tw_s =', w.s, 'gives MSE of', mse(actual[mse.int], predict[mse.int])/1e6, 'million\n')
}

cat(crayon::italic('\nchanging time-period\n'))
w.s = 0.7
initial_vals = c(S = 17.4e6-data$I[n]-data$R[n], Is = data$I[n]*w.s, 
                 A = data$I[n]*(1-w.s), R = data$R[n]) 
parameters = c('beta.s' = beta.s, 'gamma' = gamma.ridge$param.hat, 
               'beta.a' = 0.12, w.s = w.s, w.a = 1-w.s)
output <- as.data.frame(ode(y=initial_vals, func = asymptomatic.SIR, parms=parameters, times = 1:(length(test[-1])))) %>% 
  mutate(I = Is + A) %>% 
  pivot_longer(2:6, names_to = 'State')

actual = unlist(data[test[-length(test)], 'I'])
predict = unlist(subset(output, State == 'I')[, 'value'])
cat('\tperiod T1 gives MSE of', mse(actual[T1], predict[T1])/1e6, 'million\n')
cat('\tperiod T2 gives MSE of', mse(actual[T2], predict[T2])/1e6, 'million\n')
cat('\tperiod T3 gives MSE of', mse(actual[T3], predict[T3])/1e6, 'million\n')

## 5: SIRV-T model
cat(crayon::bold('\n SIRV-T\n'))
nV = min(testV)  # starting point for vaccinated 
initial_vals = c(S = 17.4e6-data$I[nV]-data$R[nV]-data$vacc[nV], I = data$I[nV], R = data$R[nV],
                 V = data$vacc[nV])
trainV =  333:450     # days used for training of vaccinated data
testV = 450:479       # days used for testing of vaccinated data

betaV = doRidge(data$beta.t[trainV], seq(0,0.01,0.00002), 3, length(testV))
gammaV = doRidge(data$gamma.t[trainV], seq(0,0.01,0.00002), 3, length(testV))

coefs = coef(lm( log(alpha.t) ~ I(-1*time), data = data[trainV,]))
alpha.0 = exp(coefs[1]); omega = coefs[2]
alpha = alpha.0 * exp(-omega * testV)

cat(crayon::italic('\nchanging inefficiency sigma \n'))
sigma.vals = c(0.24, 0.22, 0.2, 0.18, 0.16)
for(sigma.val in sigma.vals){
  parameters = c('beta' = betaV$param.hat, 'gamma' = gammaV$param.hat, 
                 'alpha' = alpha, sigma = sigma.val)
  output <- as.data.frame(ode(y=initial_vals, func = vaccinated.SIR, parms=parameters, times = 1:length(testV[-1]))) %>% 
    pivot_longer(2:5, names_to = 'State')
  actual = unlist(data[testV[-length(testV)], 'I'])
  predict = unlist(subset(output, State == 'I')[, 'value'])
  cat('\tsigma =', sigma.val, 'gives MSE of', mse(actual, predict)/1e6, 'million\n')
}

cat(crayon::italic('\nchanging J \n'))
J.vals = 5:2
gammaV = doRidge(data$gamma.t[trainV], seq(0,0.01,0.00002), 3, length(testV))
for(J in J.vals){
  betaV = doRidge(data$beta.t[trainV], seq(0,0.01,0.00002), J, length(testV))
  parameters = c('beta' = betaV$param.hat, 'gamma' = gammaV$param.hat, 
                 'alpha' = alpha, sigma = 0.20)
  output <- as.data.frame(ode(y=initial_vals, func = vaccinated.SIR, parms=parameters, times = 1:length(testV[-1]))) %>% 
    pivot_longer(2:5, names_to = 'State')
  actual = unlist(data[testV[-length(testV)], 'I'])
  predict = unlist(subset(output, State == 'I')[, 'value'])
  cat('\tJ =', J, 'gives MSE of', mse(actual, predict)/1e6, 'million\n')
}

cat(crayon::italic('\nchanging K \n'))
K.vals = 5:2
betaV = doRidge(data$beta.t[trainV], seq(0,0.01,0.00002), 3, length(testV))

for(K in K.vals){
  gammaV = doRidge(data$gamma.t[trainV], seq(0,0.01,0.00002), K, length(testV))
  parameters = c('beta' = betaV$param.hat, 'gamma' = gammaV$param.hat, 
                 'alpha' = alpha, sigma = 0.20)
  output <- as.data.frame(ode(y=initial_vals, func = vaccinated.SIR, parms=parameters, times = 1:length(testV[-1]))) %>% 
    pivot_longer(2:5, names_to = 'State')
  actual = unlist(data[testV[-length(testV)], 'I'])
  predict = unlist(subset(output, State == 'I')[, 'value'])
  cat('\tK =', K, 'gives MSE of', mse(actual, predict)/1e6, 'million\n')
}








