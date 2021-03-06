###############################################################################
#                             MAIN ANALYSIS                                   #
# This script contains the following:                                         #
#     - implementation of SIR, SIR-T, SIAR-T (I & II), SIRV-T                 #
#     - code to generate figures 6 to 10 of the thesis                        #
# This script depends on:                                                     #
#     - functions.R                                                           #
############################################################################### 
# hyperparameters 
train = 125:350       # days used for training
test = 350:450        # days used for testing
mse.int = 1:75        # range of days in test interval over which mse is calculated
J = K = 3             # orders of the two filters 

trainV =  333:450     # days used for training of vaccinated data
testV = 450:479       # days used for testing of vaccinated data

beta.seq = seq(0,0.01,0.00002)  # optimisation sequence for Ridge of beta & gamma
gamma.seq = seq(0,0.01,0.00002)

beta.a = 0.12         # guesstimate
w.s = 0.30             # guesstimate - method 1
w.s2 = 0.30            # guesstimate - method 2
sigma.val = 0.2       # vaccine efficiency

# required libraries 
library(ggplot2)
library(deSolve)
library(readr)
library(tidyverse)
library(patchwork)
library(glmnet)
#remotes::install_github('TIvanDijk/TivD')  <-- only used for themes of plots, not necessary but without it layout of plots is different than in thesis 
#library(TivD)

# 1. Load data (generated in preparation.R) & relevant functions ---------------
data <- readRDS('data.RDS')
source('functions.R')

# 2.Parameters ----------------------------------------------------------------
beta.ridge = doRidge(data$beta.t[train], beta.seq, J, length(test))
gamma.ridge = doRidge(data$gamma.t[train], gamma.seq, K, length(test))

betaV = doRidge(data$beta.t[trainV], beta.seq, J, length(testV))
gammaV = doRidge(data$gamma.t[trainV], gamma.seq, K, length(testV))

beta = mean(beta.ridge$param.hat)
gamma = mean(gamma.ridge$param.hat)

n = min(test)    # starting point 
nV = min(testV)  # starting point for vaccinated 

coefs = coef(lm( log(beta.t) ~ I(-1*time), data = data[train,]))
beta.0 = exp(coefs[1]); omega = coefs[2]
beta.s = beta.0 * exp(-omega * test)

coefs = coef(lm( log(alpha.t) ~ I(-1*time), data = data[trainV,]))
alpha.0 = exp(coefs[1]); omega = coefs[2]
alpha = alpha.0 * exp(-omega * testV)

# Verify correctness of the parameters 
betaPlot <- ggplot(data[test,], aes(x = test)) +
  geom_line( aes(y = beta.t, color = 'Actual'), size = 1) +
  geom_line( aes(y = beta.ridge$param.hat, color = 'Predicted (J = 3)'), size = 1, alpha = 0.8) +
  geom_line( aes(y = doRidge(data$beta.t[train], beta.seq, 10, length(test))$param.hat,
                 color = 'Predicted (J = 10)'), size = 1, alpha = 0.8) +
  TivD::theme_newgrey(text = element_text('mono'), legend.position = c(0.35, 0.15)) +
  labs(x = 'Time', y = expression('Transmission Rate'~~beta(t)), color = 'Legend',
       title = expression(bold('A: Accuracy of Transmission Rate'~~beta(t))))

gammaPlot <- ggplot(data[test, ], aes(x = test)) +
  geom_line( aes(y = gamma.t, color = 'Actual'), size = 1) +
  geom_line( aes(y = gamma.ridge$param.hat, color = 'Predicted (K = 3)'), size = 1, alpha = 0.8) +
  geom_line( aes(y = doRidge(data$gamma.t[train], gamma.seq, 10, length(test))$param.hat,
                 color = 'Predicted (K = 10)'), size = 1, alpha = 0.8) +
  TivD::theme_newgrey(text = element_text('mono'), legend.position = c(0.7, 0.15)) +
  labs(x = 'Time', y = expression('Recovery Rate'~~gamma(t)), color = 'Legend',
       title = expression(bold('B: Accuracy of Recovery Rate'~~gamma(t))))

betaPlot + gammaPlot
ggsave('plots/modelParams.pdf', width = 15, height = 6)

ggplot(data[test,], aes(x = test)) +
  geom_line( aes(y = beta.t, color = 'Actual'), size = 1) +
  geom_line( aes(y = beta.s, color = 'Predicted'), size = 1, alpha = 0.8) +
  TivD::theme_newgrey(text = element_text('mono'), legend.position = c(0.25, 0.15)) +
  labs(x = 'Time', y = expression('Transmission Rate'~~beta[s](t)), color = 'Legend')
ggsave('plots/beta_sFIT.pdf', height = 5, width = 8)

ggplot(data[testV, ], aes(x = testV)) +
  geom_line( aes(y = alpha.t * 100, color = 'Actual'), size = 1) +
  geom_line( aes(y = alpha * 100, color = 'Predicted'), size = 1) +
  TivD::theme_newgrey(text = element_text('mono'), legend.position = c(0.85, 0.15)) +
  labs(y = 'Vaccination Rate (in %)', color = 'Legend', x = 'Time')
ggsave('plots/alphaFIT.pdf', height = 5, width = 8)

# 3. Models --------------------------------------------------------------------
# A: Classical
initial_vals = c(S = 17.4e6-data$I[n]-data$R[n], I = data$I[n], R = data$R[n])
parameters = c(gamma = gamma, beta = beta)

# solve the model using deSolve
output <- as.data.frame(ode(y=initial_vals, func = classical.SIR, parms=parameters, times = test)) %>% 
  pivot_longer(2:4, names_to = 'State')

# Plot outcome 
# S, I and R overtime 
ggplot(output, aes(x = time, y = value, color = State, group = State)) +
  geom_line()

# predicted I vs measured cases 
fit.class <- ggplot(subset(output, State == 'I'), aes(x = test, y = value)) +
  geom_line(aes(color = 'Predicted'), size = 1) +
  geom_line(data = data[test,], aes(y = I, color = 'Actual'), size = 1) +
  TivD::theme_newgrey(text = element_text('mono')) +
  labs(title = 'A: Classical SIR', x = 'Time', y = 'Infectious Individuals I(t)',
       color = '')
fit.class

actual = unlist(data[test, 'I'])
predict = unlist(subset(output, State == 'I')[, 'value'])
cat('\nMSE of classical SIR (in millions): ', mse(actual[mse.int], predict[mse.int])/1e6)

# B: time-dependent
initial_vals = c(S = 17.4e6-data$I[n]-data$R[n], I = data$I[n], R = data$R[n])
parameters = c('beta' = beta.ridge$param.hat, 'gamma' = gamma.ridge$param.hat)
output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, parms=parameters, times = 1:length(test[-1]))) %>% 
  pivot_longer(2:4, names_to = 'State')

# Plot outcome 
# S, I and R overtime 
ggplot(output, aes(x = time, y = value, color = State, group = State)) +
  geom_line()

# predicted I vs measured cases 
fit.time <- ggplot(subset(output, State == 'I'), aes(x = test[-1], y = value)) +
  geom_line(aes(color = 'Predicted'), size = 1) +
  geom_line(data = data[test[-length(test)],], aes(y = I, color = 'Actual'), size = 1) +
  TivD::theme_newgrey(text = element_text('mono')) +
  labs(title = paste0('B: Time-Dependent SIR (J = ', J, ', K = ', K, ')'), x = 'Time', y = 'Infectious Individuals I(t)',
       color = '')
fit.time

actual = unlist(data[test[-length(test)], 'I'])
predict = unlist(subset(output, State == 'I')[, 'value'])
cat('\nMSE of SIR-T (in millions): ', mse(actual[mse.int], predict[mse.int])/1e6)

# C: asymptomatic SIR - method 1
parameters = c('beta.s' = beta.s, 'gamma' = gamma.ridge$param.hat, 
               'beta.a' = beta.a, w.s = w.s, w.a = 1-w.s)
initial_vals = c(S = 17.4e6-data$I[n]-data$R[n], Is = data$I[n]/2, 
                 A = data$I[n]/2, R = data$R[n])    

output <- as.data.frame(ode(y=initial_vals, func = asymptomatic.SIR, parms=parameters, times = 1:(length(test[-1])))) %>% 
  mutate(I = Is + A) %>% 
  pivot_longer(2:6, names_to = 'State')

# Plot outcome 
# S, I, A and R overtime 
ggplot(output, aes(x = time, y = value, color = State, group = State)) +
  geom_line() 

title = bquote(bold('C: Asympotomatic SIR I ' ~ (K == .(K)~','~ omega[a] == .(1-w.s) ~','~ beta[a] == .(beta.a))))
# predicted I vs measured cases 
fit.asymptotic <- ggplot(subset(output, State == 'I'), aes(x = test[-1], y = value)) +
  geom_line(aes(color = 'Predicted'), size = 1) +
  geom_line(data = data[test[-length(test)],], aes(y = I, color = 'Actual'), size = 1) +
  TivD::theme_newgrey(text = element_text('mono')) +
  labs(title = title, x = 'Time', y = 'Infectious Individuals I(t)',
       color = '')
fit.asymptotic

actual = unlist(data[test[-length(test)], 'I'])
predict = unlist(subset(output, State == 'I')[, 'value'])
cat('\nMSE of SIAR-T  (I, in millions): ', mse(actual[mse.int], predict[mse.int])/1e6)

# D: asymptomatic SIR - method 2
parameters = c('beta.s' = 1.5*beta.ridge$param.hat, 'gamma' = gamma.ridge$param.hat, 
               'beta.a' = 0.75*beta.ridge$param.hat, w.s = w.s2, w.a = 1-w.s2)
initial_vals = c(S = 17.4e6-data$I[n]-data$R[n], Is = data$I[n]*w.s2, 
                 A = data$I[n]*(1-w.s2), R = data$R[n])    

output <- as.data.frame(ode(y=initial_vals, func = asymptomatic2.SIR, parms=parameters, times = 1:length(test[-1]))) %>% 
  mutate(I = Is + A) %>% 
  pivot_longer(2:6, names_to = 'State')


# Plot outcome 
# S, I and R overtime 
ggplot(output, aes(x = time, y = value, color = State, group = State)) +
  geom_line() 

title = bquote(bold('D: Asympotomatic SIR II ' ~ (J == .(J)~','~K == .(K)~','~ omega[a] == .(1-w.s2))))
# predicted I vs measured cases 
fit.asymptotic2 <- ggplot(subset(output, State == 'I'), aes(x = test[-1], y = value)) +
  geom_line(aes(color = 'Predicted'), size = 1) +
  geom_line(data = data[test[-length(test)],], aes(y = I, color = 'Actual'), size = 1) +
  TivD::theme_newgrey(text = element_text('mono')) +
  labs(title = title, x = 'Time', y = 'Infectious Individuals I(t)',
       color = '')
fit.asymptotic2

actual = unlist(data[test[-length(test)], 'I'])
predict = unlist(subset(output, State == 'I')[, 'value'])
cat('\nMSE of SIAR-T  (II, in millions): ', mse(actual[mse.int], predict[mse.int])/1e6)

# E: vaccinated
initial_vals = c(S = 17.4e6-data$I[nV]-data$R[nV]-data$vacc[nV], I = data$I[nV], R = data$R[nV],
                 V = data$vacc[nV])
parameters = c('beta' = betaV$param.hat, 'gamma' = gammaV$param.hat, 
               'alpha' = alpha, sigma = sigma.val)
output <- as.data.frame(ode(y=initial_vals, func = vaccinated.SIR, parms=parameters, times = 1:length(testV[-1]))) %>% 
  pivot_longer(2:5, names_to = 'State')

# Plot outcome 
# S, I, R and V overtime 
ggplot(output, aes(x = time, y = value, color = State, group = State)) +
  geom_line()

# predicted I vs measured cases 
title = bquote(bold('E: Time-Dependent SIRV ' ~ (J == .(J)~','~K == .(K)~','~ sigma == .(sigma.val))))
fit.vacc <- ggplot(subset(output, State == 'I'), aes(x = testV[-1], y = value)) +
  geom_line(aes(color = 'Predicted'), size = 1) +
  geom_line(data = data[testV[-1],], aes(y = I, color = 'Actual'), size = 1) +
  TivD::theme_newgrey(text = element_text('mono'), legend.position = c(0.85, 0.85)) +
  labs(title = title, x = 'Time', y = 'Infectious Individuals I(t)',
       color = 'Legend')
fit.vacc
ggsave('plots/modelSIRV.pdf', width = 8, height = 5)

# 4. Model Evaluation ----------------------------------------------------------
fit.class + fit.time + fit.asymptotic + fit.asymptotic2 +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
ggsave('plots/modelFIT.pdf', width = 16, height = 10)

