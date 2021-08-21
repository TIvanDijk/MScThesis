###############################################################################
#                             POLICY EVALUATION                               #
# This script contains the following:                                         #
#     - implementation of ch 4.5                                              #
#     - code to generate figure 11                                            #
# This script depends on:                                                     #
#     - functions.R                                                           #
############################################################################### 
# 0. hyperparameters
train_size = 50
test_size = 25

J = K = 3

sigma.val = 0.20 # based on combined efficiency of different vaccines 

# used OMT reports are available at: https://www.rijksoverheid.nl/documenten
exp.curfew = 0.9 # reduction between 8% and 13% [OMT, 96] --> say 10%
des.curfew = 0.7 # WANTED: steady decline

exp.nocurfew = 1.10 # OMT 111 restates the 10% claim 
des.nocurfew = 0.85 # WANTED: vaccination-effect to continue

exp.norest = 0.9 # OMT 80: 10% for full package --> full closure more strict
des.norest = 0.70 # WANTED: 'break' exponential increase 

exp.mask = 1.05 # OMT 80: beperkt effect --> 0
des.mask = 0.8

# 1. Load data (generated in preparation.R) & relevant functions ---------------
data <- readRDS('data.RDS')
source('functions.R')

library(tidyverse)
library(patchwork)
library(deSolve)
library(glmnet)
#library(TivD)

# --- 2A --- 
# introduction of curfew (23 jan, https://www.rijksoverheid.nl/actueel/nieuws/2021/01/21/vanaf-zaterdag-23-januari-geldt-een-avondklok)
doi <- as.Date('23-01-2021', '%d-%m-%Y')
poi <- which(data$date == doi)
data.curfew = data[(poi - train_size):(poi + test_size-1), ]

# parameters 
initial_vals = c(S = 17.4e6-data$I[poi]-data$R[poi], I = data$I[poi], R = data$R[poi])

beta.ridge = doRidge(data.curfew$beta.t[1:train_size], seq(0,0.01,0.00002), J, test_size+1)
gamma.ridge = doRidge(data.curfew$gamma.t[1:train_size], seq(0,0.01,0.00002), K, test_size+1)

# benchmark case 
bm.parameters = c('beta' = beta.ridge$param.hat, 'gamma' = gamma.ridge$param.hat)
bm.output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, 
                               parms=bm.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))

# expected case 
exp.parameters = c('beta' = exp.curfew * beta.ridge$param.hat, 
                   'gamma' = gamma.ridge$param.hat)
exp.output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, 
                               parms=exp.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))

# desired case 
des.parameters = c('beta' = des.curfew * beta.ridge$param.hat, 
                   'gamma' = gamma.ridge$param.hat)
des.output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, 
                                parms=des.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))
  
p1 <- ggplot(data.curfew, aes(x = date, y = I)) +
  annotate('rect', xmin = data$date[poi-train_size]-.5, xmax = doi, ymin = 0, ymax = 120000, 
           alpha = 0.3, fill = 'steelblue') +
  annotate('text', x = data$date[poi-train_size], y = 116000, label = 'training', hjust = 0, 
           fontface = 'bold', color = 'steelblue', family = 'mono', size = 6) +
  geom_col(aes(fill = 'observed'), color = 'grey40', alpha = 0.8) +
  geom_line(data = bm.output, aes(y = value, color = 'benchmark'), size = 1.2) +
  geom_line(data = exp.output, aes(y = value, color = 'expected'), size = 1.2) +
  geom_line(data = des.output, aes(y = value, color = 'desired'), size = 1.2) +
  geom_vline(xintercept = doi, color = 'red', size = 1.2) +
  scale_x_date(date_breaks = '26 days', date_labels="%B %d, %Y", expand = c(0,0.5)) +
  scale_y_continuous(expand = c(0,0.5), limits = c(0, 125000)) +
  scale_fill_manual(values = 'grey70', name = ' ') +
  scale_color_manual(values = c('darkred', 'darkgreen', 'darkblue'), name = 'prediction type: ') +
  labs(x = '', y = 'Infectious Individuals I(t)', 
       title = 'A: Introduction of Curfew (January 23, 2021)') +
  theme_newgrey(text = element_text('mono')) 


# --- 2B --- 
#Lifting of curfew (28 Apr, https://www.rijksoverheid.nl/onderwerpen/coronavirus-covid-19/avondklok)
doi <- as.Date('28-04-2021', '%d-%m-%Y')
poi <- which(data$date == doi)

data.nocurfew = data[(poi - train_size):(poi + test_size-1), ]

# parameters 
initial_vals = c(S = 17.4e6-data$I[poi]-data$R[poi], I = data$I[poi], 
                 R = data$R[poi], V = data$vacc[poi])

coefs = coef(lm( log(alpha.t) ~ I(-1*time), data = data.nocurfew[1:train_size,]))
alpha.0 = exp(coefs[1]); omega = coefs[2]
time_vec = data.nocurfew$time[train_size] : (data.nocurfew$time[train_size] + test_size + 1)
alpha = alpha.0 * exp(-omega * time_vec)

beta.ridge = doRidge(data.nocurfew$beta.t[1:train_size], seq(0,0.01,0.00002), J, test_size+1)
gamma.ridge = doRidge(data.nocurfew$gamma.t[1:train_size], seq(0,0.01,0.00002), K, test_size+1)

# benchmark case 
bm.parameters = c('beta' = beta.ridge$param.hat, 'gamma' = gamma.ridge$param.hat, 
                  'alpha' = as.numeric(alpha), sigma = sigma.val)
bm.output <- as.data.frame(ode(y=initial_vals, func = vaccinated.SIR, 
                               parms=bm.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))

# expected case 
exp.parameters = c('beta' = exp.nocurfew * beta.ridge$param.hat, 
                   'gamma' = gamma.ridge$param.hat, 
                   'alpha' = as.numeric(alpha), sigma = sigma.val)
exp.output <- as.data.frame(ode(y=initial_vals, func = vaccinated.SIR, 
                                parms=exp.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))

# desired case 
des.parameters = c('beta' = des.nocurfew * beta.ridge$param.hat, 
                   'gamma' = gamma.ridge$param.hat, 
                   'alpha' = as.numeric(alpha), sigma = sigma.val)
des.output <- as.data.frame(ode(y=initial_vals, func = vaccinated.SIR, 
                                parms=des.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))


p2 <- ggplot(data.nocurfew, aes(x = date, y = I)) +
  annotate('rect', xmin = data$date[poi-train_size]-.5, xmax = doi, ymin = 0, ymax = 120000, 
           alpha = 0.3, fill = 'steelblue') +
  annotate('text', x = data$date[poi-train_size], y = 116000, label = 'training', hjust = 0, 
           fontface = 'bold', color = 'steelblue', family = 'mono', size = 6) +
  geom_col(aes(fill = 'observed'), color = 'grey40', alpha = 0.8) +
  geom_line(data = bm.output, aes(y = value, color = 'benchmark'), size = 1.2) +
  geom_line(data = exp.output, aes(y = value, color = 'expected'), size = 1.2) +
  geom_line(data = des.output, aes(y = value, color = 'desired'), size = 1.2) +
  geom_vline(xintercept = doi, color = 'red', size = 1.2) +
  scale_x_date(date_breaks = '25 days', date_labels="%B %d, %Y", expand = c(0,0.5)) +
  scale_y_continuous(expand = c(0,0.5), limits = c(0, 125000)) +
  scale_fill_manual(values = 'grey70', name = ' ') +
  scale_color_manual(values = c('darkred', 'darkgreen', 'darkblue'), name = 'prediction type: ') +
  labs(x = '', y = '', 
       title = 'B: Lifting of Curfew (April 28, 2021)') +
  theme_newgrey(text = element_text('mono'), axis.text.y = element_blank(), 
                axis.ticks.y  = element_blank()) 

# --- 2C --- 
# Closure of restaurants (14 Oct, https://www.rijksoverheid.nl/onderwerpen/coronavirus-tijdlijn/oktober-2020-tweede-golf-en-gedeeltelijke-lockdown)
doi <- as.Date('14-10-2020', '%d-%m-%Y')
poi <- which(data$date == doi)

data.norestaurant = data[(poi - train_size):(poi + test_size - 1), ]

# parameters 
initial_vals = c(S = 17.4e6-data$I[poi]-data$R[poi], I = data$I[poi], R = data$R[poi])

beta.ridge = doRidge(data.norestaurant$beta.t[1:train_size], seq(0,0.01,0.00002), J, test_size+1)
gamma.ridge = doRidge(data.norestaurant$gamma.t[1:train_size], seq(0,0.01,0.00002), K, test_size+1)

# benchmark case 
bm.parameters = c('beta' = beta.ridge$param.hat, 'gamma' = gamma.ridge$param.hat)
bm.output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, 
                               parms=bm.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))

# expected case 
exp.parameters = c('beta' = exp.norest * beta.ridge$param.hat, 
                   'gamma' = gamma.ridge$param.hat)
exp.output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, 
                                parms=exp.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))

# desired case 
des.parameters = c('beta' = des.norest * beta.ridge$param.hat, 
                   'gamma' = gamma.ridge$param.hat)
des.output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, 
                                parms=des.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))

p3 <- ggplot(data.norestaurant, aes(x = date, y = I)) +
  annotate('rect', xmin = data$date[poi-train_size]-.5, xmax = doi, ymin = 0, ymax = 120000, 
           alpha = 0.3, fill = 'steelblue') +
  annotate('text', x = data$date[poi-train_size], y = 116000, label = 'training', hjust = 0, 
           fontface = 'bold', color = 'steelblue', family = 'mono', size = 6) +
  geom_col(aes(fill = 'observed'), color = 'grey40', alpha = 0.8) +
  geom_line(data = bm.output, aes(y = value, color = 'benchmark'), size = 1.2) +
  geom_line(data = exp.output, aes(y = value, color = 'expected'), size = 1.2) +
  geom_line(data = des.output, aes(y = value, color = 'desired'), size = 1.2) +
  geom_vline(xintercept = doi, color = 'red', size = 1.2) +
  scale_x_date(date_breaks = '25 days', date_labels="%B %d, %Y", expand = c(0,0.5)) +
  scale_y_continuous(expand = c(0,0.5), limits = c(0, 125000)) +
  scale_fill_manual(values = 'grey70', name = ' ') +
  scale_color_manual(values = c('darkred', 'darkgreen', 'darkblue'), name = 'prediction type: ') +
  labs(x = 'Date', y = 'Infectious Individuals I(t)', 
       title = 'C: Closure of Restaurants (October 14, 2020)') +
  theme_newgrey(text = element_text('mono')) 

# --- 2D ---  
# Mandatory face masks (1 Dec, https://www.rijksoverheid.nl/onderwerpen/coronavirus-tijdlijn/december-2020-lockdown-tijdens-feestdagen-en-mutatie-van-virus-duikt-op-in-verenigd-koninkrijk)
doi <- as.Date('01-12-2020', '%d-%m-%Y')
poi <- which(data$date == doi)

data.facemask = data[(poi - train_size):(poi + test_size - 1), ]

# parameters 
initial_vals = c(S = 17.4e6-data$I[poi]-data$R[poi], I = data$I[poi], R = data$R[poi])

beta.ridge = doRidge(data.facemask$beta.t[1:train_size], seq(0,0.01,0.00002), J, test_size+1)
gamma.ridge = doRidge(data.facemask$gamma.t[1:train_size], seq(0,0.01,0.00002), K, test_size+1)

# benchmark case 
bm.parameters = c('beta' = beta.ridge$param.hat, 'gamma' = gamma.ridge$param.hat)
bm.output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, 
                               parms=bm.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))

# expected case 
exp.parameters = c('beta' = exp.mask * beta.ridge$param.hat, 
                   'gamma' = gamma.ridge$param.hat)
exp.output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, 
                                parms=exp.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))

# desired case 
des.parameters = c('beta' = des.mask * beta.ridge$param.hat, 
                   'gamma' = gamma.ridge$param.hat)
des.output <- as.data.frame(ode(y=initial_vals, func = timedependent.SIR, 
                                parms=des.parameters, times = 1:test_size)) %>% 
  pivot_longer(2:4, names_to = 'State') %>% 
  dplyr::filter(State == 'I') %>%  # store only relevant numbers for plot
  mutate(date = seq(doi, doi + test_size - 1, by = 'days'))


p4 <- ggplot(data.facemask, aes(x = date, y = I)) +
  annotate('rect', xmin = data$date[poi-train_size]-.5, xmax = doi, ymin = 0, ymax = 120000, 
           alpha = 0.3, fill = 'steelblue') +
  annotate('text', x = data$date[poi-train_size], y = 116000, label = 'training', hjust = 0, 
           fontface = 'bold', color = 'steelblue', family = 'mono', size = 6) +
  geom_col(aes(fill = 'observed'), color = 'grey40', alpha = 0.8) +
  geom_line(data = bm.output, aes(y = value, color = 'benchmark'), size = 1.2) +
  geom_line(data = exp.output, aes(y = value, color = 'expected'), size = 1.2) +
  geom_line(data = des.output, aes(y = value, color = 'desired'), size = 1.2) +
  geom_vline(xintercept = doi, color = 'red', size = 1.2) +
  scale_x_date(date_breaks = '25 days', date_labels="%B %d, %Y", expand = c(0,0.5)) +
  scale_y_continuous(expand = c(0,0.5), limits = c(0, 125000)) +
  scale_fill_manual(values = 'grey70', name = ' ') +
  scale_color_manual(values = c('darkred', 'darkgreen', 'darkblue'), name = 'prediction type: ') +
  labs(x = 'Date', y = '', 
       title = 'D: Mandatory Face Masks (December 1, 2020)') +
  theme_newgrey(text = element_text('mono'), axis.text.y = element_blank(), 
                axis.ticks.y  = element_blank()) 

(p1 + p2) / (p3 + p4) +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom', legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 12))
ggsave('plots/policies.pdf', width = 19, height = 8.5)

