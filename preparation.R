###############################################################################
#                             DATA PREPARATION                                #
# This script contains the following:                                         #
#     - creation of 'data.RDS' based on data from the RIVM                    #
#     - visualisations of the dataset                                         #
############################################################################### 
library(readr)
library(ggplot2)
library(tidyverse)
library(patchwork)

# dataset that contains reported cases
ma <- function(x, n = 7){stats::filter(c(rep(0,n), x), rep(1 / n, n), sides = 1)[-c(1:n)]}

data.cases <- read_delim('https://data.rivm.nl/covid-19/COVID-19_aantallen_gemeente_per_dag.csv', delim = ';') %>% 
  group_by(Date_of_publication) %>% 
  summarise(cases = sum(Total_reported)) %>% 
  dplyr::transmute(date = Date_of_publication, cases.reported = cases,
                   cases = round(ma(cases), 0), time = 1:n(),
                   I = cumsum(cases) - dplyr::lag(cumsum(cases), 10, default = 0),  # assuming individuals remain infectious for 7 days
                   R = dplyr::lag(cumsum(cases), 10, default = 0))  # assuming individuals recover after 7 days

# dataset (NICE) that contains hospital admissions 
data.hosp <- read_delim('https://data.rivm.nl/covid-19/COVID-19_ziekenhuisopnames.csv', delim = ';') %>% 
  group_by(Date_of_statistics) %>% 
  summarise(hosp = sum(Hospital_admission)) %>% 
  dplyr::transmute(date = Date_of_statistics, hosp = hosp)

# combine the datasets 
data <- left_join(data.hosp, data.cases)
data <- data[-1, ] %>%  #remove first observations
  mutate(time = time - 1)

data$beta.t = ma((dplyr::lead(data$I) - data$I + dplyr::lead(data$R) - data$R) / data$I) # eq. 11 of thesis
data$gamma.t  = ma((dplyr::lead(data$R) - data$R) / data$I)   # eq. 9 of thesis

# store data
saveRDS(data, 'data.RDS')

# data description 
#remotes::install_github('TIvanDijk/TivD')  <-- only used for themes of plots, not necessary but without it layout of plots is different than in thesis 
#library(TivD)
p1 <- ggplot(data, aes(x = date, y = cases.reported)) +
  annotate('rect', xmin = data$date[125], xmax = data$date[350], ymin = 0, ymax = 13000, 
           alpha = 0.3, fill = 'steelblue') +
  annotate('text', x = data$date[125], y = 12700, label = 'training', hjust = 0, 
           fontface = 'bold', color = 'steelblue', family = 'mono') +
  annotate('rect', xmin = data$date[350], xmax = data$date[450], ymin = 0, ymax = 13000, 
           alpha = 0.3, fill = 'red') +
  annotate('text', x = data$date[350], y = 12700, label = 'testing', hjust = 0, 
           fontface = 'bold', color = 'red', family = 'mono') +
  geom_col(fill = 'mistyrose3', color = 'mistyrose3') +
  scale_x_date(date_breaks = '5 months', date_labels="%B %d,%Y") +
  labs(x = 'Date', y = 'Reported Cases',
       title = 'A: Cases') +
  TivD::theme_newgrey(text = element_text('mono'))

p2 <- ggplot(data, aes(x = date, y = hosp)) +
  annotate('rect', xmin = data$date[125], xmax = data$date[350], ymin = 0, ymax = 650, 
           alpha = 0.3, fill = 'steelblue') +
  annotate('text', x = data$date[125], y = 620, label = 'training', hjust = 0, 
           fontface = 'bold', color = 'steelblue', family = 'mono') +
  annotate('rect', xmin = data$date[350], xmax = data$date[450], ymin = 0, ymax = 650, 
           alpha = 0.3, fill = 'red') +
  annotate('text', x = data$date[350], y = 620, label = 'testing', hjust = 0, 
           fontface = 'bold', color = 'red', family = 'mono') +
  geom_col(fill = 'khaki3', color = 'khaki3') +
  scale_x_date(date_breaks = '5 months', date_labels="%B %d,%Y") +
  labs(x = 'Date', y = 'Reported Hospitalisations', 
       title = 'B: Hospitalisations' ) +
  TivD::theme_newgrey(text = element_text('mono'))

p1 + p2
ggsave('plots/reportedCases.pdf', width = 16, height = 4.5)

# some 'summary' statistics used in the paper
cat('The dataset ranges from', format(min(data$date), '%B %d, %Y'), 'to', 
    format(max(data$date), '%B %d, %Y'))


