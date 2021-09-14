# Minimum replication package
rm(list = ls())

# Load libraries
library(tidyverse)
library(plm)
library(splm)
library(spdep)
library(sf)
library(lmtest)
library(spatialreg)
library(broom)
library(lfe)
library(Formula)

#### Data ####
load("data_rp.RData")

# Table: india_data
# Variable:
# "id"...district id
# "time"...year
# "Lights"...DMSP (t)
# "Lights_L1"...DMSP (t-1)
# "PVpr"...P. vivax parasite rate (t)
# "PVpr_L1"...P. vivax parasite rate (t-1)
# "PFpr"...P. falciparum parasite rate (t)
# "PFpr_L1"...P. falciparum parasite rate (t-1)
# "PVpr_L2"...P. vivax parasite rate (t-2)
# "PFpr_L2"...P. falciparum parasite rate (t-2)

# Variable ww: Object with "matrix" W from Eq. (2) and (3)

#### Definition of models (class formula) ####

model <- I(Lights-Lights_L1) ~ Lights_L1

modell <- list(
  update(model, . ~ . + I(PFpr - PFpr_L1)),
  update(model, . ~ . + I(PVpr - PVpr_L1))
)


modell_out <- list(
  update(model, . ~ . + I((PFpr - PFpr_L1)>=0.202 & (PFpr - PFpr_L1)<0.547) + I((PFpr - PFpr_L1)>=0.547 & (PFpr - PFpr_L1)<1.66) + I((PFpr - PFpr_L1)>=1.66)),
  update(model, . ~ . + I((PVpr - PVpr_L1)>=0.126 & (PVpr - PVpr_L1)<0.420) + I((PVpr - PVpr_L1)>=0.420 & (PVpr - PVpr_L1)<1.5) + I((PVpr - PVpr_L1)>=1.5))
)

modell_pfpv <- list(
  update(model, . ~ . + I(PFpr - PFpr_L1) + I(PVpr - PVpr_L1))
)

modell_out_pfpv <- list(
  update(model, . ~ . + I((PFpr - PFpr_L1)>=0.202 & (PFpr - PFpr_L1)<0.547) + I((PFpr - PFpr_L1)>=0.547 & (PFpr - PFpr_L1)<1.66) + I((PFpr - PFpr_L1)>=1.66) + 
           I((PVpr - PVpr_L1)>=0.126 & (PVpr - PVpr_L1)<0.420) + I((PVpr - PVpr_L1)>=0.420 & (PVpr - PVpr_L1)<1.5) + I((PVpr - PVpr_L1)>=1.5))
)

modell_lag2 <- list(
  update(model, . ~ . + I(PFpr - PFpr_L1) + I(PFpr_L1 - PFpr_L2)),
  update(model, . ~ . + I(PVpr - PVpr_L1) + I(PVpr_L1 - PVpr_L2))
)

modell_out_lag2 <- list(
  update(model, . ~ . + I((PFpr - PFpr_L1)>=0.202 & (PFpr - PFpr_L1)<0.547) + I((PFpr - PFpr_L1)>=0.547 & (PFpr - PFpr_L1)<1.66) + I((PFpr - PFpr_L1)>=1.66) + 
           I((PFpr_L1 - PFpr_L2)>=0.202 & (PFpr_L1 - PFpr_L2)<0.547) + I((PFpr_L1 - PFpr_L2)>=0.547 & (PFpr_L1 - PFpr_L2)<1.66) + I((PFpr_L1 - PFpr_L2)>=1.66)),
  update(model, . ~ . + I((PVpr - PVpr_L1)>=0.126 & (PVpr - PVpr_L1)<0.420) + I((PVpr - PVpr_L1)>=0.420 & (PVpr - PVpr_L1)<1.5) + I((PVpr - PVpr_L1)>=1.5) + 
           I((PVpr_L1 - PVpr_L2)>=0.126 & (PVpr_L1 - PVpr_L2)<0.420) + I((PVpr_L1 - PVpr_L2)>=0.420 & (PVpr_L1 - PVpr_L2)<1.5) + I((PVpr_L1 - PVpr_L2)>=1.5))
)

#### Results reported in Table 1 ####

# Fixed effects

modell %>% 
  map(as.Formula) %>% 
  map(
    update, . ~ . | id + time | 0 | id
  ) %>% 
  map(
    felm, data = india_data
  ) %>% 
  map(
    summary
  )

# SARAR

modell %>%
  map(
    spml, data = india_data, listw = ww, index = c("id","time"), lag = TRUE, effect = "twoways", spatial.error="kkp"
  ) %>% 
  map(
    summary
  )

# Test : Spatial lag

test_lml <- modell %>% 
  map(
    slmtest, data = india_data, index = c('id', 'time'), listw = ww, test = 'lml'
  )

print(test_lml)

# Test : Spatial error
test_lme <- modell %>% 
  map(
    slmtest, data = india_data, index = c('id', 'time'), listw = ww, test = 'lme'
  )

print(test_lme)


#### Results reported in Table 2 ####

# Fixed effects

modell_out %>% 
  map(as.Formula) %>% 
  map(
    update, . ~ . | id + time | 0 | id
  ) %>% 
  map(
    felm, data = india_data
  ) %>% 
  map(
    summary
  )

# SARAR

modell_out %>%
  map(
    spml, data = india_data, listw = ww, index = c("id","time"), lag = TRUE, effect = "twoways", spatial.error="kkp"
  ) %>% 
  map(
    summary
  )

# Test : Spatial lag

test_lml <- modell_out %>% 
  map(
    slmtest, data = india_data, index = c('id', 'time'), listw = ww, test = 'lml'
  )

print(test_lml)

# Test : Spatial error
test_lme <- modell_out %>% 
  map(
    slmtest, data = india_data, index = c('id', 'time'), listw = ww, test = 'lme'
  )

print(test_lme)


#### Results reported in Table 3 ####

# Fixed effects

c(modell_pfpv,
  modell_out_pfpv) %>% 
  map(as.Formula) %>% 
  map(
    update, . ~ . | id + time | 0 | id
  ) %>% 
  map(
    felm, data = india_data
  ) %>% 
  map(
    summary
  )

# SARAR

c(modell_pfpv,
  modell_out_pfpv) %>%
  map(
    spml, data = india_data, listw = ww, index = c("id","time"), lag = TRUE, effect = "twoways", spatial.error="kkp"
  ) %>% 
  map(
    summary
  )

#### Results reported in Table A.4 ####

# Remove missing observations (b/c SARAR)
india_data_l2 <- india_data %>% drop_na()

# Fixed effects

modell_lag2 %>% 
  map(as.Formula) %>% 
  map(
    update, . ~ . | id + time | 0 | id
  ) %>% 
  map(
    felm, data = india_data_l2
  ) %>% 
  map(
    summary
  )

# SARAR

modell_lag2 %>%
  map(
    spml, data = india_data_l2, listw = ww, index = c("id","time"), lag = TRUE, effect = "twoways", spatial.error="kkp"
  ) %>% 
  map(
    summary
  )

#### Results reported in Table A.5 ####

# Fixed effects

c(modell_lag2,
  modell_out_lag2) %>% 
  map(as.Formula) %>% 
  map(
    update, . ~ . | id + time | 0 | id
  ) %>% 
  map(
    felm, data = india_data_l2
  ) %>% 
  map(
    summary
  )

# SARAR

c(modell_lag2,
  modell_out_lag2) %>%
  map(
    spml, data = india_data_l2, listw = ww, index = c("id","time"), lag = TRUE, effect = "twoways", spatial.error="kkp"
  ) %>% 
  map(
    summary
  )

