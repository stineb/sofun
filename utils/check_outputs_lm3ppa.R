library(readr)
library(ggplot2)
library(dplyr)

##--------------------------------
## Comparing Original BiomeE-Allocation with implementation in SOFUN
##--------------------------------
## Annual outputs
df_orig   <- read_csv("~/cagibi/BiomeE-Allocation/model/output/Ecosystem_yearlytest.csv")
df_sofun  <- read_csv("~/sofun/output/Ecosystem_yearly_test.csv")

df_test <- select(df_orig, year, LAI_orig = LAI) %>% 
  left_join(select(df_sofun, year, LAI_sofun = LAI), by = "year")

## slight differences for soilN
## sudden difference for SapwoodC in some years
ggplot() +
  geom_line(aes(x = year, y = NSN), data = df_orig, color = "black", size = 2) +
  geom_line(aes(x = year, y = NSN), data = df_sofun, color = "red", size = 1) +
  xlim(c(0,200))

ggplot() +
  geom_line(aes(x = year, y = NSN), data = df_orig, color = "black", size = 2) +
  geom_line(aes(x = year, y = NSN), data = df_sofun, color = "red", size = 1) +
  xlim(c(1600,1809))

ggplot() +
  geom_line(aes(x = year, y = mineralN), data = df_sofun, color="red") +
  geom_line(aes(x = year, y = mineralN), data = df_orig) + 
  ylim( c(0,0.3) )

ggplot() +
  geom_point(aes(x = LAI_orig, y = LAI_sofun), data = df_test) +
  geom_abline(intercept=0, slope=1, linetype="dotted")
  

##--------------------------------
## Comparing latest development with reference
##--------------------------------
## Annual outputs
df     <- read_csv("~/sofun/output/Annual_tile_test.csv")
df_ref <- read_csv("~/sofun/output/lm3ppa_pmodel_runlm3ppa_pmodel/Annual_tile_test.csv")

df_test <- select(df_ref, year, LAI_ref = LAI) %>% 
  left_join(select(df, year, LAI = LAI), by = "year")

ggplot() +
  geom_point(aes(x = LAI_ref, y = LAI), data = df_test) +
  geom_abline(intercept=0, slope=1, linetype="dotted")


##--------------------------------
## Comparing SOFUN branches lm3ppa_pmodel with lm3ppa
##--------------------------------
## Annual outputs
df_lm3ppa        <- read_csv("~/sofun/output/lm3ppa_runlm3ppa/Annual_tile_test.csv")
df_lm3ppa_pmodel <- read_csv("~/sofun/output/lm3ppa_pmodel_runlm3ppa/Annual_tile_test.csv")

df_test <- select(df_lm3ppa, year, LAI_lm3ppa = LAI) %>% 
  left_join(select(df_lm3ppa_pmodel, year, LAI_lm3ppa_pmodel = LAI), by = "year")

ggplot() +
  geom_point(aes(x = LAI_lm3ppa, y = LAI_lm3ppa_pmodel), data = df_test) +
  geom_abline(intercept=0, slope=1, linetype="dotted")

## ok, identical

##--------------------------------
## Comparing LM3-PPA with different photosynthesis: gs_Leuning() vs. P-model
##--------------------------------
## Annual outputs
df_lm3ppa        <- read_csv("~/sofun/output/lm3ppa_runlm3ppa/Annual_tile_test.csv")
df_lm3ppa_pmodel <- read_csv("~/sofun/output/lm3ppa_runlm3ppa_pmodel/Annual_tile_test.csv")

df_test <- select(df_lm3ppa, year, GPP_lm3ppa = GPP) %>% 
  left_join(select(df_lm3ppa_pmodel, year, GPP_lm3ppa_pmodel = GPP), by = "year")

ggplot() +
  geom_line(aes(x = year, y = GPP_lm3ppa), data = df_test) +
  geom_line(aes(x = year, y = GPP_lm3ppa_pmodel), data = df_test, color = 'royalblue')


##--------------------------------
## Checking daily and hourly consistency (using hack: constant LUE)
##--------------------------------
## Annual outputs
df_d <- read_csv("~/sofun/output/lm3ppa_hackdaily/Annual_tile_test.csv")
df_h <- read_csv("~/sofun/output/lm3ppa_hackhourly/Annual_tile_test.csv")

df_test <- select(df_d, year, GPP_d = GPP) %>% 
  left_join(select(df_h, year, GPP_h = GPP), by = "year")

ggplot() +
  geom_point(aes(x = GPP_h, y = GPP_d), data = df_test) +
  geom_abline(intercept=0, slope=1, linetype="dotted")


##--------------------------------
## Checking SOFUN implementation of BiomeE-Allocation agains P-model in LM3-PPA
##--------------------------------
## Annual outputs
df_orig <- read_csv("~/sofun/output/lm3ppa_runlm3ppa/Annual_tile_test.csv")
df_lm3ppa_pmodel <- read_csv("~/sofun/output/lm3ppa_lm3ppa_pmodel/Annual_tile_test.csv")

df_test <- select(df_orig, year, GPP_orig = GPP) %>% 
  left_join(select(df_lm3ppa_pmodel, year, GPP_lm3ppa_pmodel = GPP), by = "year")

ggplot() +
  geom_line(aes(x = year, y = GPP_orig), data = df_test) +
  geom_line(aes(x = year, y = GPP_lm3ppa_pmodel), data = df_test, color = 'royalblue')
  

##--------------------------------
## Checking daily hack consistency between original and pmodel implementation
##--------------------------------
## Annual outputs
df_d <- read_csv("~/sofun/output/lm3ppa_hackdaily/Annual_tile_test.csv")
df_h <- read_csv("~/sofun/output/lm3ppa_hackdaily_pmodel/Annual_tile_test.csv")

df_test <- select(df_d, year, GPP_d = GPP) %>% 
  left_join(select(df_h, year, GPP_h = GPP), by = "year")

ggplot() +
  geom_point(aes(x = GPP_h, y = GPP_d), data = df_test) +
  geom_abline(intercept=0, slope=1, linetype="dotted")

ggplot() +
  geom_line(aes(x = year, y = GPP_d), data = df_test) +
  geom_line(aes(x = year, y = GPP_h), data = df_test, color = 'royalblue')
