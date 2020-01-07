library(readr)
library(ggplot2)

df_orig   <- read_csv("~/cagibi/BiomeE-Allocation/model/output/Ecosystem_yearlytest.csv")
df_noleap <- read_csv("~/cagibi/BiomeE-Allocation/model/output/Ecosystem_yearlytest_NOLEAP.csv")
df_sofun  <- read_csv("~/sofun/output/Ecosystem_yearlytest.csv")

df_test <- select(df_orig, year, LAI_orig = LAI) %>% 
  left_join(select(df_sofun, year, LAI_sofun = LAI), by = "year")

ggplot() +
  geom_line(aes(x = year, y = LAI), data = df_sofun, color="red") +
  geom_line(aes(x = year, y = LAI), data = df_orig)

ggplot() +
  geom_line(aes(x = year, y = mineralN), data = df_sofun, color="red") +
  geom_line(aes(x = year, y = mineralN), data = df_orig) + 
  ylim( c(0,0.3) )

ggplot() +
  geom_point(aes(x = LAI_orig, y = LAI_sofun), data = df_test)
  