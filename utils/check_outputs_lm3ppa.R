library(readr)
library(ggplot2)
library(dplyr)

##--------------------------------
## Annual outputs
##--------------------------------
df_orig   <- read_csv("~/cagibi/BiomeE-Allocation/model/output/Ecosystem_yearlytest.csv")
#df_noleap <- read_csv("~/cagibi/BiomeE-Allocation/model/output/Ecosystem_yearlytest_NOLEAP.csv")
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
  

<<<<<<< HEAD
#############################
library(ggplot2)

# Files: Ecosystem_yearlytest.csv; Ecosystem_dailytest.csv; Annual_cohortstest.csv, Cohorts_dailytest.csv; PhotosynthesisDynamicstest.csv
outBiomeE <- read.csv("/Users/lmarques/polybox/sofun_plots/outputs_BiomeE/Ecosystem_yearlytest.csv",sep=",")
str(outBiomeE)
dim(outBiomeE)
outSOFUN <- read.csv("/Users/lmarques/polybox/sofun_plots/outputs_SOFUN/Ecosystem_yearly_test.csv",sep=",")
=======




#########3


library(ggplot2)
library(readr)

# Files: Ecosystem_yearlytest.csv; Ecosystem_dailytest.csv; Annual_cohortstest.csv, Cohorts_dailytest.csv; PhotosynthesisDynamicstest.csv
outBiomeE <- read_csv("~/cagibi/BiomeE-Allocation/model/output/Ecosystem_yearlytest.csv")
str(outBiomeE)
dim(outBiomeE)

outSOFUN <- read_csv("/Users/bestocke/sofun/output/Ecosystem_yearly_test.csv")
>>>>>>> 0b032c4c78c37b7a4a2556782fa8c9d92abe224d
str(outSOFUN)
dim(outSOFUN)

# Subset data Ecosystem_daily_test in outputs_SOFUN from year = 1801
outSOFUN <- subset(outSOFUN, year>= 1801)
<<<<<<< HEAD

outSOFUN <- subset(outSOFUN, PFT== 1)
outBiomeE <- subset(outBiomeE, PFT== 1)


ggplot() +
  geom_line(data=outBiomeE, aes(x=year, y=Evap), color='blue') +
  geom_line(data=outSOFUN, aes(x=year, y=Evap), color='red')

ggplot() +
  geom_line(data=outBiomeE, aes(x=year, y=leafN), color='blue') +
=======
outSOFUN <- subset(outSOFUN, PFT== 1)
outBiomeE <- subset(outBiomeE, PFT== 1)

ggplot() +
  geom_line(data=outBiomeE, aes(x=year, y=GPP), color='blue', size = 2) +
  geom_line(data=outSOFUN, aes(x=year, y=GPP), color='red') + 
  xlim( c(0,200) )

ggplot() +
  geom_line(data=outBiomeE, aes(x=year, y=leafN), color='blue', size = 2) +
>>>>>>> 0b032c4c78c37b7a4a2556782fa8c9d92abe224d
  geom_line(data=outSOFUN, aes(x=year, y=leafN), color='red')

## Loop for columns in two data.frames
plot_dfs <- lapply(names(outBiomeE),function(nm)data.frame(col1 = outBiomeE[,1],col2 = outBiomeE[,nm], col3 = outSOFUN[,nm]))
<<<<<<< HEAD

=======
>>>>>>> 0b032c4c78c37b7a4a2556782fa8c9d92abe224d
for (idx in seq_along(plot_dfs)){
  tiff(file = paste(names(outBiomeE)[idx], '.tiff', sep = ''))
  print(
    ggplot()+geom_smooth(data = plot_dfs[[idx]], aes(x=col1, y=col2),color="blue") +
      geom_smooth(data = plot_dfs[[idx]], aes(x=col1, y=col3),color="red")+
      ggtitle(names(outBiomeE)[idx]))
<<<<<<< HEAD
  dev.off()
}

=======
  #ggsave("path.tiff")
  dev.off()
}
>>>>>>> 0b032c4c78c37b7a4a2556782fa8c9d92abe224d
#for (df in plot_dfs){
# print(
#  ggplot()+geom_line(data = df, aes(x=col1, y=col2),color="blue") +
#           geom_line(data = df, aes(x=col1, y=col3),color="red"))}
<<<<<<< HEAD



=======
>>>>>>> 0b032c4c78c37b7a4a2556782fa8c9d92abe224d
