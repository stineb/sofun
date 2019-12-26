df <- read_delim(file = "~/cagibi/ForestESS/input/Temperate_forcing.txt", delim = "\t")

## write annual files for sofun
write_oneyear <- function(year, df){
  df %>% 
    mutate(date = ymd(paste0(as.character(YEAR), "-01-01")) + days(DOY-1)) %>% 
    mutate(moy = month(date), dom = mday(date)) %>% 
    filter(!(moy == 2 & dom == 29)) %>% 
    dplyr::filter(YEAR==year) %>% 
    write_delim(path = paste0("~/sofun/input/Temperate_forcing_", year, ".txt"), delim = "\t")
}

purrr::map(
  as.list(df$YEAR %>% unique()),
  ~write_oneyear(., df)
)

## write a file for testing where Feb 29 are removed
library(lubridate)
df_test <- df %>% 
  mutate(date = ymd(paste0(as.character(YEAR), "-01-01")) + days(DOY-1)) %>% 
  mutate(moy = month(date), dom = mday(date)) %>% 
  filter(!(moy == 2 & dom == 29))

df_test %>% write_delim(path = paste0("~/cagibi/ForestESS/input/Temperate_forcing_TEST.txt"), delim = "\t")
