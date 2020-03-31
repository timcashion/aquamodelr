

library(aquamodelr)
library(tidyverse)
#clean_fao_aqua(path = "./raw_data/FAO_Global_Production_2019.1.0.zip", output_path="./data/FAO_Global_Production_2019.1.0.csv") #Only need to run the first time
fao_aqua <- read_csv("../FishmealDemand/data/FAO_Global_Production_2019.1.0.csv")
fm_data <- readr::read_csv("../FishmealDemand/data/FishmealModeling_Oct062018.csv")
fm_distributions <- estimate_fmfo_distributions(data= fm_data)
