

library(aquamodelr)
library(tidyverse)

if(file.exists("./data/FAO_Global_Production_2019.1.0.csv")==F){
  clean_fao_aqua(path_to_zipfile = "./raw_data/GlobalProduction_2019.1.0.zip", output_path="./data/FAO_Global_Production_2019.1.0.csv") #Only need to run the first time
}
fao_aqua <- readr::read_csv("./data/FAO_Global_Production_2019.1.0.csv") #Load FAO Aquaculture data 
fm_data <- readr::read_csv("../FishmealDemand/data/FishmealModeling_Oct062018.csv") #Load key diet Data
fm_distributions <- estimate_fmfo_distributions(data = fm_data) #Estimate distributions of key diet parameters 
readr::write_csv(fm_distributions, "./outputs/FishmealModeling_Output.csv")
