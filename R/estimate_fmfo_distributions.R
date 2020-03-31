#' Estimates distributions of species and functional group key feed parameters
#'
#' \code{estimate_fmfo_distribution} #Takes in a dataframe of feed information of aquaculture species (see example table)
#' Outputs a tidy format of main feed parameters (FM, FO, FCR, etc.) with
#' Mean, Standard Deviation, Lower bound, and Upper Bound for each country for each taxon and functional group

#' @param data dataframe
#' @param output_format string of either \code{"tidy"} or "NA"
#' @return tibble with bootstrapped data
#' @export
#' @examples
#' estimate_fmfo_distribution(fmfo_data)
#' estimate_fmfo_distribution(fmfo_data, output_format=NA)
#' @importFrom dplyr filter select bind_rows pull
#' @importFrom tidyr pivot_longer
#' @importFrom boot boot
#' @importFrom boot boot.ci
#' @importFrom tibble tibble
#' @importFrom magrittr `%>%`

####Description ####

#Original by Tim Cashion (t.cashion@oceans.ubc.ca; trcashion@gmail.com)
#March 30, 2020

estimate_fmfo_distributions <- function(data=NA, output_format="tidy"){
  set.seed(0)
  data$Species[which(is.na(data$Species)==TRUE)] <- data$SpeciesGroup[which(is.na(data$Species))]
  data <- data %>% filter(SpeciesGroup != "Unfed species") #Drop these as they are not useful for this stage of analysis and will create erroneous values

  #Alter species groups to have fewer groups. Original groups from reported sources were maintained in raw data
  data$SpeciesGroup[grepl(data$SpeciesGroup, pattern='shrimp', ignore.case = T)] <- "Shrimps"
  data$SpeciesGroup[grepl(data$Species, pattern='shrimp', ignore.case = T)] <- "Shrimps"
  data$SpeciesGroup[grepl(data$Species, pattern='prawn', ignore.case = T)] <- "Shrimps"
  data$SpeciesGroup[grepl(data$SpeciesGroup, pattern='Crustaceans', ignore.case = T)] <- "Crustaceans"

  #Current method only uses different weighted reported values and cannot intake lower bounds, upper bounds, or standard deviations
  #To include the data from these, these are treated as 'reported' values and their rows are added
  # variables <- c("On_Feed", "FCR", "FM", "FO")
  # distributions <- c("_SD", "_LL", "_UL")
  # new_vars <- character()
  # for (var in variables) {
  #   x <- paste(var, distributions, sep="")
  #   new_vars <- c(new_vars, x)
  # }
  # extra_data <- data[0,]
  # var <- new_vars[11]
  # for (var in new_vars){
  #   df <- data %>% filter(is.na(var)==F)
  #   df <- data %>% filter(var>=0)
  #   df$DataType <- "Average"
  #
  #   extra_data <- rbind(extra_data, df)
  # }
  #
  # rbind(data, extra_data)

  #Expected output calculations:
  # test <- data %>% select(Species, SpeciesGroup, Area) %>% unique()
  # nrow(test) * 4 #4 variables of interest for each unqiue combo

  #Weights should be re-defined based on accuracy of weight, time period, and predicted values
  #Score out of 5
  #Base score is 1
  #Is it country specifc? #This could be modified to region as if it is no, a near neighbour is given the same weight as a different continent
  #Is it species specific?
  #Is it average for the sector, or one study?
  #Is it recent +/- 5 years?
  #Assumption: Predicted values are treated equally to observed values
  #Assumption: Our aim is to weight value for each species-country more specifically than global averages
  #Assumption: Single studies are less descriptive of the entire aquaculture sector for a country than an average.
  #Assumption: Technological changes outside the species group have no effect or bearing on values within a different species group.
  #Weight is then determined by converting these to decimal fractions out of 1.00

  #In this way, each species-country would have its own weighting applied to it.

  #Bootstrap looped groups
  groups <- unique(data$SpeciesGroup) #Define unique species groups to search for
  variables <- c("On_Feed", "FCR", "FM", "FO") #Variables to estimate at present

  new_data <- tibble(group = NA,
                     country = NA,
                     taxon = NA,
                     variable = NA,
                     mean = NA,
                     sd = NA,
                     ll = NA,
                     ul = NA) %>%
    na.omit()

  #Loop to create new data for all groups and their species and countries ####

  for (group in groups){
    df <- data %>% filter(SpeciesGroup==group)

    #Within each species group, find where we need to estimate for different countries and taxa
    taxa <- df %>% pull(Species) %>% na.omit() %>% unique()
    countries <- df %>% pull(Area) %>% na.omit() %>% unique()
    for (taxon in taxa){
      for (country in countries){

        #Weight relevant data based on specificity to taxon and country
        df$Weight <- 1
        df$Weight[which(df$Species==as.character(taxon) & df$Species!=df$SpeciesGroup)] <- df$Weight[which(df$Species==taxon & df$Species!=df$SpeciesGroup)] + 1
        df$Weight[which(df$Area==country & df$Area != "World")] <- df$Weight[which(df$Area==country & df$Area != "World")] + 1
        df$Weight[which(df$DataType=="Average")] <- df$Weight[which(df$DataType=="Average")] + 1
        df$Weight[which(df$Year >=2010 & df$Year <=2020)] <- df$Weight[which(df$Year >=2010 & df$Year <=2020)] + 1
        df$NormWeight <- df$Weight/ sum(df$Weight)

        #Loop over variables:
        for (variable in variables){
          sample <- df %>% select(variable, NormWeight) %>% na.omit()
          if (nrow(sample)==0){ #If no data for that variable within the sample, give an estimate of 0
            sample <- df %>% select(variable, NormWeight)
            sample[,1] <- 0
          }

          #If data are present...
          if (nrow(sample) > 1){
            #And if observations have variance (n> 2 where they are not equal to each other)
            if (var(sample %>% pull(variable), na.rm=TRUE)!=0){
              boot = boot(sample %>% pull(variable),
                          statistic=boot_mean_var,
                          R=1000,
                          weights=sample %>% pull(NormWeight))
              ci <- boot.ci(boot, conf=0.95, type="stud")
              mean <- mean(boot$t[,1])
              sd <- sqrt(var(boot$t[,1]))
              ll <- ci$student[4]
              ul <- ci$student[5]

            } else { #If there is no variance, cannot bootstrap
              #Mean = existing value with sd of 0.

              mean <- mean(sample %>% pull(variable), na.rm=TRUE)
              sd <- 0
              ll <- mean(sample %>% pull(variable), na.rm=TRUE)
              ul <- mean(sample %>% pull(variable), na.rm=TRUE)
            }
            #Fill in country and taxon information if not given previously
            if(length(country)==0){
              country <- "World"
            }
            if(length(taxon)==0){
              taxon <- ""
            }
            #And if there is only one sample, the same rules apply if there is no variance
          } else {
              mean <- mean(sample %>% pull(variable), na.rm=TRUE)
              sd <- 0
              ll <- mean(sample %>% pull(variable), na.rm=TRUE)
              ul <- mean(sample %>% pull(variable), na.rm=TRUE)
          }

          #Create new entry with boostrapped estimates
          x <- tibble("group" = group,
                      "country" = country,
                      "taxon" = taxon,
                      "variable" = variable,
                      "mean" = mean,
                      "sd" = sd,
                      "ll" = ll,
                      "ul" = ul)

          #Bind bootstrap estimate to output dataframe
          new_data <- bind_rows(new_data, x)
        }
      }
    }
  }
  if (output_format=="tidy"){
    new_data_tidy <- new_data %>%
      pivot_longer(cols=c(mean, sd, ll, ul), names_to = "measure", "value")
    return(new_data_tidy)
  } else{
    return(new_data)
  }
}
