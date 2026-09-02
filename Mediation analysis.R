############################ CMAverse Mediation Analysis ############################

# Clean environment and load libraries
rm(list = ls())
suppressMessages({
  library(data.table)
  library(dplyr)
  library(CMAverse)
  library(survival)
})

setwd("/path/to/your/project")

# Load data
covariate_df <- fread('./data/covariates.txt', header = TRUE)
metabolite_df <- fread('./data/metabolites_normalized.txt', header = TRUE) %>% 
  select(-FID)
exposure_df <- fread("./data/exposures.tsv", header = TRUE)
sample_ids <- fread('./data/sample_ids.txt', header = FALSE)$V1
triplets <- fread("./data/significant_triplets.txt", header = TRUE) %>% 
  as.data.frame()

# Preprocess covariates
covariate_df <- covariate_df %>%
  mutate(across(c("Sex", "Ethnicity", "Month", "Center"), as.factor))

# Load disease information
disease_info <- fread("./data/disease_endpoints.csv", 
                      select = c('NAME', 'LONGNAME'), header = TRUE)
disease_list <- disease_info$NAME

# Set seed for reproducibility
set.seed(123)

#------------------------------------------------------------------------------
# Helper function: Extract CMAverse results into named vector
#------------------------------------------------------------------------------
extract_cmest_results <- function(cmest_obj) {
  summarydf <- summary(cmest_obj)$summarydf
  
  summarydf <- summarydf %>%
    mutate(
      Estimate_fmt = sprintf("%.3f", round(Estimate, 3)),
      SE_fmt = sprintf("%.3f", round(`Std.error`, 3)),
      CIL_fmt = sprintf("%.3f", round(`95% CIL`, 3)),
      CIU_fmt = sprintf("%.3f", round(`95% CIU`, 3)),
      `95% CI` = paste0(Estimate_fmt, ' [', CIL_fmt, '-', CIU_fmt, ']')
    ) %>%
    select(Estimate = Estimate_fmt, SE = SE_fmt, `95% CI`, P = `P.val`)
  
  # Flatten to named vector
  result_vector <- c()
  for(row in rownames(summarydf)) {
    for(col in colnames(summarydf)) {
      name <- paste0(row, "_", col)
      result_vector <- c(result_vector, summarydf[row, col])
      names(result_vector)[length(result_vector)] <- name
    }
  }
  return(result_vector)
}

# Define result column names
info_cols <- c(
  "Rcde_Estimate", "Rcde_SE", "Rcde_95% CI", "Rcde_P",
  "Rpnde_Estimate", "Rpnde_SE", "Rpnde_95% CI", "Rpnde_P",
  "Rtnde_Estimate", "Rtnde_SE", "Rtnde_95% CI", "Rtnde_P",
  "Rpnie_Estimate", "Rpnie_SE", "Rpnie_95% CI", "Rpnie_P",
  "Rtnie_Estimate", "Rtnie_SE", "Rtnie_95% CI", "Rtnie_P",
  "Rte_Estimate", "Rte_SE", "Rte_95% CI", "Rte_P",
  "pm_Estimate", "pm_SE", "pm_95% CI", "pm_P"
)

#------------------------------------------------------------------------------
# Main mediation analysis loop
#------------------------------------------------------------------------------
for (i in 1:nrow(triplets)) {
  exposure <- as.character(triplets$Exposure[i])
  metabolite <- triplets$Metabolite[i]
  outcome <- triplets$Outcome[i]
  
  tryCatch({
   outcome_path <- paste0('./data/disease/', outcome, '.csv')
   outcome_df <- fread(outcome_path, 
                       select = c('eid', 'target_y', 'BL2Target_yrs'), 
                       header = TRUE) %>%
     mutate(eid = as.integer(eid), 
            target_y = as.numeric(target_y), 
            BL2Target_yrs = as.numeric(BL2Target_yrs)) %>%
     na.omit() %>%
     filter(BL2Target_yrs >= 0)
    
    # Merge all data
    data <- exposure_df %>%
      select(eid, Exposure = !!sym(exposure)) %>%
      filter(eid %in% sample_ids) %>%
      left_join(metabolite_df %>% 
                  select(IID, Metabolite = !!sym(metabolite)), 
                by = c('eid' = 'IID')) %>%
      left_join(outcome_df, by = 'eid') %>%
      left_join(covariate_df, by = c('eid' = 'IID')) %>%
      na.omit()

    # Run mediation analysis
    cmest_result <- cmest(
      data = data,
      exposure = 'Exposure',
      mediator = 'Metabolite',
      outcome = "BL2Target_yrs",
      event = "target_y",
      basec = valid_covars,
      mreg = list("linear"),
      yreg = 'coxph',
      model = "rb",
      EMint = FALSE,
      astar = 0, 
      a = 1, 
      mval = list(mediator = mean(data$Metabolite)),
      estimation = "imputation",
      inference = "bootstrap",
      nboot = 500
    )
    
    # Extract and store results
    triplets[i, info_cols] <- extract_cmest_results(cmest_result)
    
  }, error = function(e) {
    message(paste("Error in iteration", i, ":", e$message))
  })
}

# Save results
fwrite(triplets, "./results/mediation_results.txt", sep = '\t', row.names = FALSE, col.names = TRUE)

