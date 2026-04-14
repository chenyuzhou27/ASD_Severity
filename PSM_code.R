library(openxlsx)
library(MatchIt)

# Read and prepare data
read_and_prepare_data <- function() {
  # Read datasets
  ct <- read.xlsx("./HATCH_metadata.xlsx", sheet = 'ct')
  asd <- read.xlsx("./HATCH_metadata.xlsx", sheet = 'asd')
  
  # Select columns and add ASD indicator
  columns <- c('Subject.ID', 'Cluster', 'Age', 'Gender', 'Height.(cm)', 'Weight.(kg)', 'BMI')
  ct <- ct[, columns]
  ct$ASD <- 0
  
  asd <- asd[, columns]
  asd$ASD <- 1
  
  # Merge and prepare data
  merged <- rbind(asd, ct)
  merged$Age <- as.numeric(merged$Age)
  
  # Scale numeric variables
  scaled_vars <- as.data.frame(lapply(merged[, c(3, 5:7)], scale))
  names(scaled_vars) <- c("Age", "Height", "Weight", "BMI")
  
  # Combine all data
  merged0 <- cbind(scaled_vars, 
                   Gender = factor(merged$Gender),
                   Subject.ID = merged$Subject.ID,
                   Cluster = merged$Cluster,
                   ASD = merged$ASD)
  
  return(merged0)
}

# Perform matching for a specific group
perform_matching <- function(data, cluster_type, gender, ratio) {
  # Filter by gender
  gender_data <- data[data$Gender == gender, ]
  
  # Filter by cluster type
  if(cluster_type == "Mild") {
    filtered_data <- gender_data[gender_data$Cluster %in% c("Mild", "control"), ]
  } else if(cluster_type == "Moderate") {
    filtered_data <- gender_data[gender_data$Cluster %in% c("Moderate", "control"), ]
  } else { # Severe
    filtered_data <- gender_data[gender_data$Cluster %in% c("Severe", "control"), ]
  }
  
  # Perform matching
  matching <- matchit(ASD ~ Age + Height + Weight + BMI, 
                      data = filtered_data, 
                      method = "nearest", 
                      ratio = ratio, 
                      replace = TRUE)
  
  return(match.data(matching))
}

# Extract control IDs that appear in multiple matches
get_control_ids <- function(matched_data_list, original_ct) {
  # Extract control IDs from each matching result
  id_lists <- lapply(matched_data_list, function(matched_data) {
    matched_data[matched_data$Cluster == "control", "Subject.ID"]
  })
  
  # Find IDs that appear at least twice
  all_ids <- unlist(id_lists)
  at_least_twice <- unique(all_ids[duplicated(all_ids) | duplicated(all_ids, fromLast = TRUE)])
  
  # Subset original control data
  control_subset <- subset(original_ct, Subject.ID %in% at_least_twice)
  
  return(control_subset)
}

# Main execution
merged_data <- read_and_prepare_data()

# Separate by gender
merged_ma <- merged_data[merged_data$Gender == '1', ]
merged_fe <- merged_data[merged_data$Gender == '2', ]

# Perform matching for males
matched_mild_ma <- perform_matching(merged_data, "Mild", '1', 1)
matched_moder_ma <- perform_matching(merged_data, "Moderate", '1', 1)
matched_seve_ma <- perform_matching(merged_data, "Severe", '1', 2)

# Perform matching for females
matched_mild_fe <- perform_matching(merged_data, "Mild", '2', 1)
matched_moder_fe <- perform_matching(merged_data, "Moderate", '2', 1)
matched_seve_fe <- perform_matching(merged_data, "Severe", '2', 5)

# Get control IDs and save results
original_ct <- read.xlsx("./HATCH_metadata.xlsx", sheet = 'ct')

# Male controls
male_matches <- list(matched_mild_ma, matched_moder_ma, matched_seve_ma)
male_controls <- get_control_ids(male_matches, original_ct)
write.xlsx(male_controls, './HATCH_control_ma.xlsx')

# Female controls
female_matches <- list(matched_mild_fe, matched_moder_fe, matched_seve_fe)
female_controls <- get_control_ids(female_matches, original_ct)
write.xlsx(female_controls, './HATCH_control_fe.xlsx')
