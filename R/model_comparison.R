# Compare two models
compare_models <- function(models_df, model1_name, model2_name) {
  
  # Filter models
  model1_data <- models_df %>% filter(model_name == model1_name)
  model2_data <- models_df %>% filter(model_name == model2_name)
  
  if (nrow(model1_data) == 0) {
    stop(paste("Model not found:", model1_name))
  }
  if (nrow(model2_data) == 0) {
    stop(paste("Model not found:", model2_name))
  }
  
  # Get protein lists
  proteins_m1 <- model1_data$protein_name
  proteins_m2 <- model2_data$protein_name
  
  # Find overlaps
  overlap_proteins <- intersect(proteins_m1, proteins_m2)
  unique_m1 <- setdiff(proteins_m1, proteins_m2)
  unique_m2 <- setdiff(proteins_m2, proteins_m1)
  
  # Create overlap data
  overlap_data <- data.frame(
    protein_name = overlap_proteins,
    coeff_model1 = model1_data$coefficient[model1_data$protein_name %in% overlap_proteins],
    coeff_model2 = model2_data$coefficient[model2_data$protein_name %in% overlap_proteins]
  )
  
  # Create detailed comparison
  all_proteins <- unique(c(proteins_m1, proteins_m2))
  detailed <- data.frame(
    protein_name = all_proteins,
    in_model1 = all_proteins %in% proteins_m1,
    in_model2 = all_proteins %in% proteins_m2
  )
  
  detailed <- detailed %>%
    left_join(
      model1_data %>% select(protein_name, coefficient) %>% 
        rename(coeff_model1 = coefficient),
      by = "protein_name"
    ) %>%
    left_join(
      model2_data %>% select(protein_name, coefficient) %>% 
        rename(coeff_model2 = coefficient),
      by = "protein_name"
    )
  
  return(list(
    overlap_proteins = overlap_proteins,
    unique_count_model1 = length(unique_m1),
    unique_count_model2 = length(unique_m2),
    overlap_data = overlap_data,
    detailed_comparison = detailed
  ))
}