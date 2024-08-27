# Leave-One-Population-Out Cross-Validation for Gradient Forest

library(gradientForest)
library(dplyr)
library(ggplot2)
library(reshape2)
library(viridis)


# import data
data <- readRDS("GFinput.rds")

# get env and snp names
envVars <- c("bio02", "bio03", "bio06", "bio09", "bio14", "bio15")
candSNPs <- setdiff(colnames(data), c("site", "X", "Y", envVars))

# function to run gradient forest model
run_gf_model <- function(train_data, test_data, envVars, candSNPs) {
  maxLevel <- log2(0.368 * nrow(train_data) / 2)

  gf_model <- gradientForest(
    cbind(train_data[, envVars], train_data[, candSNPs]),
    predictor.vars = envVars,
    response.vars = candSNPs,
    ntree = 500,
    maxLevel = maxLevel,
    trace = FALSE,
    corr.threshold = 0.50
  )
  

  return(gf_model)
}

# run full model
full_model <- run_gf_model(data, data, envVars, candSNPs)

# run Leave-One-Population-Out Cross-Validation
populations <- unique(data$site)
lopo_results <- list()

for (pop in populations) {
  train_data <- data[data$site != pop, ]
  test_data <- data[data$site == pop, ]
  
  # run gf model
  gf_results <- run_gf_model(train_data, test_data, envVars, candSNPs)
  

  lopo_results[[pop]] <- gf_results
}


# function to compare LOPO model with full model
compare_models <- function(lopo_result, full_model) {
  # compare R2
  lopo_r2 <- mean(lopo_result$result)
  full_r2 <- mean(full_model$result)
  
  # compare variable importance
  lopo_imp <- lopo_result$overall.imp
  full_imp <- full_model$overall.imp
  
  # predictions for LOPO model
  lopo_preds <- predict(lopo_result)
  
  # predictions for full model
  full_preds <- predict(full_model)

  # Calculate RMSE difference
  rmse_diff <- sapply(names(lopo_preds), function(var) {
    # corrected RMSE for LOPO model
    lopo_rmse <- sqrt(mean((lopo_preds[[var]] - mean(full_preds[[var]], na.rm = TRUE))^2, na.rm = TRUE))
    
    # RMSE for full model
    full_rmse <- sqrt(mean((full_preds[[var]] - mean(full_preds[[var]], na.rm = TRUE))^2, na.rm = TRUE))
    
    return(lopo_rmse - full_rmse)
  })

  
  return(list(
    r2_diff = lopo_r2 - full_r2,
    imp_diff = lopo_imp - full_imp,
    rmse_diff = rmse_diff
  ))
}


# run sensitivity analyses
sensitivity_results <- lapply(names(lopo_results), function(pop) {
  compare_models(lopo_results[[pop]], full_model)
})
names(sensitivity_results) <- names(lopo_results)


# summarize results
sensitivity_summary <- do.call(rbind, lapply(names(sensitivity_results), function(pop) {
  res <- sensitivity_results[[pop]]
  data.frame(
    Population = pop,
    R2_Difference = res$r2_diff,
    RMSE_Difference = mean(res$rmse_diff, na.rm = TRUE),
    Importance_Difference = mean(abs(res$imp_diff), na.rm = TRUE)
  )
}))


# Plot sensitivity results
lopo_sensitivity <- ggplot(sensitivity_summary, aes(x = Population)) +
  geom_bar(aes(y = R2_Difference, fill = "R2"), stat = "identity") +
  geom_point(aes(y = RMSE_Difference, color = "RMSE", shape = "RMSE")) +
  geom_line(aes(y = RMSE_Difference, color = "RMSE", group = 1)) +
  geom_point(aes(y = Importance_Difference, color = "Importance", shape = "Importance")) +
  geom_line(aes(y = Importance_Difference, color = "Importance", group = 1)) +
  scale_fill_manual(name = "", values = "grey", labels = c("R2")) +
  scale_color_manual(name = "", values = c("Importance" = "blue", "RMSE" = "red"), labels = c("Importance", "RMSE")) +
  scale_shape_manual(name = "", values = c("Importance" = 16, "RMSE" = 16), labels = c("Importance", "RMSE")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(x="Catchment", y = "Difference from Full Model") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


ggsave(plot = lopo_sensitivity, "lopo_sensitivity.pdf", width = 12, height = 8)


# heat map of RMSE per model, per environmental variable
rmse_data <- do.call(rbind, lapply(names(sensitivity_results), function(pop) {
  data.frame(
    Population = pop,
    t(sensitivity_results[[pop]]$rmse_diff)
  )
}))

rmse_melted <- melt(rmse_data, id.vars = "Population", variable.name = "EnvVar", value.name = "RMSE_diff")

# calculate mean absolute RMSE difference per population
population_sensitivity <- rmse_melted %>%
  group_by(Population) %>%
  summarize(mean_abs_rmse_diff = mean(abs(RMSE_diff), na.rm = TRUE)) %>%
  arrange(desc(mean_abs_rmse_diff))

top_3 <- head(population_sensitivity, 3)
bottom_3 <- tail(population_sensitivity, 3)



rmse_plot <- ggplot(rmse_melted, aes(x = EnvVar, y = Population, fill = RMSE_diff)) +
  geom_tile() +
  scale_fill_viridis(option = "plasma", direction = -1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Environmental Variable",
       y = "Catchment",
       fill = "RMSE Difference")

ggsave(plot = rmse_plot, "rmse_heatmap.pdf", width = 8, height = 8)


# heat map of variable importance across populations
imp_data <- do.call(rbind, lapply(names(sensitivity_results), function(pop) {
  data.frame(
    Population = pop,
    t(sensitivity_results[[pop]]$imp_diff)
  )
}))

imp_melted <- melt(imp_data, id.vars = "Population", variable.name = "EnvVar", value.name = "Importance_diff")

imp_plot <- ggplot(imp_melted, aes(x = EnvVar, y = Population, fill = Importance_diff)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis", direction = -1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Environmental Variable",
       y = "Catchment",
       fill = "Importance Difference")

ggsave(plot = imp_plot, "imp_heatmap.pdf", width = 8, height = 8)




# combined GF analyses
mod_list <- c(list(full_model = full_model), lopo_results)

# Use do.call to pass the arguments to combinedGradientForest
combined.mods <- do.call(combinedGradientForest, mod_list)

pdf(file = "lopo_predictor_density.pdf", height = 6, width = 8)
plot(combined.mods,plot.type="Predictor.Density")
dev.off()


