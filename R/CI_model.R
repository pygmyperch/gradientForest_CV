# Leave-One-Population-Out Cross-Validation for Gradient Forest with Confidence Intervals

library(gradientForest)
library(dplyr)
library(ggplot2)
library(reshape2)
library(viridis)
library(parallel)



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

# function to compare LOPO model with full model
compare_models <- function(lopo_result, full_model) {
  lopo_r2 <- mean(lopo_result$result)
  full_r2 <- mean(full_model$result)
  
  lopo_imp <- lopo_result$overall.imp
  full_imp <- full_model$overall.imp
  
  lopo_preds <- predict(lopo_result)
  full_preds <- predict(full_model)
  
  rmse_diff <- sapply(names(lopo_preds), function(var) {
    lopo_rmse <- sqrt(mean((lopo_preds[[var]] - mean(full_preds[[var]], na.rm = TRUE))^2, na.rm = TRUE))
    full_rmse <- sqrt(mean((full_preds[[var]] - mean(full_preds[[var]], na.rm = TRUE))^2, na.rm = TRUE))
    return(lopo_rmse - full_rmse)
  })
  
  return(list(
    r2_diff = lopo_r2 - full_r2,
    imp_diff = lopo_imp - full_imp,
    rmse_diff = rmse_diff
  ))
}

# function to run the LOPO model
run_lopo_iteration <- function(data, envVars, candSNPs) {
  populations <- unique(data$site)
  full_model <- run_gf_model(data, data, envVars, candSNPs)
  
  lopo_results <- lapply(populations, function(pop) {
    train_data <- data[data$site != pop, ]
    test_data <- data[data$site == pop, ]
    gf_results <- run_gf_model(train_data, test_data, envVars, candSNPs)
    compare_models(gf_results, full_model)
  })
  
  names(lopo_results) <- populations
  return(lopo_results)
}

# run n iterations
n_iterations <- 1000
#n_cores <- detectCores() - 1
n_cores <- 2

cat("Running", n_iterations, "iterations on", n_cores, "cores...\n")

results <- mclapply(1:n_iterations, function(i) {
  if (i %% 100 == 0) cat("Iteration", i, "\n")
  run_lopo_iteration(data, envVars, candSNPs)
}, mc.cores = n_cores)

# estimate mean and 95% CIs for R2, RMSE, and variable importance across all iterations
process_results <- function(results) {
  populations <- names(results[[1]])
  
  processed <- lapply(populations, function(pop) {
    pop_results <- lapply(results, `[[`, pop)
    
    r2_diffs <- sapply(pop_results, function(x) x$r2_diff)
    rmse_diffs <- do.call(rbind, lapply(pop_results, function(x) x$rmse_diff))
    imp_diffs <- do.call(rbind, lapply(pop_results, function(x) x$imp_diff))
    
    list(
      r2 = list(
        mean = mean(r2_diffs),
        ci = quantile(r2_diffs, c(0.025, 0.975))
      ),
      rmse = apply(rmse_diffs, 2, function(x) {
        c(mean = mean(x), ci_lower = quantile(x, 0.025), ci_upper = quantile(x, 0.975))
      }),
      imp = apply(imp_diffs, 2, function(x) {
        c(mean = mean(x), ci_lower = quantile(x, 0.025), ci_upper = quantile(x, 0.975))
      })
    )
  })
  
  names(processed) <- populations
  return(processed)
}

processed_results <- process_results(results)

# plotting functions
plot_rmse_differences <- function(processed_results) {
  rmse_data <- do.call(rbind, lapply(names(processed_results), function(pop) {
    pop_data <- processed_results[[pop]]$rmse
    data.frame(
      Population = pop,
      Variable = colnames(pop_data),
      RMSE_diff = pop_data["mean", ]
    )
  }))
  
  ggplot(rmse_data, aes(x = Variable, y = Population, fill = RMSE_diff)) +
    geom_tile() +
    scale_fill_viridis(option = "plasma", direction = -1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Environmental Variable", y = "Population", fill = "Mean RMSE Difference")
}

plot_importance_differences <- function(processed_results) {
  imp_data <- do.call(rbind, lapply(names(processed_results), function(pop) {
    pop_data <- processed_results[[pop]]$imp
    data.frame(
      Population = pop,
      Variable = colnames(pop_data),
      Importance_diff = pop_data["mean", ]
    )
  }))
  
  ggplot(imp_data, aes(x = Variable, y = Population, fill = Importance_diff)) +
    geom_tile() +
    scale_fill_viridis(option = "viridis", direction = -1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Environmental Variable", y = "Population", fill = "Mean Importance Difference")
}

# plot sensitivity stats
plot_sensitivity_with_ci <- function(processed_results) {
  sensitivity_data <- do.call(rbind, lapply(names(processed_results), function(pop) {
    res <- processed_results[[pop]]
    data.frame(
      Population = pop,
      R2_Difference = res$r2$mean,
      R2_Lower = res$r2$ci["2.5%"],
      R2_Upper = res$r2$ci["97.5%"],
      RMSE_Difference = mean(res$rmse["mean", ]),
      RMSE_Lower = mean(res$rmse["ci_lower.2.5%", ]),
      RMSE_Upper = mean(res$rmse["ci_upper.97.5%", ]),
      Importance_Difference = mean(res$imp["mean", ]),
      Importance_Lower = mean(res$imp["ci_lower.2.5%", ]),
      Importance_Upper = mean(res$imp["ci_upper.97.5%", ])
    )
  }))
  
  # set pop factor to original order
  sensitivity_data$Population <- factor(sensitivity_data$Population, levels = unique(sensitivity_data$Population))
  
  sensitivity_data_long <- reshape2::melt(sensitivity_data, id.vars = "Population")
  sensitivity_data_long$Metric <- gsub("_.*", "", sensitivity_data_long$variable)
  sensitivity_data_long$Type <- gsub(".*_", "", sensitivity_data_long$variable)
  
  sensitivity_data_wide <- reshape2::dcast(sensitivity_data_long, Population + Metric ~ Type, value.var = "value")
  
  # set colors
  color_values <- c("R2" = "grey", "RMSE" = "red", "Importance" = "blue")
  
  ggplot(sensitivity_data_wide, aes(x = Population, y = Difference, group = Metric)) +
    geom_col(data = subset(sensitivity_data_wide, Metric == "R2"), 
             aes(fill = Metric), alpha = 0.5, position = position_dodge(width = 0.5)) +
    geom_point(data = subset(sensitivity_data_wide, Metric != "R2"),
               aes(color = Metric),
               position = position_dodge(width = 0.5)) +
    geom_errorbar(data = subset(sensitivity_data_wide, Metric != "R2"),
                  aes(ymin = Lower, ymax = Upper, color = Metric), 
                  position = position_dodge(width = 0.5), width = 0.2) +
    geom_errorbar(data = subset(sensitivity_data_wide, Metric == "R2"),
                  aes(ymin = Lower, ymax = Upper), 
                  position = position_dodge(width = 0.5), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    scale_color_manual(values = color_values, 
                       breaks = c("R2", "RMSE", "Importance"),
                       labels = c(expression(R^2), "RMSE", "Importance")) +
    scale_fill_manual(values = color_values,
                      breaks = c("R2", "RMSE", "Importance"),
                      labels = c(expression(R^2), "RMSE", "Importance")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.title = element_blank()) +
    labs(x = "Catchment", y = "Difference from Full Model") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(fill = guide_legend(order = 1),
           color = guide_legend(order = 2))
}




# save plots
rmse_plot <- plot_rmse_differences(processed_results)
ggsave("rmse_differences_heatmap.pdf", rmse_plot, width = 12, height = 8)

imp_plot <- plot_importance_differences(processed_results)
ggsave("importance_differences_heatmap.pdf", imp_plot, width = 12, height = 8)

sensitivity_plot <- plot_sensitivity_with_ci(processed_results)
ggsave("sensitivity_with_ci.pdf", sensitivity_plot, width = 12, height = 8)

# save result object
saveRDS(processed_results, "processed_lopo_results.rds")
