# gradientForest_CV

Code to run the Leave-One-Population-Out (LOPO) cross-validation to assess sensitivity of gradient forest model to the exclusion of individual catchments for:

> Booth, E, Brauer, CJ, Sandoval-Castillo, J, Harrisson, K, Rourke, ML, Attard, CM, Gilligan, DM, Tonkin, Z, Thiem, J, Zampatti, B, Unmack, PJ, & Beheregaray, LB. (2024). Genomic vulnerability to climate change of an Australian migratory freshwater fish, the golden perch (Macquaria ambigua). *Molecular Ecology*. Under review.

* `GF_LOPOCV.R`: quick check for outlying populations/catchments/regions
* `CI_model.R`: run n iterations of the LOPO models to estimate mean and 95% confidence intervals for R^2, RMSE, and variable importance

Plot results as below:

<table>
  <tr>
    <td>
      <p align="center">
        <img src="../main/images/sensitivity.jpg" alt="Figure 1">
      </p>
      <p align="left">
        <em><strong>Figure 1.</strong> Sensitivity of the full gradient forest model to the exclusion of individual catchments based on 1000 Leave-One-Population-Out (LOPO) cross-validation models, assessed by difference in R^2 values (and 95% CIs) between the LOPO model and the full model. Positive R^2 suggests improved model fit while negative R^2 suggests decreased model fit when the population is excluded. RMSE (Root Mean Square Error; red points) shows the mean difference in RMSE between LOPO model and full model predictions (the environmental variables). Importance (blue points) represents the mean (across all environmental variables) absolute difference in variable importance between each LOPO model and the full model.</em>
      </p>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td>
      <p align="center">
        <img src="../main/images/rmse_heatmap.jpg" alt="Figure 2">
      </p>
      <p align="left">
        <em><strong>Figure 2.</strong> Root Mean Square Error (RMSE) differences by catchment and environmental variable. Differences in mean RMSE (and 95% CIs) based on 1000 Leave-One-Population-Out (LOPO) models and the full model for each combination of catchment and environmental variable. Each row is the catchment that was held out in the LOPO model training. Darker cells suggest predictions of the environmental variable are more sensitive to the exclusion of that population. Patterns of dark cells across rows indicate catchments that strongly influence the full model's predictive accuracy when excluded. Patterns down columns suggest environmental variables that are consistently sensitive across different population exclusions.</em>
      </p>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td>
      <p align="center">
        <img src="../main/images/importance_heatmap.jpg" alt="Figure 3">
      </p>
      <p align="left">
        <em><strong>Figure 3.</strong> Mean (and 95% CIs) variable importance differences based on 1000 Leave-One-Population-Out (LOPO) models and the full model for each combination of catchment and environmental variable. Each row is the catchment that was held out in the LOPO model training. Darker cells highlight combinations where the importance of an environmental variable changes substantially when a specific catchment is excluded from the model. Lighter cells indicate that the variable's importance is relatively insensitive to the exclusion of that catchment. Horizontal patterns indicate catchments that influence the importance of multiple environmental variables. Vertical patterns identify environmental variables that are sensitive to which catchments are excluded.</em>
      </p>
    </td>
  </tr>
</table>
