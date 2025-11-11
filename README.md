# GEE Plotting Toolkit

**Author:** Dr. Fariborz Aref  
**Scope:** geepack models with binary or Gaussian outcomes  
**License:** MIT  

### What it gives you
- Forest plot of robust coefficients with 95 percent CI  
- Residuals vs fitted diagnostic  
- Marginal effect curve with confidence band and optional facet  
- ROC AUC for binomial models  
- Cluster influence plot using leave one cluster out

### Structure

### Quick use
```r
source("GEE_Toolkit/gee_plots.R")

library(geepack)
m <- geeglm(y ~ x + z + age, id = id, data = df,
            family = binomial("logit"), corstr = "exchangeable")

gee_forest_plot(m, "figs/gee_forest.png")
gee_resid_plot(m,  "figs/gee_residuals.png")
gee_marginal_curve(m, df, var = "x", by = "z",
                   save_path = "figs/gee_marginal.png")
gee_roc_plot(m, df, yvar = "y", save_path = "figs/gee_roc.png")
gee_cluster_influence(m, df, cluster_id = "id",
                      target_term = "x", save_path = "figs/gee_influence.png")
