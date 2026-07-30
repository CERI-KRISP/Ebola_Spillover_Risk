library(biomod2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(fastshap)
library(shapviz)
library(ggrepel)

dir.create("shap", recursive = TRUE, showWarnings = FALSE)


# BACKGROUND / EXPLANATORY DATA

X <- get_formal_data(myBiomodModels, "expl.var")
X <- as.data.frame(X)

cat("Background data for SHAP:", nrow(X), "rows x", ncol(X), "predictors\n")
print(names(X))


# UNIFIED, SCALE-ROBUST PREDICTION WRAPPER

predict_biomod <- function(object, newdata) {
  p <- as.numeric(predict(object, newdata))
  if (max(p, na.rm = TRUE) > 1) p <- p / 1000
  p
}

# RETRIEVE MODELS 
ens_all_names   <- get_built_models(myBiomodEM)
ens_wmean_name  <- grep("EMwmean", ens_all_names, value = TRUE)[1]
stopifnot(!is.na(ens_wmean_name))
cat("Ensemble model selected:", ens_wmean_name, "\n")

base_all_names  <- get_built_models(myBiomodModels)
algo_names      <- c("GLM", "GBM", "RF", "MAXNET")

rep_names <- sapply(algo_names, function(a) {
  hit <- grep(paste0("allRun_", a, "$"), base_all_names, value = TRUE)
  if (length(hit) == 0) hit <- grep(paste0("_", a, "$"), base_all_names, value = TRUE)[1]
  hit[1]
})
names(rep_names) <- algo_names
cat("Base models selected for robustness check:\n")
print(rep_names)


# SHAP 

BIOMOD_LoadModels(myBiomodEM, full.name = ens_wmean_name)
ens_model <- get(ens_wmean_name)

set.seed(42)
message("Computing ensemble SHAP values (this may take a few minutes)...")
shap_ens <- fastshap::explain(
  object       = ens_model,
  X            = X,
  pred_wrapper = predict_biomod,
  nsim         = 200     # increase to 150-200 for the final manuscript run
)
shap_ens_mat <- as.matrix(shap_ens)
write.csv(shap_ens_mat, "shap/shap_values_ensemble.csv", row.names = FALSE)

# Global importance table (mean |SHAP|), ranked
imp_ens <- data.frame(
  variable = colnames(shap_ens_mat),
  mean_abs_shap = colMeans(abs(shap_ens_mat))
) %>% arrange(desc(mean_abs_shap))
write.csv(imp_ens, "shap/importance_ensemble.csv", row.names = FALSE)
print(imp_ens)

baseline <- mean(predict_biomod(ens_model, X))
reconstructed <- baseline + rowSums(shap_ens_mat)
actual <- predict_biomod(ens_model, X)

cat("Max reconstruction error:", max(abs(reconstructed - actual)), "\n")
cat("Correlation (should be ~1):", cor(reconstructed, actual), "\n")


# SHAP FOR EACH BASE ALGORITHM 
shap_by_algo <- list()

for (a in algo_names) {
  message("Computing SHAP for: ", a)
  BIOMOD_LoadModels(myBiomodModels, full.name = rep_names[[a]])
  mod <- get(rep_names[[a]])
  
  shap_a <- tryCatch(
    fastshap::explain(object = mod, X = X, pred_wrapper = predict_biomod, nsim = 100),
    error = function(e) { message("  skipped ", a, ": ", e$message); NULL }
  )
  if (!is.null(shap_a)) {
    shap_by_algo[[a]] <- as.matrix(shap_a)
  }
}

imp_by_algo <- bind_rows(lapply(names(shap_by_algo), function(a) {
  m <- shap_by_algo[[a]]
  data.frame(algo = a, variable = colnames(m), mean_abs_shap = colMeans(abs(m)))
}))

imp_by_algo_ens <- bind_rows(
  imp_by_algo,
  imp_ens %>% mutate(algo = "Ensemble")
)
write.csv(imp_by_algo_ens, "shap/importance_by_algorithm.csv", row.names = FALSE)

#PLOTS

p_imp <- ggplot(imp_ens, aes(x = reorder(variable, mean_abs_shap), y = mean_abs_shap)) +
  geom_col(fill = "#932667") +
  coord_flip() +
  labs(x = NULL, y = "Mean |SHAP value|")+
  theme_minimal(base_size = 12)
p_imp
ggsave("shap/fig_ensemble_importance.png", p_imp, width = 7, height = 5, dpi = 300)


baseline_ens <- mean(predict_biomod(ens_model, X))
sv_ens <- shapviz(shap_ens_mat, X = X, baseline = baseline_ens)
p_summary <- sv_importance(sv_ens, kind = "beeswarm") +
  theme_minimal(base_size = 12)
p_summary
ggsave("shap/fig_ensemble_shap_summary.png", p_summary, width = 7, height = 5, dpi = 300)


top_vars <- head(imp_ens$variable, 4)
for (v in top_vars) {
  p_dep <- sv_dependence(sv_ens, v = v) +
    labs(title = paste("SHAP dependence:", v))
  ggsave(paste0("shap/fig_dependence_", v, ".png"), p_dep,
         width = 6, height = 4.5, dpi = 300)
}


p_cross <- ggplot(imp_by_algo_ens,
                  aes(x = reorder(variable, mean_abs_shap, FUN = median),
                      y = mean_abs_shap, fill = algo)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = NULL, y = "Mean |SHAP value|", fill = "Model",
       title = "SHAP importance: ensemble vs. base algorithms") +
  theme_minimal(base_size = 12)
ggsave("shap/fig_cross_algorithm_consistency.png", p_cross,
       width = 8, height = 6, dpi = 300)

#CROSS-CHECK AGAINST BIOMOD2's NATIVE PERMUTATION IMPORTANCE

raw_perm <- get_variables_importance(myBiomodEM)
cat("\nget_variables_importance() output structure (for column-name check):\n")
print(head(raw_perm))

name_col  <- intersect(c("full.name", "merged.by.algo", "algo"), names(raw_perm))[1]
var_col   <- intersect(c("expl.var", "var", "variable"), names(raw_perm))[1]
score_col <- intersect(c("var.imp", "value", "imp", "importance"), names(raw_perm))[1]
stopifnot(!is.na(name_col), !is.na(var_col), !is.na(score_col))

raw_perm_wmean <- raw_perm[grepl("EMwmean", raw_perm[[name_col]]), ]

perm_imp <- raw_perm_wmean %>%
  group_by(.data[[var_col]]) %>%
  summarise(mean_perm_imp = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
  rename(variable = all_of(var_col))

comparison <- imp_ens %>%
  left_join(perm_imp, by = "variable") %>%
  arrange(desc(mean_abs_shap))
write.csv(comparison, "shap/shap_vs_permutation_comparison.csv", row.names = FALSE)

spearman_rho <- cor(rank(comparison$mean_abs_shap), rank(comparison$mean_perm_imp),
                    method = "spearman", use = "complete.obs")
cat("\nSpearman rank correlation between SHAP and permutation importance:",
    round(spearman_rho, 3), "\n")

p_compare <- ggplot(comparison, aes(x = mean_perm_imp, y = mean_abs_shap, label = variable)) +
  geom_point(size = 2.5, colour = "#9e2a2b") +
  geom_text_repel(size = 3) +
  labs(x = "biomod2 permutation importance", y = "Mean |SHAP value|",
       title = paste0("SHAP vs. permutation importance (Spearman rho = ",
                      round(spearman_rho, 2), ")")) +
  theme_minimal(base_size = 12)
p_compare
ggsave("outputs/shap/fig_shap_vs_permutation.png", p_compare, width = 6.5, height = 5.5, dpi = 300)

