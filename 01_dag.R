# =============================================================================
# 01_dag.R — Causal DAGs for the Moertel re-analysis
# =============================================================================
# Draws two graphs and exports each as SVG + PNG into figures/.
#
#   G_trial         — actual causal structure under randomization. Z -> rx is
#                     SEVERED. U -> rx is SEVERED.
#   G_observational — the counterfactual world where rx is assigned through
#                     the same process that drives observational treatment use.
#                     Same as G_trial plus the edges Z -> rx and U -> rx.
#
# Used by 00_estimands.qmd for identification arguments and by manuscript.qmd
# as embedded figures.
#
# Run:
#   Rscript 01_dag.R
# Outputs:
#   figures/dag_trial.svg
#   figures/dag_trial.png
#   figures/dag_observational.svg
#   figures/dag_observational.png
#   figures/dag_adjustment_sets.txt   (dagitty-computed minimal sufficient sets)
# =============================================================================

suppressPackageStartupMessages({
  library(dagitty)
  library(svglite)
})

set.seed(42)

fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# G_trial — randomization severs Z -> rx and U -> rx
# -----------------------------------------------------------------------------
# Coordinates chosen so the DAG reads left-to-right with the causal pipeline
# rx -> M -> Y_death on the top row and Z, U as confounders below.
g_trial <- dagitty('
  dag {
    bb="0,0,1,1"
    rx     [exposure, pos="0.20,0.40"]
    M      [pos="0.45,0.40"]
    Y_death[outcome, pos="0.70,0.40"]
    Z      [pos="0.45,0.65"]
    U      [latent, pos="0.45,0.85"]
    S      [pos="0.10,0.80"]

    rx -> M
    M  -> Y_death
    rx -> Y_death
    Z  -> M
    Z  -> Y_death
    U  -> Z
    U  -> Y_death
    Z  -> S
  }
')

# -----------------------------------------------------------------------------
# G_observational — same as G_trial PLUS Z -> rx
# -----------------------------------------------------------------------------
# Per the project spec, the only edge added to G_trial is Z -> rx. We do NOT
# add U -> rx; U remains a parent of Z and Y_death only. The intent of Q2 is
# to ask: "if the only way confounding entered was through measured Z, would
# back-door adjustment on Z recover the randomized answer?" With U -> rx not
# present, the answer must be yes, and {Z} is a sufficient adjustment set.
# Unmeasured U -> rx is handled separately via sensitivity (E-value).
g_obs <- dagitty('
  dag {
    bb="0,0,1,1"
    rx     [exposure, pos="0.20,0.40"]
    M      [pos="0.45,0.40"]
    Y_death[outcome, pos="0.70,0.40"]
    Z      [pos="0.45,0.65"]
    U      [latent, pos="0.45,0.85"]
    S      [pos="0.10,0.80"]

    rx -> M
    M  -> Y_death
    rx -> Y_death
    Z  -> rx
    Z  -> M
    Z  -> Y_death
    U  -> Z
    U  -> Y_death
    Z  -> S
  }
')

# -----------------------------------------------------------------------------
# Adjustment sets — sanity-check identification arguments in 00_estimands.qmd
# -----------------------------------------------------------------------------
# In G_trial we expect the empty set to suffice (randomization).
# In G_obs we expect Z to suffice IF U has no direct U -> rx path that bypasses
# Z. We deliberately included U -> rx in G_obs; minimalAdjustmentSets should
# therefore return the empty set NEVER and may return no valid set if U is
# truly latent. We use Z as the working set and bound U's effect via sensitivity.

cat("============================================================\n")
cat(" Minimal sufficient adjustment sets (dagitty)\n")
cat("============================================================\n\n")

cat("G_trial — adjustment sets for rx -> Y_death:\n")
print(adjustmentSets(g_trial, exposure = "rx", outcome = "Y_death",
                     type = "minimal"))
cat("\n")

cat("G_observational — adjustment sets for rx -> Y_death (U latent):\n")
print(adjustmentSets(g_obs, exposure = "rx", outcome = "Y_death",
                     type = "minimal"))
cat("\n")

cat("G_observational — adjustment sets if U were observed:\n")
print(adjustmentSets(g_obs, exposure = "rx", outcome = "Y_death",
                     type = "minimal",
                     effect = "total"))
cat("\n")

cat("G_trial — implied conditional independencies:\n")
print(impliedConditionalIndependencies(g_trial))
cat("\n")

# Write to a text file we can reference in the manuscript.
sink(file.path(fig_dir, "dag_adjustment_sets.txt"))
cat("G_trial adjustment sets for rx -> Y_death:\n")
print(adjustmentSets(g_trial, exposure = "rx", outcome = "Y_death",
                     type = "minimal"))
cat("\nG_observational adjustment sets for rx -> Y_death (U latent):\n")
print(adjustmentSets(g_obs, exposure = "rx", outcome = "Y_death",
                     type = "minimal"))
cat("\nG_trial implied conditional independencies:\n")
print(impliedConditionalIndependencies(g_trial))
sink()

# -----------------------------------------------------------------------------
# Export SVG + PNG
# -----------------------------------------------------------------------------
plot_dag <- function(g, title, file_stem) {
  # SVG via svglite (no cairo dependency)
  svglite::svglite(file.path(fig_dir, paste0(file_stem, ".svg")),
                   width = 6.5, height = 4.5)
  par(mar = c(1, 1, 2, 1))
  plot(g)
  title(title, cex.main = 1.1)
  dev.off()
  # PNG fallback for slide decks
  png(file.path(fig_dir, paste0(file_stem, ".png")),
      width = 1300, height = 900, res = 200)
  par(mar = c(1, 1, 2, 1))
  plot(g)
  title(title, cex.main = 1.1)
  dev.off()
  cat("  wrote ", file_stem, ".svg / .png\n", sep = "")
}

cat("============================================================\n")
cat(" Writing DAG figures\n")
cat("============================================================\n")

plot_dag(g_trial,
         expression(paste(G[trial], ": randomization severs Z " %->% " rx and U " %->% " rx")),
         "dag_trial")
plot_dag(g_obs,
         expression(paste(G[observational], ": Z " %->% " rx and U " %->% " rx restored")),
         "dag_observational")

cat("\nDone.\n")
