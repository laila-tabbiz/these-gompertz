# -----------------------------------------------------------------------------
# Simulation du processus de diffusion de Gompertz stochastique (cas homogène)
# par le schéma d'Euler-Maruyama.
#
# Modèle (EDS d'Itô) :
#   dX(t) = [ alpha * X(t) - beta * X(t) * ln X(t) ] dt + sigma * X(t) dW(t),
#   X(0) = x0 > 0.
#
# Le script :
#   1. simule N trajectoires par Euler-Maruyama ;
#   2. compare la moyenne empirique à la moyenne théorique E[X(t)] ;
#   3. enregistre les données simulées (CSV) et une figure (PNG).
#
# À lancer depuis la racine du dépôt :
#   Rscript code/R/simulation_gompertz.R
# Uniquement R de base : aucun paquet à installer.
# -----------------------------------------------------------------------------

# Paramètres (tous enregistrés pour la reproductibilité)
seed       <- 2026   # graine aléatoire
x0         <- 1.0    # valeur initiale
alpha      <- 1.0    # taux de croissance intrinsèque
beta       <- 0.2    # facteur de décélération
sigma      <- 0.1    # intensité du bruit
T_final    <- 30     # horizon de simulation
dt         <- 0.01   # pas de temps
n_paths    <- 200    # nombre de trajectoires simulées
n_saved    <- 20     # nombre de trajectoires enregistrées et tracées
save_every <- 10     # on enregistre un point sur 10 (pas de 0.1)

out_data <- file.path("data", "simulated", "gompertz_simulation_R.csv")
out_fig  <- file.path("figures", "gompertz_trajectoires_R.png")

# Simulation d'Euler-Maruyama : matrice (n_steps + 1) x n_paths
simulate_euler_maruyama <- function() {
  n_steps <- round(T_final / dt)
  t <- seq(0, T_final, length.out = n_steps + 1)
  x <- matrix(NA_real_, nrow = n_steps + 1, ncol = n_paths)
  x[1, ] <- x0
  for (k in seq_len(n_steps)) {
    xk <- x[k, ]
    drift <- alpha * xk - beta * xk * log(xk)
    dw <- sqrt(dt) * rnorm(n_paths)
    # pmax(..., 1e-12) : garde-fou pour rester strictement positif
    x[k + 1, ] <- pmax(xk + drift * dt + sigma * xk * dw, 1e-12)
  }
  list(t = t, x = x)
}

# E[X(t)] pour x0 déterministe :
# Y = ln X suit un processus d'Ornstein-Uhlenbeck, donc Y(t) est gaussien
# de moyenne m(t) et de variance v(t), et E[X(t)] = exp(m + v/2).
theoretical_mean <- function(t) {
  mu_inf <- (alpha - sigma^2 / 2) / beta
  m <- mu_inf + (log(x0) - mu_inf) * exp(-beta * t)
  v <- sigma^2 * (1 - exp(-2 * beta * t)) / (2 * beta)
  exp(m + v / 2)
}

set.seed(seed)
sim <- simulate_euler_maruyama()
t <- sim$t
x <- sim$x
mean_emp <- rowMeans(x)
mean_th <- theoretical_mean(t)

# Données simulées
idx <- seq(1, length(t), by = save_every)
df <- data.frame(t = t[idx], x[idx, 1:n_saved])
names(df)[-1] <- sprintf("x_%02d", 1:n_saved)
df$moyenne_empirique <- mean_emp[idx]
df$moyenne_theorique <- mean_th[idx]
dir.create(dirname(out_data), recursive = TRUE, showWarnings = FALSE)
write.csv(round(df, 6), out_data, row.names = FALSE)

# Figure
dir.create(dirname(out_fig), recursive = TRUE, showWarnings = FALSE)
png(out_fig, width = 1200, height = 750, res = 150)
matplot(t, x[, 1:n_saved], type = "l", lty = 1,
        col = adjustcolor("steelblue", alpha.f = 0.35),
        xlab = "Temps t", ylab = "X(t)",
        main = "Processus de Gompertz stochastique (Euler-Maruyama)")
lines(t, mean_emp, col = "black", lwd = 2)
lines(t, mean_th, col = "red3", lwd = 2, lty = 2)
legend("bottomright",
       legend = c(sprintf("Moyenne empirique (%d trajectoires)", n_paths),
                  "Moyenne th\u00e9orique E[X(t)]"),
       col = c("black", "red3"), lwd = 2, lty = c(1, 2), bty = "n")
invisible(dev.off())

erreur <- max(abs(mean_emp - mean_th) / mean_th)
cat("Données enregistrées : ", out_data, "\n", sep = "")
cat("Figure enregistrée   : ", out_fig, "\n", sep = "")
cat(sprintf("Écart relatif maximal moyenne empirique / théorique : %.2f%%\n", 100 * erreur))
