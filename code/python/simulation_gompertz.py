"""
Simulation du processus de diffusion de Gompertz stochastique (cas homogène)
par le schéma d'Euler-Maruyama.

Modèle (EDS d'Itô) :
    dX(t) = [ alpha * X(t) - beta * X(t) * ln X(t) ] dt + sigma * X(t) dW(t),
    X(0) = x0 > 0.

Le script :
  1. simule N trajectoires par Euler-Maruyama ;
  2. compare la moyenne empirique à la moyenne théorique E[X(t)] ;
  3. enregistre les données simulées (CSV) et une figure (PNG).

À lancer depuis la racine du dépôt :
    python code/python/simulation_gompertz.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # pas de fenêtre : la figure est enregistrée dans un fichier
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paramètres (tous enregistrés pour la reproductibilité)
# ---------------------------------------------------------------------------
SEED = 2026          # graine aléatoire
X0 = 1.0             # valeur initiale
ALPHA = 1.0          # taux de croissance intrinsèque
BETA = 0.2           # facteur de décélération
SIGMA = 0.1          # intensité du bruit
T = 30.0             # horizon de simulation
DT = 0.01            # pas de temps
N_PATHS = 200        # nombre de trajectoires simulées
N_SAVED = 20         # nombre de trajectoires enregistrées et tracées
SAVE_EVERY = 10      # on enregistre un point sur 10 (pas de 0.1)

ROOT = Path(__file__).resolve().parents[2]
OUT_DATA = ROOT / "data" / "simulated" / "gompertz_simulation_python.csv"
OUT_FIG = ROOT / "figures" / "gompertz_trajectoires_python.png"


def simulate_euler_maruyama(rng):
    """Retourne la grille de temps et une matrice (n_pas+1, N_PATHS) de trajectoires."""
    n_steps = int(round(T / DT))
    t = np.linspace(0.0, T, n_steps + 1)
    x = np.empty((n_steps + 1, N_PATHS))
    x[0] = X0
    sqrt_dt = np.sqrt(DT)
    for k in range(n_steps):
        xk = x[k]
        drift = ALPHA * xk - BETA * xk * np.log(xk)
        dw = sqrt_dt * rng.standard_normal(N_PATHS)
        # max(..., 1e-12) : garde-fou pour rester strictement positif
        x[k + 1] = np.maximum(xk + drift * DT + SIGMA * xk * dw, 1e-12)
    return t, x


def theoretical_mean(t):
    """E[X(t)] pour x0 déterministe.

    Y = ln X vérifie dY = (alpha - sigma^2/2 - beta Y) dt + sigma dW (Ornstein-Uhlenbeck),
    donc Y(t) est gaussien de moyenne m(t) et de variance v(t), et E[X(t)] = exp(m + v/2).
    """
    mu_inf = (ALPHA - SIGMA**2 / 2) / BETA
    m = mu_inf + (np.log(X0) - mu_inf) * np.exp(-BETA * t)
    v = SIGMA**2 * (1 - np.exp(-2 * BETA * t)) / (2 * BETA)
    return np.exp(m + v / 2)


def main():
    rng = np.random.default_rng(SEED)
    t, x = simulate_euler_maruyama(rng)
    mean_emp = x.mean(axis=1)
    mean_th = theoretical_mean(t)

    # Données simulées : temps, N_SAVED trajectoires, moyennes empirique et théorique
    idx = slice(None, None, SAVE_EVERY)
    df = pd.DataFrame({"t": t[idx]})
    for j in range(N_SAVED):
        df[f"x_{j + 1:02d}"] = x[idx, j]
    df["moyenne_empirique"] = mean_emp[idx]
    df["moyenne_theorique"] = mean_th[idx]
    OUT_DATA.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_DATA, index=False, float_format="%.6f")

    # Figure
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, x[:, :N_SAVED], color="steelblue", alpha=0.35, linewidth=0.8)
    ax.plot(t, mean_emp, color="black", linewidth=2,
            label=f"Moyenne empirique ({N_PATHS} trajectoires)")
    ax.plot(t, mean_th, color="crimson", linestyle="--", linewidth=2,
            label="Moyenne théorique E[X(t)]")
    ax.set_xlabel("Temps t")
    ax.set_ylabel("X(t)")
    ax.set_title("Processus de Gompertz stochastique (Euler-Maruyama)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=150)

    erreur = np.max(np.abs(mean_emp - mean_th) / mean_th)
    print(f"Données enregistrées : {OUT_DATA.relative_to(ROOT)}")
    print(f"Figure enregistrée   : {OUT_FIG.relative_to(ROOT)}")
    print(f"Écart relatif maximal moyenne empirique / théorique : {erreur:.2%}")


if __name__ == "__main__":
    main()
