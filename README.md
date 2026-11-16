# Modélisation stochastique de systèmes gompertziens corrélés

Dépôt de travail de la thèse de doctorat de **Laila Tabbiz** : code de simulation et d'estimation, données simulées, bibliographie et plan de gestion des données.

- **Doctorante** : Laila Tabbiz — ORCID [0009-0006-4969-4455](https://orcid.org/0009-0006-4969-4455)
- **Direction de thèse** : Pr El Kettani Moummou
- **Laboratoire** : LaR2A — Analyse Numérique, Analyse Non Linéaire et Applications (UAE/U04FS)
- **Établissement** : Université Abdelmalek Essaâdi, Faculté des Sciences de Tétouan, Département de Mathématiques
- **Financement** : bourse doctorale PASS, CNRST

## Objet

La thèse étudie des extensions multivariées et corrélées du processus de diffusion de Gompertz stochastique, ainsi que l'inférence statistique associée (estimation par maximum de vraisemblance).

Le modèle de base est l'équation différentielle stochastique d'Itô :

```
dX(t) = [ α X(t) − β X(t) ln X(t) ] dt + σ X(t) dW(t),   X(0) = x0 > 0
```

où α est le taux de croissance intrinsèque, β le facteur de décélération, σ l'intensité du bruit et W un mouvement brownien.

## Arborescence

```
these-gompertz/
├── README.md                 présentation du dépôt (ce fichier)
├── LICENSE                   licence MIT (code)
├── .gitignore                fichiers exclus de Git (dont data/raw/)
├── requirements.txt          dépendances Python
├── code/
│   ├── R/                    scripts R
│   └── python/               scripts Python
├── data/
│   ├── raw/                  données réelles confidentielles — JAMAIS versionnées
│   └── simulated/            données simulées (partageables)
├── figures/                  figures produites par les scripts
├── docs/                     plan de gestion des données et documentation
└── references/               bibliographie Zotero vérifiée (.bib)
```

## Reproduire les résultats

Toutes les commandes se lancent depuis la racine du dépôt.

**Python** (3.9 ou plus) :

```bash
pip install -r requirements.txt
python code/python/simulation_gompertz.py
```

**R** (4.0 ou plus, R de base uniquement) :

```bash
Rscript code/R/simulation_gompertz.R
```

Chaque script simule 200 trajectoires par le schéma d'Euler-Maruyama (graine `2026`, pas `dt = 0.01`, horizon `T = 30`), compare la moyenne empirique à la moyenne théorique E[X(t)], puis enregistre :

- les données dans `data/simulated/` (CSV, 20 trajectoires + moyennes) ;
- une figure dans `figures/`.

## Données

| Dossier | Contenu | Versionné | Partage |
|---|---|---|---|
| `data/simulated/` | Trajectoires simulées, régénérables par les scripts | Oui | Public (CC BY 4.0) |
| `data/raw/` | Données agrégées fournies par Inwi, sous accord de confidentialité | **Non** (`.gitignore`) | Jamais partagées |

La gestion complète des données (stockage, sauvegarde, aspects juridiques, archivage) est décrite dans le plan de gestion des données (`docs/`).

## Bibliographie

`references/references_these.bib` : 20 références exportées de Zotero. Chaque DOI a été vérifié (existence, métadonnées, absence de rétractation).

## Licences

- **Code** : licence MIT (voir `LICENSE`).
- **Données simulées** : CC BY 4.0.
- **Données d'Inwi** : propriété d'Inwi, non incluses dans ce dépôt.

## Citer ce travail

Tabbiz, L. (2026). *Modélisation stochastique de systèmes gompertziens corrélés* — dépôt de travail de thèse. Université Abdelmalek Essaâdi. https://github.com/laila-tabbiz/these-gompertz

## Contact

Laila Tabbiz — [ORCID](https://orcid.org/0009-0006-4969-4455) · [LinkedIn](https://www.linkedin.com/in/laila-tabbiz) · [ResearchGate](https://www.researchgate.net/profile/Laila-Tabbiz)
