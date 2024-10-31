README
================
Mitchell Aaron Schepps, Jeremy Seurat, France Mentre, and Weng Kee
Wong <jeremy.seurat@gmail.com>
2024-10-16

# Manuscript Introduction

We provide code used to support our findings in the manuscript “Design
optimization of longitudinal studies using metaheuristics: application
to lithium pharmacokinetics”.

Our zip file contains code and data for the designs, tbales, and figures
of the optimization of a physician constrained pharmacokinetics
pharmacodynamics nonlinear mixed effects model model of
sustained-release lithium using flexible metaheuristic algorithms. We
provide code used to produce the $\text{D}$-, $\text{D}_{s}$-, and
multiobjective optimal designs based on the Fisher Information Matrix
(FIM). Here, we showcase the usage of the the $\texttt{PFIM}$ R code and
$\texttt{ppso}$ and $\texttt{ecr}$ R packages for single objective
optimization and multiobjective optimization with $\texttt{ecr}$. Plots
are created with $\texttt{ggplot2}$.

# Table 1 code

- The code to compute the D-optimal designs for different amount of time
  points in one group by the simplex algorithm implemented in PFIM can
  be found in the PFIM_analyses.R file.

# Table 2 and 3 code

- The code to compute $\text{D}$-optimal designs from the simplex
  algorithm implemented in PFIM for the one, two, and five group
  scenarios with and without a genetic covariate can also be found in
  the table_PFIM.R file.

- File table_ecr.R computes $\text{D}$- and $\text{D}_s$-optimal designs
  for the one and two group scenarios with and without a genetic
  covariate for the $\texttt{ecr}$ package’s single objective
  optimization functionality.

- File table_ppso.R computes $\text{D}$ and $\text{D}_s$-optimal designs
  for the one, two, and five group scenarios with and without a genetic
  covariate for the $\texttt{ppso}$ package’s single objective
  optimization functionality.

# Figure 2 code

- We analyze the convergence properties of the design points. The 100
  iteration search will show the long term convergence of the optim_pso
  function. The code figure_2_code.R uses the data from ppso_long.log to
  create Figure 2, the convergence of each particle to the respective
  design point in PSO.

# Figure 3 code

- Physicians are interested in best estimating the fixed and random
  effects for the clearance parameters as well as the genetic covariate
  in the sustained-release lithium model. We create a multiobjective
  optimization problem and optimize using the ECR package and ecr()
  function to simultaneously maximize $\text{D}_{s3}$- and
  $\text{D}_{s8}$-defficiencies. The “figure_3_code.R” has code to
  create and analyze the Pareto front for Figure 3.

# Supplemental Appendix Table 1 code

The code for Supplementary Table 1, i.e. the five group designs, can be
found within each of the table_PFIM.R, table_ecr.R, and table_ppso.R
files.

# Supplemental Material Figure 1 Code

- We have a single set of nominal values from the earlier lithium study
  with 17 patients and there were no other sets of nominal values
  available. In our case, physicians were also not willing to use other
  sets of nominal values. We implement a hypercube D-optimal design
  (HCD-optimal design) that optimises a pseudo-Bayesian robust
  criterion. To find such a design, we usede bootstrapped confidence
  intervals for each fixed effect as prior information. The criterion
  uses every combination of the 2.5th and 97.5th percentiles from the
  intervals and with 5 fixed effects, there are $2^5=32$ summands. We
  optimized the criterion and found that the HCD-optimal designs are
  similar to our locally D-optimal designs, suggesting that for our
  problem, the latter designs are relatively robust to misspecifications
  in the nominal values. Web Figure 1 shows the D-optimal design time
  points from each of the 32 designs. Provided R code
  “Supplementary_Figure_1.R” contains code to run the HCD experiment and
  presents the figure from precalculated results found in the data
  “HCD_results.csv”.

# References

- Francke T (2020). ppso: Particle Swarm Optimization and Dynamically
  Dimensioned Search, optionally using parallel computing based on Rmpi.
  R package version 0.9-99991. <https://github.com/TillF/ppso>

- Bossek, J. (2023). ecr: Evolutionary Computation in R. R package
  version 2.1.1. <https://github.com/jakobbossek/ecr2>

- Wickahm, H. (2023) ggplot2: Create Elegant Data Visualisations Using
  the Grammar of Graphics. R package version 3.4.4.
  <https://github.com/cran/ggplot2>
