# Inference

OLS/WLS cluster scores are x_i w_i e_i, pre/postmultiplied by (X'WX)^(-1).
The CR1 correction is G/(G-1)*(N-1)/(N-K), with time effects counted in K.
Use Student t(G-1). RA estimates the untreated regression and treated residual
mean jointly: score_i = [D_i(r_i-theta) -(1-D_i) r_i X_i B mean_T(X) N_T]/N_T.
Sum scores within original clusters and square; CR0, normal reference.
These conventions matched official lpdid on the tested fixtures.

Jackknife deletes every original inference cluster, rebuilding everything:
V_J=(G-1)/G sum(theta_-g - mean(theta_-g))^2. Degrees of freedom G-1.
A failed deletion aborts rather than silently reducing G.

RI preserves the entire cohort-count vector and excludes observed/duplicate
assignments. Assignment columns correspond to sorted cluster IDs; zero=never.
Each cluster must have a single cohort. Supplied assignments provide language
equivalence. Python/R can also sample uniquely without replacement, seeded;
R and Python RNG streams are not intended to match. Python enumerates small
spaces when requested reps cover the space. Stata requires supplied assignments.
Exact enumeration in all languages is possible via the supplied full space.

## November 2026 support revision

The Stata delete-cluster calculation now recomputes the estimator on surviving support. Disappearance of an adoption cohort by itself does not invalidate a deletion. Every required deletion must yield an estimate for a jackknife standard error.
