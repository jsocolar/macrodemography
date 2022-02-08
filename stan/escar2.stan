// generated with brms 2.16.3
functions {
 /* Return the log probability of a proper conditional autoregressive (CAR) 
  * prior with a sparse representation for the adjacency matrix
  * Full credit to Max Joseph (https://github.com/mbjoseph/CARstan)
  * Args:
  *   phi: Vector containing the CAR parameters for each location
  *   car: Dependence (usually spatial) parameter for the CAR prior
  *   sdcar: Standard deviation parameter for the CAR prior
  *   Nloc: Number of locations
  *   Nedges: Number of edges (adjacency pairs)
  *   Nneigh: Number of neighbors for each location
  *   eigenW: Eigenvalues of D^(-1/2) * W * D^(-1/2)
  *   edges1, edges2: Sparse representation of adjacency matrix
  * Details:
  *   D = Diag(Nneigh)
  * Returns:
  *   Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real car, real sdcar, int Nloc,
                       int Nedges, data vector Nneigh, data vector eigenW,
                       int[] edges1, int[] edges2) {
    real tau;  // precision parameter
    row_vector[Nloc] phit_D;  // phi' * D
    row_vector[Nloc] phit_W;  // phi' * W
    vector[Nloc] ldet;
    tau = inv_square(sdcar);
    phit_D = (phi .* Nneigh)';
    phit_W = rep_row_vector(0, Nloc);
    for (i in 1:Nedges) {
      phit_W[edges1[i]] = phit_W[edges1[i]] + phi[edges2[i]];
      phit_W[edges2[i]] = phit_W[edges2[i]] + phi[edges1[i]];
    }
    for (i in 1:Nloc) {
      ldet[i] = log1m(car * eigenW[i]);
    }
    return 0.5 * (Nloc * log(tau) + sum(ldet) - 
           tau * (phit_D * phi - car * (phit_W * phi)));
  }
 /* Return the log probability of an intrinsic conditional autoregressive 
  * (ICAR) prior with a sparse representation for the adjacency matrix
  * Full credit to Max Joseph (https://github.com/mbjoseph/CARstan)
  * Args:
  *   phi: Vector containing the CAR parameters for each location
  *   sdcar: Standard deviation parameter for the CAR prior
  *   Nloc: Number of locations
  *   Nedges: Number of edges (adjacency pairs)
  *   Nneigh: Number of neighbors for each location
  *   eigenW: Eigenvalues of D^(-1/2) * W * D^(-1/2)
  *   edges1, edges2: Sparse representation of adjacency matrix
  * Details:
  *   D = Diag(Nneigh)
  * Returns:
  *   Log probability density of CAR prior up to additive constant
  */
  real sparse_icar_lpdf(vector phi, real sdcar, int Nloc, 
                        int Nedges, data vector Nneigh, data vector eigenW, 
                        int[] edges1, int[] edges2) {
    real tau;  // precision parameter
    row_vector[Nloc] phit_D;  // phi' * D
    row_vector[Nloc] phit_W;  // phi' * W
    tau = inv_square(sdcar);
    phit_D = (phi .* Nneigh)';
    phit_W = rep_row_vector(0, Nloc);
    for (i in 1:Nedges) {
      phit_W[edges1[i]] = phit_W[edges1[i]] + phi[edges2[i]];
      phit_W[edges2[i]] = phit_W[edges2[i]] + phi[edges1[i]];
    }
    return 0.5 * ((Nloc - 1) * log(tau) - 
           tau * (phit_D * phi - (phit_W * phi)));
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
   vector[N] Y_var;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for the CAR structure
  int<lower=1> Nloc;
  int<lower=1> Jloc[N];
  int<lower=0> Nedges;
  int<lower=1> edges1[Nedges];
  int<lower=1> edges2[Nedges];
  vector[Nloc] Nneigh;
  vector[Nloc] eigenMcar;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sdcar;  // SD of the CAR structure
  vector[Nloc] rcar;
  real<lower=0, upper=1> car;
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
}
 model {
   vector[N] sigma_new = sqrt(sigma^2 + Y_var);
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += rcar[Jloc[n]];
    }
    target += normal_lpdf(Y | mu, sigma_new);
  }
  // priors including constants
  target += student_t_lpdf(Intercept | 3, 0, 2.5);
  target += student_t_lpdf(sdcar | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += sparse_car_lpdf(
    rcar | car, sdcar, Nloc, Nedges, Nneigh, eigenMcar, edges1, edges2
  );
  target += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}

