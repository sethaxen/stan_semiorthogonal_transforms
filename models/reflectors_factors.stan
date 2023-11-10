functions {
  #include ../transforms/semiorthogonal_reflector_transform.stan
}
data {
  int<lower=1> K;
  int<lower=1,upper=K> N;
  int<lower=0,upper=1> special;
}
transformed data {
  int M = semiorthogonal_reflector_num_params(K, N, special);
}
parameters {
  vector[M] alpha;
}
transformed parameters {
  tuple(matrix[K,N],vector[N]) beta_factors = semiorthogonal_reflector_factors_lp(alpha, K, N, special);
}
generated quantities {
  matrix[K,N] beta = factors_to_Q(beta_factors.1, beta_factors.2);
}
