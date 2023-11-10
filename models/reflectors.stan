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
  matrix[K, N] beta = semiorthogonal_reflector_constrain_lp(alpha, K, N, special);
}
generated quantities {
  real d_beta = not_a_number();
  if (N == K) d_beta = determinant(beta);  // should be 1 if special, else +/- 1
  matrix[N, N] I_beta = beta' * beta;  // should be identity
}
