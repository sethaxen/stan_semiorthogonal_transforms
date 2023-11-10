data {
  int<lower=1> K;
  int<lower=1,upper=K> N;
  int<lower=0,upper=1> special;
}
parameters {
  matrix[K, N] alpha;
}
transformed parameters {
  matrix[K, N] beta = qr_thin_Q(alpha);
  if (special && N == K && determinant(beta) < 0) {
    vector[N] temp = beta[ , 1];
    beta[ , 1] = beta[ , 2];
    beta[ , 2] = temp;
  }
}
model {
  to_vector(alpha) ~ std_normal();
}
generated quantities {
  real d_beta = not_a_number();
  if (N == K) d_beta = determinant(beta);  // should be 1 if special, else +/- 1
  matrix[N, N] I_beta = beta' * beta;  // should be identity
}
