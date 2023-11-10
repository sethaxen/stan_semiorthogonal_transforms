/**
 * Transform an unconstrained vector to a N x K semiorthogonal matrix Q.
 * 
 * See: https://doi.org/10.1137/0717034 and http://proceedings.mlr.press/v97/nirwan19a.html
 * 
 * @param y unconstrained vector with length computed by `semiorthogonal_reflector_num_params`
 * @param N number of rows of Q
 * @param K number of columns of Q (must be <= N)
 * @param special whether Q is special orthogonal (i.e. det(Q)=+1). Must be in (0,1). Only
 *                checked if N==K.
 * @return Q semiorthogonal matrix
 */
matrix semiorthogonal_reflector_constrain_lp(vector y, int N, int K, int special){
  tuple(matrix[N,K],vector[K]) Q_fact = semiorthogonal_reflector_factors_lp(y, N, K, special);
  return factors_to_Q(Q_fact.1, Q_fact.2);
}

/**
 * Transform an unconstrained vector to a LAPACK-style factors of an N x K semiorthogonal matrix Q.
 * 
 * @param y unconstrained vector with length computed by `semiorthogonal_reflector_num_params`
 * @param N number of rows of Q
 * @param K number of columns of Q (must be <= N)
 * @param special whether Q is special orthogonal (i.e. det(Q)=+1). Must be in (0,1). Only
 *                checked if N==K.
 * @return Q_fact LAPACK-style factors of Q (V, tau)
 */
tuple(matrix,vector) semiorthogonal_reflector_factors_lp(vector y, int N, int K, int special){
  target += -dot_self(y) / 2;  // std_normal_lupdf(y)
  return reflector_factors(y, N, K, special);
}

/**
 * Return the number of unconstrained parameters needed for a N x K semiorthogonal matrix Q.
 * 
 * @param N number of rows of Q
 * @param K number of columns of Q (must be <= N)
 * @param special whether Q is special orthogonal (i.e. det(Q)=+1). Must be in (0,1). Only
 *                checked if N==K.
 * @return nparams Number of unconstrained parameters
 */
int semiorthogonal_reflector_num_params(int N, int K, int special){
  return N * K - (K * (K - 1)) %/% 2 - (N==K && special);
}

// convert unconstrained parameters to LAPACK-style factorized form of Q
// output is effectively `(Eigen::HouseholderQR::matrixQR(), Eigen::HouseholderQR::hCoeffs())`
tuple(matrix,vector) reflector_factors(vector y, int N, int K, int special){
  matrix[N,K] factors = rep_matrix(0, N, K);
  vector[K] tau;
  int Nlower = N - 1;
  int iy = 1;
  int maxcols = K - (special == 1 && N == K);
  for (j in 1:maxcols) {
    real x1 = y[iy];
    if (Nlower > 0) {
      vector[Nlower] xtail = segment(y, iy + 1, Nlower);
      tuple(real,vector[Nlower],real) reflector = get_reflector(x1, xtail);
      tau[j] = reflector.1;
      factors[j+1:N, j] = reflector.2;
    } else { // avoid length 0 segment
      tuple(real,real) reflector = get_reflector(x1);
      tau[j] = reflector.1;
    }
    factors[j, j] = 1; // set R[j,j] to 1 and discard actual R[j,j]
    iy += Nlower + 1;
    Nlower -= 1;
  }
  if (special == 1 && N == K)
    tau[N] = determinant_from_tau(head(diagonal(factors), N-1));
  return (factors, tau);
}

// convert LAPACK-style factorized form of Q to a dense Q matrix
// adapted from LAPACK DORGRQ
matrix factors_to_Q(matrix factors, vector tau){
  int N = rows(factors);
  int K = cols(factors);
  matrix[N,K] Q = rep_matrix(0, N, K);
  int Nlower = N - K + 1;
  for (j in reverse(linspaced_int_array(K, 1, K))) {
    vector[Nlower] v = factors[j:N, j];
    Q[j, j] = 1;
    Q[j:N, j:K] = apply_reflector(tau[j], v, Q[j:N, j:K]);
    Nlower += 1;
  }
  return Q;
}

// compute parameters of reflector H(tau, v, beta)=(I - tau [1; v] [1; v]') so H*[x1; xtail]=beta*e1
// adapted from LAPACK DLARFGP, see also https://doi.org/10.1137/080725763
// ensures the first entry of the resulting vector is positive while avoiding numerical issues
tuple(real,vector,real) get_reflector(real x1, vector xtail){
  real tau;
  vector[num_elements(xtail)] v;
  real beta;
  real xtail_norm = norm2(xtail);
  if (xtail_norm == 0) {
    tau = 2 * (x1 <= 0);
    v = xtail;
    beta = abs(x1);
    return (tau, v, beta);
  }
  real xnorm = hypot(x1, xtail_norm);
  beta = -xnorm;
  // TODO: needs a possibly_rescale function call here to avoid underflow below
  if (x1 < 0)
    beta = -beta;
  real eta;
  if (beta >= 0) {
    eta = x1 - beta;
  } else {
    beta = -beta;
    real gamma = x1 + beta;
    eta = -xtail_norm * (xtail_norm / gamma);
  }
  tau = -eta / beta;
  v = xtail / eta;
  return (tau, v, beta);
}
tuple(real,real) get_reflector(real x1){ // case where v has length 0
  real tau = 2 * (x1 <= 0);
  real beta = abs(x1);
  return (tau, beta);
}

// apply elementary reflector (I - tau v v') to matrix B
matrix apply_reflector(real tau, vector v, matrix B) {
  if (tau == 0) return B;
  return B - (tau * v) * (v' * B);
}

real determinant_from_tau(vector tau){
  int num_negative_det = 0;
  for (j in 1:num_elements(tau)) {
    // determinant of reflector is -1 iff tau[j]!=0 and 1 iff tau[j]==0
    if (tau[j] != 0)
      num_negative_det += 1;
  }
  return 1 - 2 * fmod(num_negative_det, 2);
}
