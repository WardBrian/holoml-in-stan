functions {
  /**
   * Return M1 x M2 matrix of 1 values with blocks in corners set to
   * 0, where the upper left is (r x r), the upper right is (r x r-1),
   * the lower left is (r-1 x r), and the lower right is (r-1 x r-1).
   * This corresponds to zeroing out the lowest-frequency portions of
   * an FFT.
   * @param M1 number of rows
   * @param M2 number of cols
   * @param r block dimension
   * @return matrix of 1 values with 0-padded corners
   */
  matrix beamstop_gen(int M1, int M2, int r) {
    matrix[M1, M2] B_cal = rep_matrix(1, M1, M2);
    if (r == 0) {
      return B_cal;
    }
    // upper left
    B_cal[1 : r, 1 : r] = rep_matrix(0, r, r);
    // upper right
    B_cal[1 : r, M2 - r + 2 : M2] = rep_matrix(0, r, r - 1);
    // lower left
    B_cal[M1 - r + 2 : M1, 1 : r] = rep_matrix(0, r - 1, r);
    // lower right
    B_cal[M1 - r + 2 : M1, M2 - r + 2 : M2] = rep_matrix(0, r - 1, r - 1);
    return B_cal;
  }
  
  /**
   * Return the matrix corresponding to the fast Fourier
   * transform of Z after it is padded with zeros to size
   * N by M
   * When N by M is larger than the dimensions of Z,
   * this computes an oversampled FFT.
   *
   * @param Z matrix of values
   * @param N number of rows desired (must be >= rows(Z))
   * @param M number of columns desired (must be >= cols(Z))
   * @return the FFT of Z padded with zeros
   */
  complex_matrix fft2(complex_matrix Z, int N, int M) {
    int r = rows(Z);
    int c = cols(Z);
    if (r > N) {
      reject("N must be at least rows(Z)");
    }
    if (c > M) {
      reject("M must be at least cols(Z)");
    }
    
    complex_matrix[N, M] pad = rep_matrix(0, N, M);
    pad[1 : r, 1 : c] = Z;
    
    return fft2(pad);
  }
}
data {
  int<lower=0> N; // image dimension
  matrix<lower=0, upper=1>[N, N] R; // registration image
  int<lower=0, upper=N> d; // separation between sample and registration image
  int<lower=N> M1; // rows of padded matrices
  int<lower=2 * N + d> M2; // cols of padded matrices
  int<lower=0, upper=M1> r; // beamstop radius. replaces omega1, omega2 in paper
  
  real<lower=0> N_p; // avg number of photons per pixel
  array[M1, M2] int<lower=0> Y_tilde; // observed number of photons
  
  real<lower=0> sigma; // standard deviation of pixel prior.
}
transformed data {
  matrix[M1, M2] B_cal = beamstop_gen(M1, M2, r);
  matrix[d, N] separation = rep_matrix(0, d, N);
}
parameters {
  matrix<lower=0, upper=1>[N, N] X;
}
model {
  // prior - penalizing L2 on adjacent pixels
  for (i in 1 : rows(X) - 1) {
    X[i] ~ normal(X[i + 1], sigma);
  }
  for (j in 1 : cols(X) - 1) {
    X[ : , j] ~ normal(X[ : , j + 1], sigma);
  }
  
  // likelihood
  // object representing specimen and reference together
  matrix[N, 2 * N + d] X0R = append_col(X, append_col(separation, R));
  // signal - squared magnitude of the (oversampled) FFT
  matrix[M1, M2] Y = abs(fft2(X0R, M1, M2)) .^ 2;
  
  real N_p_over_Y_bar = N_p / mean(Y);
  matrix[M1, M2] lambda = N_p_over_Y_bar * Y;
  
  for (m1 in 1 : M1) {
    for (m2 in 1 : M2) {
      if (B_cal[m1, m2]) {
        Y_tilde[m1, m2] ~ poisson(lambda[m1, m2]);
      }
    }
  }
}