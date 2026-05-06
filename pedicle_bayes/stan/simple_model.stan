// ============================================================
// simple_model.stan
//
// Gemini 버전과의 차이점:
//   - mu를 그룹별로 분리: mu_sPT[6,2], mu_nonsPT[6,2]
//   - data 블록에 structural[N] 추가
//   - generated quantities에서 그룹 간 차이(mu_diff) 계산
// ============================================================

data {
  int<lower=0> N;
  array[N] int<lower=1, upper=6> vertebra_level_numeric;
  array[N] int<lower=1, upper=2> side_numeric;  // 1: Left, 2: Right
  array[N] real width;
  array[N] int<lower=0, upper=1> structural;    // 1=sPT, 0=non-sPT
}

parameters {
  matrix[6, 2] mu_sPT;      // sPT 그룹: 6레벨 × 2측면
  matrix[6, 2] mu_nonsPT;   // non-sPT 그룹: 6레벨 × 2측면
  real<lower=0> sigma;
}

model {
  // 사전분포
  to_vector(mu_sPT)    ~ normal(5, 3);
  to_vector(mu_nonsPT) ~ normal(5, 3);
  sigma ~ exponential(1);

  // 가능도
  for (i in 1:N) {
    real mu_group = (structural[i] == 1)
      ? mu_sPT[vertebra_level_numeric[i], side_numeric[i]]
      : mu_nonsPT[vertebra_level_numeric[i], side_numeric[i]];
    width[i] ~ normal(mu_group, sigma);
  }
}

generated quantities {
  // 1. 사후예측 (모델 검증용)
  array[N] real y_rep;
  for (n in 1:N) {
    real mu_group = (structural[n] == 1)
      ? mu_sPT[vertebra_level_numeric[n], side_numeric[n]]
      : mu_nonsPT[vertebra_level_numeric[n], side_numeric[n]];
    y_rep[n] = normal_rng(mu_group, sigma);
  }

  // 2. 그룹별 2mm 미만 확률
  matrix[6, 2] prob_narrow_sPT;
  matrix[6, 2] prob_narrow_nonsPT;
  for (k in 1:6) {
    for (s in 1:2) {
      prob_narrow_sPT[k, s]    = (mu_sPT[k, s]    < 2.0);
      prob_narrow_nonsPT[k, s] = (mu_nonsPT[k, s] < 2.0);
    }
  }

  // 3. 그룹 간 차이 (sPT - non-sPT)
  // 음수 = sPT가 더 좁음 (논문의 핵심 발견)
  matrix[6, 2] mu_diff;
  for (k in 1:6) {
    for (s in 1:2) {
      mu_diff[k, s] = mu_sPT[k, s] - mu_nonsPT[k, s];
    }
  }

  // 4. LOO-CV용 log_lik
  vector[N] log_lik;
  for (n in 1:N) {
    real mu_group = (structural[n] == 1)
      ? mu_sPT[vertebra_level_numeric[n], side_numeric[n]]
      : mu_nonsPT[vertebra_level_numeric[n], side_numeric[n]];
    log_lik[n] = normal_lpdf(width[n] | mu_group, sigma);
  }
}
