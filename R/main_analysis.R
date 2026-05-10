# ============================================================
# R/main_analysis.R
#
# Gemini 버전과의 차이점:
#   - mu_sPT / mu_nonsPT 그룹별 결과 확인
#   - mu_diff (그룹 간 차이) 사후분포 시각화
#   - 그룹별 P(width < 2mm) 히트맵 추가
# ============================================================

library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ----------------------------------------------------------
# 1. 레시피 호출 (데이터 준비)
# ----------------------------------------------------------
source(here("R", "data.R"))

# ----------------------------------------------------------
# 2. Stan 모델 실행
# ----------------------------------------------------------
stan_file <- here("stan", "simple_model.stan")

if (!file.exists(stan_file)) {
  stop(paste("Stan 파일을 찾을 수 없습니다.\n확인 경로:", stan_file))
}

cat("Stan 실행 중 (약 5~15분 소요)...\n\n")

fit <- stan(
  file   = stan_file,
  data   = stan_data,
  iter   = 2000,
  chains = 4,
  seed   = 42
)

# 결과 저장 (재실행 없이 불러올 수 있도록)
# output 폴더가 없으면 자동 생성
dir.create(here("output"), showWarnings = FALSE, recursive = TRUE)

saveRDS(fit, here("output", "fit_simple.rds"))
cat("fit 저장 완료 → output/fit_simple.rds\n\n")

# ----------------------------------------------------------
# 3. 기본 결과 확인 (Gemini 버전과 동일)
# ----------------------------------------------------------
summary_fit <- summary(fit)$summary

# mu_sPT만 골라내기
mu_sPT_summary <- summary_fit[grep("mu_sPT\\[", rownames(summary_fit)), ]
sorted_sPT <- mu_sPT_summary[order(mu_sPT_summary[, "mean"]),
                              c("mean", "sd", "2.5%", "97.5%", "Rhat")]

cat("=== sPT 그룹 mu 결과 (오름차순) ===\n")
print(round(sorted_sPT, 3))

# mu_nonsPT만 골라내기
mu_nsPT_summary <- summary_fit[grep("mu_nonsPT\\[", rownames(summary_fit)), ]
sorted_nsPT <- mu_nsPT_summary[order(mu_nsPT_summary[, "mean"]),
                                c("mean", "sd", "2.5%", "97.5%", "Rhat")]

cat("\n=== non-sPT 그룹 mu 결과 (오름차순) ===\n")
print(round(sorted_nsPT, 3))

# mu_diff 확인 (sPT - non-sPT, 음수 = sPT가 더 좁음)
mu_diff_summary <- summary_fit[grep("mu_diff\\[", rownames(summary_fit)), ]
cat("\n=== mu_diff (sPT - non-sPT) — 음수 = sPT가 더 좁음 ===\n")
print(round(mu_diff_summary[, c("mean", "sd", "2.5%", "97.5%")], 3))

# ----------------------------------------------------------
# 4. 시각화
# ----------------------------------------------------------
posterior <- as.matrix(fit)
CONCAVE   <- 2   # Right = concave (left-sided PT curve)
CONVEX    <- 1

# 레벨 이름
level_names <- paste0("T", 1:6)

## (A) sPT vs non-sPT Concave μ 비교 그래프
mu_plot_data <- data.frame()
for (l in 1:6) {
  for (grp in c("sPT", "nonsPT")) {
    col  <- paste0("mu_", grp, "[", l, ",", CONCAVE, "]")
    samp <- posterior[, col]
    mu_plot_data <- rbind(mu_plot_data, data.frame(
      Level    = paste0("T", l),
      Group    = grp,
      Mean     = mean(samp),
      CI_low   = quantile(samp, 0.025),
      CI_high  = quantile(samp, 0.975),
      CI80_low = quantile(samp, 0.10),
      CI80_high= quantile(samp, 0.90)
    ))
  }
}

p_mu_compare <- ggplot(
  mu_plot_data,
  aes(x = Level, y = Mean, color = Group, group = Group)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low,   ymax = CI_high),
                width = 0.15, linewidth = 0.7) +
  geom_errorbar(aes(ymin = CI80_low, ymax = CI80_high),
                width = 0, linewidth = 2, alpha = 0.4) +
  geom_hline(yintercept = 2, linetype = "dashed",
             color = "darkred", linewidth = 0.8) +
  annotate("text", x = 0.7, y = 2.12,
           label = "2 mm 안전 기준", color = "darkred",
           size = 3.5, hjust = 0) +
  scale_color_manual(
    values = c("sPT" = "#D7191C", "nonsPT" = "#2C7BB6"),
    labels = c("sPT" = "sPT (구조적)", "nonsPT" = "non-sPT (비구조적)")
  ) +
  labs(
    title    = "Concave μ 사후분포: sPT vs non-sPT",
    subtitle = "굵은 선: 80% CI / 가는 선: 95% CI",
    x = "Vertebral Level", y = "Posterior Mean (mm)",
    color = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(here("output", "mu_compare_concave.png"),
       p_mu_compare, width = 9, height = 6, dpi = 300)
cat("그래프 저장 → output/mu_compare_concave.png\n")

## (B) mu_diff 사후분포 (Concave, sPT - non-sPT)
diff_data <- data.frame()
for (l in 1:6) {
  col  <- paste0("mu_diff[", l, ",", CONCAVE, "]")
  samp <- posterior[, col]
  diff_data <- rbind(diff_data, data.frame(
    Level          = paste0("T", l),
    Diff_Mean      = mean(samp),
    CI_low         = quantile(samp, 0.025),
    CI_high        = quantile(samp, 0.975),
    P_sPT_narrower = mean(samp < 0)
  ))
}

cat("\n=== Concave: sPT가 non-sPT보다 좁을 사후확률 ===\n")
print(diff_data %>%
  mutate(P_pct = paste0(round(P_sPT_narrower * 100, 1), "%")) %>%
  select(Level, Diff_Mean, CI_low, CI_high, P_pct))

p_diff <- ggplot(diff_data,
                 aes(x = Level, y = Diff_Mean,
                     fill = P_sPT_narrower)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.8) +
  geom_text(aes(
    label = paste0("P=", round(P_sPT_narrower * 100, 0), "%"),
    y     = ifelse(Diff_Mean < 0, CI_low - 0.08, CI_high + 0.08)
  ), size = 4, fontface = "bold") +
  scale_fill_gradient(low  = "#FFFFCC", high = "#D7191C",
                      name = "P(sPT < non-sPT)") +
  labs(
    title    = "μ 차이(sPT − non-sPT) 사후분포 — Concave",
    subtitle = "논문 p=0.002(T3), p=0.003(T4)에 대응하는 베이즈 결과\n음수 = sPT가 더 좁음",
    x = "Vertebral Level", y = "μ_sPT − μ_non-sPT (mm)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

ggsave(here("output", "mu_diff_concave.png"),
       p_diff, width = 8, height = 5, dpi = 300)
cat("그래프 저장 → output/mu_diff_concave.png\n")

## (C) P(width < 2mm) 히트맵 — 그룹별
prob_data <- data.frame()
for (l in 1:6) {
  for (grp in c("sPT", "nonsPT")) {
    col  <- paste0("mu_", grp, "[", l, ",", CONCAVE, "]")
    samp <- posterior[, col]
    prob_data <- rbind(prob_data, data.frame(
      Level      = paste0("T", l),
      Group      = grp,
      P_below2mm = mean(samp < 2) * 100
    ))
  }
}

p_heatmap <- prob_data %>%
  mutate(Group = factor(Group,
           levels = c("sPT", "nonsPT"),
           labels = c("sPT (구조적)", "non-sPT (비구조적)"))) %>%
  ggplot(aes(x = Level, y = Group, fill = P_below2mm)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(round(P_below2mm, 0), "%")),
            color = "white", fontface = "bold", size = 6) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#FEE090", high = "#D7191C",
    midpoint = 50, limits = c(0, 100),
    name = "P(μ < 2mm) %"
  ) +
  labs(
    title    = "Concave μ < 2mm 사후확률 — 그룹별",
    subtitle = "빨간색: 나사 삽입 위험. sPT의 T3/T4가 특히 위험",
    x = "Vertebral Level", y = ""
  ) +
  theme_minimal(base_size = 13)

ggsave(here("output", "prob_heatmap.png"),
       p_heatmap, width = 9, height = 4, dpi = 300)
cat("그래프 저장 → output/prob_heatmap.png\n")

## (D) 전체 사후분포 밀도 그래프 (Gemini 버전 유지)
target_names_sPT <- c(
  "T1_Concave_sPT", "T1_Convex_sPT",
  "T2_Concave_sPT", "T2_Convex_sPT",
  "T3_Concave_sPT", "T3_Convex_sPT",
  "T4_Concave_sPT", "T4_Convex_sPT",
  "T5_Concave_sPT", "T5_Convex_sPT",
  "T6_Concave_sPT", "T6_Convex_sPT"
)

posterior_sPT <- posterior[, grep("mu_sPT\\[", colnames(posterior))]
colnames(posterior_sPT) <- target_names_sPT

p_areas_sPT <- mcmc_areas(posterior_sPT,
                           pars = target_names_sPT,
                           prob = 0.95) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
  ggtitle("sPT 그룹 사후분포 (T1-T6)") +
  xlab("Width (mm)") +
  theme_minimal()

ggsave(here("output", "mcmc_areas_sPT.png"),
       p_areas_sPT, width = 8, height = 8, dpi = 300)
cat("그래프 저장 → output/mcmc_areas_sPT.png\n")

## (E) PPC (모델 검증)
y_rep_matrix <- as.matrix(fit, pars = "y_rep")
p_ppc <- ppc_dens_overlay(stan_data$width, y_rep_matrix[1:100, ]) +
  labs(title    = "PPC: Density Overlay",
       subtitle = "파란 선(예측)이 검은 선(실제)을 잘 감싸면 모델 적합",
       x = "Pedicle Width (mm)") +
  theme_minimal(base_size = 12)

ggsave(here("output", "ppc_density.png"),
       p_ppc, width = 8, height = 5, dpi = 300)
cat("그래프 저장 → output/ppc_density.png\n")

cat("\n[완료] 모든 분석 및 그래프 저장 완료.\n")
cat("output/ 폴더를 확인하세요.\n")
