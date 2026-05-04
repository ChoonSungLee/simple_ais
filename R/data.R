# ============================================================
# R/data.R
#
# Gemini 버전과의 차이점:
#   stan_data에 structural 벡터 추가 (1=sPT, 0=non-sPT)
# ============================================================

library(readxl)
library(tidyverse)
library(here)

# here() 경로 확인 (pedicle_bayes/ 가 나와야 함)
cat("=== 프로젝트 루트 경로 ===\n")
cat(here(), "\n\n")

# ----------------------------------------------------------
# 1. 엑셀 파일 로드
# ----------------------------------------------------------
data_path <- here("data", "kbpark_modify.xlsm")

if (!file.exists(data_path)) {
  stop(paste(
    "파일을 찾을 수 없습니다.\n",
    "확인 경로:", data_path, "\n",
    "kbpark_modify.xlsm을 data/ 폴더에 넣어주세요."
  ))
}

raw_df <- read_excel(data_path, skip = 1)
cat("데이터 로드 완료:", nrow(raw_df), "행\n\n")

# ----------------------------------------------------------
# 2. Wide → Long 변환 (Gemini 버전과 동일)
# ----------------------------------------------------------
clean_df <- raw_df %>%
  select(ID, 33:66) %>%
  rename_with(~str_remove(., "\\.\\.\\..*"), -ID) %>%
  pivot_longer(cols = -ID, names_to = "key", values_to = "width") %>%
  mutate(
    Level      = str_extract(key, "[TL]\\d+"),
    Side_Label = ifelse(str_detect(key, "Lt"), "Left", "Right")
  ) %>%
  filter(!is.na(width), width > 0) %>%
  filter(Level %in% c("T1", "T2", "T3", "T4", "T5", "T6")) %>%
  mutate(
    level_idx = as.integer(factor(Level,
                  levels = c("T1", "T2", "T3", "T4", "T5", "T6")))
  )

# ----------------------------------------------------------
# 3. sPT/non-sPT 그룹 정보 합치기
# ----------------------------------------------------------
# ★ 아래 'structural_col'을 엑셀의 실제 컬럼명으로 수정하세요.
# 일단 그대로 실행하면 오류 메시지와 함께 전체 컬럼 목록이
# 출력되므로 그 중에서 찾으시면 됩니다.
structural_col <- "structural"

if (structural_col %in% colnames(raw_df)) {

  group_df <- raw_df %>%
    select(ID, all_of(structural_col)) %>%
    rename(is_structural = all_of(structural_col)) %>%
    mutate(is_structural = as.integer(is_structural))

  clean_df <- clean_df %>%
    left_join(group_df, by = "ID")

  cat("그룹 정보 합치기 완료\n")

} else {

  cat("사용 가능한 컬럼 목록:\n")
  print(colnames(raw_df))
  stop(paste(
    "\n'", structural_col, "' 컬럼을 찾을 수 없습니다.\n",
    "위 목록에서 sPT 여부를 나타내는 컬럼명을 확인하고\n",
    "structural_col <- '실제컬럼명' 으로 수정하세요."
  ))

}

# 그룹 분포 확인
cat("\n=== 환자별 그룹 분포 ===\n")
clean_df %>%
  distinct(ID, is_structural) %>%
  count(is_structural) %>%
  mutate(Group = ifelse(is_structural == 1, "sPT", "non-sPT")) %>%
  print()

# ----------------------------------------------------------
# 4. Stan 데이터 리스트 생성
# ----------------------------------------------------------
stan_data <- list(
  N                      = nrow(clean_df),
  width                  = as.numeric(clean_df$width),
  vertebra_level_numeric = clean_df$level_idx,
  side_numeric           = ifelse(clean_df$Side_Label == "Left", 1, 2),
  structural             = as.integer(clean_df$is_structural)
)

cat("\n=== stan_data 구조 확인 ===\n")
cat("총 관측치(N)     :", stan_data$N, "\n")
cat("sPT 관측치       :", sum(stan_data$structural == 1), "\n")
cat("non-sPT 관측치   :", sum(stan_data$structural == 0), "\n")

cat("\n[완료] stan_data 준비 완료.\n")
cat("다음 단계: main_analysis.R 실행\n")
