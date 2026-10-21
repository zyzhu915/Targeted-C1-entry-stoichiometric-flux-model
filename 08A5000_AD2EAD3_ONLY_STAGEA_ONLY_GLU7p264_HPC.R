# ============================================================
# 08_STAGEA_ONLY_ONE_SCENARIO_HPC.R
#
# HPC production version: run Stage A only, without Stage B.
# Production settings for 2026-07-05: run 6 scenarios on separate nodes; each scenario = 5000 baseline starts + 8 parameter sets x 5000 starts; 100 cores.
# Run only 1 scenario per job; run the 6 scenarios in parallel using 6 jobs or a job array.
# Save a separate RDS for each scenario for subsequent merging or Stage B processing.
#
# Important principles:
#   1) Keep the liquid-volume basis unchanged: liquid_volume_L = 1.0 L
#   2) Keep the feed-volume basis unchanged: feed_volume_L = 0.25 L/d
#   3) Keep the gas conversion factor unchanged: 25.45 mL/mmol
#   4) Keep VFA/formate/lactate dilution factors unchanged: all are 1 in the corrected new 24h file.
#   5) Do not introduce free background_CH4 / background_IC.
#   6) Run only 1 scenario per job: Full or one of the 5 pathway deletion/restriction scenarios.
#
# Outputs:
#   08A5000_8PARAM_GLU7p264_STAGEA_<scenario>_results.xlsx
#   08A5000_8PARAM_GLU7p264_STAGEA_<scenario>_all_candidates_with_par.rds
#   08A5000_8PARAM_GLU7p264_STAGEA_<scenario>_converted_constraints_from_raw24h.xlsx
#
# Note: this script terminates after Stage A; it neither reads the omics capacity file nor enters Stage B.
# ============================================================

# ============================================================
# 0. Load packages
# ============================================================
packages_needed <- c("readxl", "dplyr", "tidyr", "stringr", "openxlsx", "tibble", "purrr", "parallel")
missing_pkgs <- packages_needed[!vapply(packages_needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  stop(
    "缺少 R 包：", paste(missing_pkgs, collapse = ", "),
    "\n请先在登录节点安装到个人 R library，或让超算管理员加载对应 R module。"
  )
}

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(tibble)
library(purrr)
library(parallel)

# ============================================================
# 1. Parameters: modify paths and verified experimental settings only here
# ============================================================
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit[1] == length(args)) stop("参数 ", flag, " 后面缺少取值。")
  args[hit[1] + 1]
}

work_dir <- get_arg("--work_dir", getwd())
work_dir <- normalizePath(work_dir, winslash = "/", mustWork = FALSE)

raw_24h_file <- get_arg("--raw_24h_file", file.path(work_dir, "24h2.xlsx"))
raw_24h_file <- normalizePath(raw_24h_file, winslash = "/", mustWork = FALSE)

raw_24h_sheet <- get_arg("--raw_24h_sheet", "Sheet2")

groups_keep <- c("AD2", "EAD3")

# ------------------------
# The fixed-parameter table you just provided
# ------------------------
liquid_volume_L <- 1.0000       # Use 1 L to convert liquid-phase concentrations from mg/L to mmol.
feed_volume_L <- 0.2500         # Daily feed volume
COD_volume_L <- 1.0000          # Use 1 L to convert COD from mgO2/L to mmol O2.
gas_molar_volume_mL_per_mmol <- 25.4500

VFA_DILUTION_FACTOR <- 1.0000
FORMATE_DILUTION_FACTOR <- 1.0000
LACTATE_DILUTION_FACTOR <- 1.0000

# Production setting: influent glucose concentration = 7.264705882 g/L.
glucose_feed_g_L <- 7.264705882
sodium_acetate_feed_g_L <- 13.0000
NaHCO3_feed_g_L <- 2.0000
NH4Cl_feed_g_L <- 1.1687
KH2PO4_feed_g_L <- 0.2678
tryptone_feed_g_L <- 0.3068

MW_glucose <- 180.156
MW_sodium_acetate_anhydrous <- 82.0343
MW_NaHCO3 <- 84.0066
MW_NH4Cl <- 53.491
MW_KH2PO4 <- 136.086
MW_acetate <- 60.052
MW_propionate <- 74.079
MW_butyrate <- 88.106
MW_valerate <- 102.133
MW_formate <- 46.025
MW_lactate <- 90.078
MW_O2 <- 32.000
MW_C <- 12.011

# ------------------------
# Stage A run size
# ------------------------
# Use 150 for an initial quick run; this production version uses 5000.
n_starts_arg <- as.integer(get_arg("--n_starts", "5000"))
if (!is.finite(n_starts_arg) || n_starts_arg < 1) stop("--n_starts 必须是正整数。")
n_starts_baseline <- n_starts_arg
n_starts_turnover <- n_starts_arg
n_starts_range <- as.integer(get_arg("--n_starts_range", "120"))

# The debugging version runs only one main parameter combination and temporarily skips range analysis.
# This production version sets RUN_MAIN_ONLY=FALSE, retains only weak/medium turnover (8 sets total), and skips range analysis.
RUN_MAIN_ONLY <- FALSE
RUN_RANGE_ANALYSIS <- FALSE

# Production run: each job runs only 1 of the following 6 scenarios.
# Full: complete network; H2 / formate / PFL / CO2-formate cycling are all allowed.
# No_H2_entry: disable H2-mediated C1 entry.
# No_formate_entry: disable formate-mediated C1 entry.
# H2_only: allow only H2 entry; disable formate entry, CO2-to-formate, formate oxidation, FHL, and PFL formate supply.
# No_FDH_cycle: disable CO2/formate interconversion and FHL-like formate conversion.
# No_PFL: disable PFL-derived formate supply.
scenario_choices <- c("Full", "No_H2_entry", "No_formate_entry", "H2_only", "No_FDH_cycle", "No_PFL")
scenario_arg <- get_arg("--scenario", Sys.getenv("SCENARIO", "Full"))
if (!(scenario_arg %in% scenario_choices)) {
  stop("--scenario 必须是以下之一：", paste(scenario_choices, collapse = ", "),
       "\n当前输入为：", scenario_arg)
}
scenarios_to_run <- scenario_arg

maxit_baseline <- 1500
maxit_turnover <- 2500
maxit_range <- 2500

# ------------------------
# QC and metabolite objective: retain the 05G settings
# ------------------------
qc_tol_CH4 <- 3
qc_tol_acetate <- 3
qc_tol_IC <- 5
qc_tol_CO2 <- 5
qc_tol_constraint <- 1e-4

tol_CH4_obj <- 3
tol_acetate_obj <- 3
tol_IC_obj <- 5
tol_CO2_obj <- 5
tol_H2_obj <- 0.02
tol_formate_obj <- 0.05
tol_lactate_obj <- 0.50

w_CH4 <- 1.0
w_acetate <- 1.0
w_IC <- 1.0
w_CO2 <- 1.0
w_H2 <- 0.20
w_formate <- 0.50
w_lactate <- 0.20

lambda_constraint <- 1000
lambda_flux_parsimony <- 1e-5
lambda_glucose_use_soft <- 0.01

# After baseline metabolite fitting, turnover-aware optimization permits a slightly worse metabolite objective.
# This does not leave the model unconstrained; it prevents the minimum-objective solution from suppressing H2/formate turnover.
objective_guard_factor <- 1.20
objective_guard_delta <- 0.50
lambda_objective_guard <- 2000

# ------------------------
# New: Stage A turnover-aware parameters
# ------------------------
turnover_scale_H2 <- 10
turnover_scale_formate <- 10
C1_entry_scale <- 5
acetate_excess_scale <- 5

# Use multiple weight combinations rather than a single combination to avoid cherry-picking results.
# Parameter descriptions:
#   turnover_weight: encourage H2/formate production-consumption turnover
#   C1_weight: encourage H2/formate flux into methanogenic entry pathways
#   acetate_excess_weight: penalize excessive acetate cycling; none means no penalty
if (RUN_MAIN_ONLY) {
  # Run only the main model: 1 parameter set to first check the direction and representative solution.
  param_grid <- tidyr::crossing(
    turnover_level = c("medium"),
    C1_level = c("medium"),
    acetate_level = c("weak")
  )
} else {
  # Production run with 8 parameter sets: remove turnover_strong and retain only weak / medium.
  param_grid <- tidyr::crossing(
    turnover_level = c("weak", "medium"),
    C1_level = c("weak", "medium"),
    acetate_level = c("none", "weak")
  )
}

param_grid <- param_grid %>%
  mutate(
    w_turnover = case_when(
      turnover_level == "weak" ~ 0.05,
      turnover_level == "medium" ~ 0.15,
      turnover_level == "strong" ~ 0.40,
      TRUE ~ 0.15
    ),
    w_C1_entry = case_when(
      C1_level == "weak" ~ 0.05,
      C1_level == "medium" ~ 0.15,
      TRUE ~ 0.05
    ),
    w_acetate_excess = case_when(
      acetate_level == "none" ~ 0.00,
      acetate_level == "weak" ~ 0.10,
      TRUE ~ 0.00
    ),
    acetate_slack_mmol = case_when(
      acetate_level == "none" ~ 1e6,
      acetate_level == "weak" ~ 10,
      TRUE ~ 10
    ),
    param_id = paste0("turnover_", turnover_level, "__C1_", C1_level, "__acetate_", acetate_level)
  )

main_param_id <- "turnover_medium__C1_medium__acetate_weak"

set.seed(123)

# HPC parallel settings:
# Here, n_cores is the number of workers optimizing starts in parallel within one scenario.
# For SLURM, use --ntasks=1 --cpus-per-task=100; read SLURM_CPUS_PER_TASK=100 here.
cores_default <- Sys.getenv("SLURM_CPUS_PER_TASK", "100")
n_cores <- as.integer(get_arg("--cores", cores_default))
if (!is.finite(n_cores) || n_cores < 1) stop("--cores 必须是正整数。")
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)
cat("Parallel cores for starts =", n_cores, "\n")

safe_scenario_name <- gsub("[^A-Za-z0-9_]+", "_", scenario_arg)
output_prefix <- paste0("08A5000_AD2EAD3_ONLY_8PARAM_GLU7p264_STAGEA_", safe_scenario_name)
output_file <- file.path(work_dir, paste0(output_prefix, "_results.xlsx"))
rds_file <- file.path(work_dir, paste0(output_prefix, "_all_candidates_with_par.rds"))
converted_constraints_file <- file.path(work_dir, paste0(output_prefix, "_converted_constraints_from_raw24h.xlsx"))

cat("============================================================\n")
cat("08A5000 8PARAM GLU7p264 Stage A-only one-scenario HPC turnover-aware model\n")
cat("raw_24h_file =", raw_24h_file, "\n")
cat("glucose_feed_g_L =", glucose_feed_g_L, "\n")
cat("liquid_volume_L =", liquid_volume_L, "\n")
cat("feed_volume_L =", feed_volume_L, "\n")
cat("gas_molar_volume_mL_per_mmol =", gas_molar_volume_mL_per_mmol, "\n")
cat("scenario =", scenario_arg, "\n")
cat("groups_keep =", paste(groups_keep, collapse = ","), "\n")
cat("n_starts_baseline =", n_starts_baseline, "\n")
cat("n_starts_turnover =", n_starts_turnover, "\n")
cat("param_grid_n =", nrow(param_grid), "\n")
cat("param_id list =", paste(param_grid$param_id, collapse = "; "), "\n")
cat("============================================================\n\n")

# ============================================================
# 2. Utility functions
# ============================================================
fix_names <- function(df) {
  df <- as.data.frame(df, check.names = FALSE)
  names(df) <- stringr::str_squish(names(df))
  names(df) <- make.unique(names(df), sep = "_dup")
  tibble::as_tibble(df)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))
zero_if_na <- function(x) { x <- safe_num(x); x[is.na(x)] <- 0; x }
cap_value <- function(x, lo, hi) pmin(pmax(x, lo), hi)
scale_abs_resid <- function(pred, obs, tol) (pred - obs) / tol
vec_from_mat <- function(M) as.vector(t(M))
names_from_groups_vars <- function(groups, vars) as.vector(t(outer(groups, vars, paste, sep = "__")))

share_by_sum <- function(x) {
  x <- zero_if_na(x)
  s <- sum(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 1e-12) return(rep(NA_real_, length(x)))
  x / s
}

# ============================================================
# 3. Generate constraints from the latest raw 24h file
#    Note: perform unit conversions only; do not change the 05G model structure.
# ============================================================
read_raw_24h_wide <- function(path, sheet = raw_24h_sheet) {
  raw <- readxl::read_excel(path, sheet = sheet, col_names = FALSE) %>% as.data.frame(check.names = FALSE)
  if (nrow(raw) < 4) stop("raw 24h 表至少需要 4 行：指标行、组名行、0h、24h。")

  header1 <- as.character(unlist(raw[1, ], use.names = FALSE))
  header2 <- as.character(unlist(raw[2, ], use.names = FALSE))

  # Handle merged Excel cells by filling analyte names to the right.
  for (j in seq_along(header1)) {
    if (is.na(header1[j]) || header1[j] == "") header1[j] <- if (j > 1) header1[j - 1] else ""
  }
  header1 <- stringr::str_squish(header1)
  header2 <- stringr::str_squish(header2)

  dat <- raw[-c(1, 2), , drop = FALSE]
  names(dat) <- paste(header1, header2, sep = "__")
  names(dat)[1] <- "Time_h"
  dat <- tibble::as_tibble(dat)
  dat <- dat %>% mutate(Time_h = safe_num(Time_h))
  for (nm in names(dat)[-1]) dat[[nm]] <- safe_num(dat[[nm]])
  dat
}

pick_col <- function(df, analyte_pattern, group) {
  nm <- names(df)
  hit <- nm[stringr::str_detect(nm, regex(analyte_pattern, ignore_case = TRUE)) &
              stringr::str_detect(nm, paste0("__", group, "$"))]
  if (length(hit) == 0) {
    stop("找不到列：", analyte_pattern, " / ", group,
         "。请检查 Sheet2 表头。现有列名：\n", paste(nm, collapse = "\n"))
  }
  hit[1]
}

value_at_time <- function(df, col, time_h) {
  row <- which(abs(df$Time_h - time_h) < 1e-9)[1]
  if (is.na(row)) stop("找不到 Time(h)=", time_h)
  as.numeric(df[[col]][row])
}

mgL_delta_to_mmol <- function(start, end, MW_g_mol, volume_L = liquid_volume_L, factor = 1) {
  # Return the end-minus-start change in mmol; use factor if raw concentrations require a dilution multiplier.
  ((end * factor) - (start * factor)) * volume_L / MW_g_mol
}

mgC_L_delta_to_Cmmol <- function(start, end, volume_L = liquid_volume_L) {
  (end - start) * volume_L / MW_C
}

mgCOD_L_removed_to_mmolO2 <- function(start, end, volume_L = COD_volume_L) {
  pmax((start - end) * volume_L / MW_O2, 0)
}

gas_mL_to_mmol <- function(delta_mL) delta_mL / gas_molar_volume_mL_per_mmol

convert_raw_24h_to_constraints <- function(path, sheet = raw_24h_sheet) {
  dat <- read_raw_24h_wide(path, sheet)

  out <- lapply(groups_keep, function(g) {
    H2_0 <- value_at_time(dat, pick_col(dat, "H2 cumulative", g), 0)
    H2_24 <- value_at_time(dat, pick_col(dat, "H2 cumulative", g), 24)
    CH4_0 <- value_at_time(dat, pick_col(dat, "CH4 cumulative", g), 0)
    CH4_24 <- value_at_time(dat, pick_col(dat, "CH4 cumulative", g), 24)
    CO2_0 <- value_at_time(dat, pick_col(dat, "CO2 cumulative", g), 0)
    CO2_24 <- value_at_time(dat, pick_col(dat, "CO2 cumulative", g), 24)

    ace0 <- value_at_time(dat, pick_col(dat, "乙酸", g), 0)
    ace24 <- value_at_time(dat, pick_col(dat, "乙酸", g), 24)
    pro0 <- value_at_time(dat, pick_col(dat, "丙酸", g), 0)
    pro24 <- value_at_time(dat, pick_col(dat, "丙酸", g), 24)
    ib0 <- value_at_time(dat, pick_col(dat, "异丁酸", g), 0)
    ib24 <- value_at_time(dat, pick_col(dat, "异丁酸", g), 24)
    nb0 <- value_at_time(dat, pick_col(dat, "正丁酸", g), 0)
    nb24 <- value_at_time(dat, pick_col(dat, "正丁酸", g), 24)
    iv0 <- value_at_time(dat, pick_col(dat, "异戊酸", g), 0)
    iv24 <- value_at_time(dat, pick_col(dat, "异戊酸", g), 24)
    nv0 <- value_at_time(dat, pick_col(dat, "正戊酸", g), 0)
    nv24 <- value_at_time(dat, pick_col(dat, "正戊酸", g), 24)
    form0 <- value_at_time(dat, pick_col(dat, "甲酸", g), 0)
    form24 <- value_at_time(dat, pick_col(dat, "甲酸", g), 24)
    IC0 <- value_at_time(dat, pick_col(dat, "IC", g), 0)
    IC24 <- value_at_time(dat, pick_col(dat, "IC", g), 24)
    COD0 <- value_at_time(dat, pick_col(dat, "COD", g), 0)
    COD24 <- value_at_time(dat, pick_col(dat, "COD", g), 24)
    lac0 <- value_at_time(dat, pick_col(dat, "乳酸", g), 0)
    lac24 <- value_at_time(dat, pick_col(dat, "乳酸", g), 24)

    acetate_start_mmol <- ace0 * VFA_DILUTION_FACTOR * liquid_volume_L / MW_acetate
    lactate_start_mmol <- lac0 * LACTATE_DILUTION_FACTOR * liquid_volume_L / MW_lactate

    dacetate_mmol <- mgL_delta_to_mmol(ace0, ace24, MW_acetate, factor = VFA_DILUTION_FACTOR)
    dformate_mmol <- mgL_delta_to_mmol(form0, form24, MW_formate, factor = FORMATE_DILUTION_FACTOR)
    dIC_mmol <- mgC_L_delta_to_Cmmol(IC0, IC24)
    dlactate_mmol <- mgL_delta_to_mmol(lac0, lac24, MW_lactate, factor = LACTATE_DILUTION_FACTOR)

    pro_cons <- pmax(-mgL_delta_to_mmol(pro0, pro24, MW_propionate, factor = VFA_DILUTION_FACTOR), 0)
    ib_cons <- pmax(-mgL_delta_to_mmol(ib0, ib24, MW_butyrate, factor = VFA_DILUTION_FACTOR), 0)
    nb_cons <- pmax(-mgL_delta_to_mmol(nb0, nb24, MW_butyrate, factor = VFA_DILUTION_FACTOR), 0)
    iv_cons <- pmax(-mgL_delta_to_mmol(iv0, iv24, MW_valerate, factor = VFA_DILUTION_FACTOR), 0)
    nv_cons <- pmax(-mgL_delta_to_mmol(nv0, nv24, MW_valerate, factor = VFA_DILUTION_FACTOR), 0)

    # Retain the 05G/06D approximation: VFA consumption supplies acetate / IC / redox equivalents.
    VFA_acetate_supply_fixed <- pro_cons * 1 + (ib_cons + nb_cons) * 2 + (iv_cons + nv_cons) * 2
    VFA_IC_supply_fixed <- pro_cons * 1 + (iv_cons + nv_cons) * 1
    VFA_redox_supply_fixed <- pro_cons * 1 + (ib_cons + nb_cons) * 2 + (iv_cons + nv_cons) * 2
    VFA_formation_Cmmol_fixed <- 0
    VFA_consumption_Cmmol_fixed <- pro_cons * 3 + (ib_cons + nb_cons) * 4 + (iv_cons + nv_cons) * 5

    glucose_input_mmol <- glucose_feed_g_L * feed_volume_L / MW_glucose * 1000
    sodium_acetate_feed_mmol <- sodium_acetate_feed_g_L * feed_volume_L / MW_sodium_acetate_anhydrous * 1000
    NaHCO3_input_mmol <- NaHCO3_feed_g_L * feed_volume_L / MW_NaHCO3 * 1000

    tibble(
      group = g,
      dCH4_mmol = gas_mL_to_mmol(CH4_24 - ifelse(is.na(CH4_0), 0, CH4_0)),
      dCO2_mmol = gas_mL_to_mmol(CO2_24 - ifelse(is.na(CO2_0), 0, CO2_0)),
      dH2_mmol = gas_mL_to_mmol(H2_24 - ifelse(is.na(H2_0), 0, H2_0)),
      dformate_mmol = dformate_mmol,
      dacetate_mmol = dacetate_mmol,
      dIC_mmol = dIC_mmol,
      dlactate_mmol = dlactate_mmol,
      glucose_input_mmol = glucose_input_mmol,
      glucose_CH4eq_mmol = 3 * glucose_input_mmol,
      acetate_start_mmol = acetate_start_mmol,
      lactate_start_mmol = lactate_start_mmol,
      COD_removed_mmolO2 = mgCOD_L_removed_to_mmolO2(COD0, COD24),
      VFA_acetate_supply_fixed = VFA_acetate_supply_fixed,
      VFA_IC_supply_fixed = VFA_IC_supply_fixed,
      VFA_redox_supply_fixed = VFA_redox_supply_fixed,
      VFA_formation_Cmmol_fixed = VFA_formation_Cmmol_fixed,
      VFA_consumption_Cmmol_fixed = VFA_consumption_Cmmol_fixed,
      sodium_acetate_feed_mmol_audit = sodium_acetate_feed_mmol,
      NaHCO3_input_mmol_audit = NaHCO3_input_mmol,
      tryptone_feed_g_per_day_audit = tryptone_feed_g_L * feed_volume_L,
      raw_acetate_0_mgL = ace0,
      raw_acetate_24_mgL = ace24,
      raw_COD_0_mgL = COD0,
      raw_COD_24_mgL = COD24
    )
  }) %>% bind_rows() %>% arrange(factor(group, levels = groups_keep))

  out
}

constraints <- convert_raw_24h_to_constraints(raw_24h_file, raw_24h_sheet)

# Safety check: prevent an incorrect liquid volume from again yielding only 20 mmol of acetate.
acetate_consumed_check <- pmax(-constraints$dacetate_mmol, 0)
if (mean(acetate_consumed_check, na.rm = TRUE) < 30) {
  stop("乙酸净消耗均值 < 30 mmol，疑似液相体积或稀释倍数错误。请先检查 converted constraints。")
}

wb_constraints <- createWorkbook()
addWorksheet(wb_constraints, "02_converted_constraints")
writeData(wb_constraints, "02_converted_constraints", constraints)
addWorksheet(wb_constraints, "00_parameters")
writeData(wb_constraints, "00_parameters", tibble(
  parameter = c("liquid_volume_L", "feed_volume_L", "COD_volume_L", "gas_molar_volume_mL_per_mmol",
                "VFA_DILUTION_FACTOR", "FORMATE_DILUTION_FACTOR", "LACTATE_DILUTION_FACTOR",
                "glucose_feed_g_L", "sodium_acetate_feed_g_L", "NaHCO3_feed_g_L", "tryptone_feed_g_L"),
  value = c(liquid_volume_L, feed_volume_L, COD_volume_L, gas_molar_volume_mL_per_mmol,
            VFA_DILUTION_FACTOR, FORMATE_DILUTION_FACTOR, LACTATE_DILUTION_FACTOR,
            glucose_feed_g_L, sodium_acetate_feed_g_L, NaHCO3_feed_g_L, tryptone_feed_g_L)
))
saveWorkbook(wb_constraints, converted_constraints_file, overwrite = TRUE)

cat("converted constraints 审计：\n")
print(constraints %>% select(group, dCH4_mmol, dCO2_mmol, dacetate_mmol, dIC_mmol, dformate_mmol,
                             glucose_input_mmol, glucose_CH4eq_mmol, acetate_start_mmol, COD_removed_mmolO2,
                             sodium_acetate_feed_mmol_audit, NaHCO3_input_mmol_audit))
cat("\n")

# ============================================================
# 4. Model variables, obs, and bounds: retain the 05G structure
# ============================================================
var_names <- c(
  "v_glucose_to_pyruvate",
  "v_PFL_formate_supply",
  "v_POR_reduced_Fd_supply",
  "v_acetylCoA_to_acetate_pta_ack",
  "v_acetate_to_CH4",
  "v_H2_entry",
  "v_formate_entry",
  "v_CO2_to_formate",
  "v_formate_to_CO2",
  "v_FHL_formate_to_H2_CO2",
  "v_IC_to_gas_CO2_export",
  "v_redox_to_H2",
  "v_lactate_formation",
  "v_lactate_consumption"
)

G <- length(groups_keep)
V <- length(var_names)
idx <- setNames(seq_along(var_names), var_names)
joint_names <- names_from_groups_vars(groups_keep, var_names)

obs <- list(
  dCH4 = constraints$dCH4_mmol,
  dCO2 = constraints$dCO2_mmol,
  dH2 = constraints$dH2_mmol,
  dformate = constraints$dformate_mmol,
  dacetate = constraints$dacetate_mmol,
  dIC = constraints$dIC_mmol,
  dlactate = constraints$dlactate_mmol,
  glucose_input = constraints$glucose_input_mmol,
  glucose_CH4eq = constraints$glucose_CH4eq_mmol,
  acetate_start = constraints$acetate_start_mmol,
  lactate_start = constraints$lactate_start_mmol,
  COD_removed = constraints$COD_removed_mmolO2,
  VFA_acetate = constraints$VFA_acetate_supply_fixed,
  VFA_IC = constraints$VFA_IC_supply_fixed,
  VFA_redox = constraints$VFA_redox_supply_fixed,
  VFA_formation_C = constraints$VFA_formation_Cmmol_fixed,
  VFA_consumption_C = constraints$VFA_consumption_Cmmol_fixed
)

xmat <- function(x) matrix(x, nrow = G, ncol = V, byrow = TRUE, dimnames = list(groups_keep, var_names))

make_group_bounds_numeric <- function(g_idx, scenario = "Full") {
  ch4_scale <- max(abs(obs$dCH4[g_idx]), 1)
  glucose_max <- max(obs$glucose_input[g_idx], 0)
  pyruvate_max <- max(2 * glucose_max + obs$lactate_start[g_idx] + 5, 1)
  acetate_avail_max <- max(obs$acetate_start[g_idx] + obs$VFA_acetate[g_idx] + pyruvate_max + 5, 1)
  redox_max <- max(pyruvate_max + obs$VFA_redox[g_idx] + 5, 1)
  formate_cycle_max <- max(4 * ch4_scale + abs(obs$dIC[g_idx]) + 10, 10)

  lb <- rep(0, V); names(lb) <- var_names
  ub <- c(
    v_glucose_to_pyruvate = glucose_max,
    v_PFL_formate_supply = pyruvate_max,
    v_POR_reduced_Fd_supply = pyruvate_max,
    v_acetylCoA_to_acetate_pta_ack = pyruvate_max + obs$VFA_acetate[g_idx] + 5,
    v_acetate_to_CH4 = min(acetate_avail_max, ch4_scale + 10),
    v_H2_entry = ch4_scale + 5,
    v_formate_entry = ch4_scale + 5,
    v_CO2_to_formate = formate_cycle_max,
    v_formate_to_CO2 = formate_cycle_max,
    v_FHL_formate_to_H2_CO2 = formate_cycle_max,
    v_IC_to_gas_CO2_export = max(obs$dCO2[g_idx] + abs(obs$dIC[g_idx]) + 10, 10),
    v_redox_to_H2 = redox_max,
    v_lactate_formation = pyruvate_max,
    v_lactate_consumption = max(obs$lactate_start[g_idx] + 5, 5)
  )
  # ===== Disable specified pathways according to the network scenario =====
  # Note: change only bounds, not the underlying stoichiometric equations; disabling a pathway sets its flux upper bound to 0.
  if (scenario == "No_H2_entry") {
    ub["v_H2_entry"] <- 0
  } else if (scenario == "No_formate_entry") {
    ub["v_formate_entry"] <- 0
  } else if (scenario == "H2_only") {
    # Allow only H2 as the C1 entry route; disable formate entry, CO2/formate cycling, FHL, and PFL formate supply.
    ub[c("v_formate_entry", "v_CO2_to_formate", "v_formate_to_CO2",
         "v_FHL_formate_to_H2_CO2", "v_PFL_formate_supply")] <- 0
  } else if (scenario == "No_FDH_cycle") {
    # Disable CO2 <-> formate interconversion cycling; retain PFL formate production and formate entry.
    ub[c("v_CO2_to_formate", "v_formate_to_CO2")] <- 0
  } else if (scenario == "No_PFL") {
    ub["v_PFL_formate_supply"] <- 0
  } else if (scenario == "Full") {
    # Complete network; no additional pathways disabled
  } else {
    stop("未知 scenario: ", scenario)
  }

  list(lb = lb, ub = ub)
}

make_joint_bounds <- function(scenario = "Full") {
  lb_mat <- matrix(0, nrow = G, ncol = V, dimnames = list(groups_keep, var_names))
  ub_mat <- matrix(0, nrow = G, ncol = V, dimnames = list(groups_keep, var_names))
  for (i in seq_len(G)) {
    b <- make_group_bounds_numeric(i, scenario)
    lb_mat[i, ] <- b$lb
    ub_mat[i, ] <- b$ub
  }
  lb <- vec_from_mat(lb_mat); ub <- vec_from_mat(ub_mat)
  names(lb) <- joint_names; names(ub) <- joint_names
  list(lb = lb, ub = ub, lb_mat = lb_mat, ub_mat = ub_mat)
}

# ============================================================
# 5. Predictions, loss, and QC: retain the 05G equations and add turnover diagnostics
# ============================================================
predict_numeric <- function(x) {
  X <- xmat(x)
  g2p <- X[, "v_glucose_to_pyruvate"]
  pfl <- X[, "v_PFL_formate_supply"]
  por <- X[, "v_POR_reduced_Fd_supply"]
  pta <- X[, "v_acetylCoA_to_acetate_pta_ack"]
  aceCH4 <- X[, "v_acetate_to_CH4"]
  h2e <- X[, "v_H2_entry"]
  forme <- X[, "v_formate_entry"]
  co2form <- X[, "v_CO2_to_formate"]
  formco2 <- X[, "v_formate_to_CO2"]
  fhl <- X[, "v_FHL_formate_to_H2_CO2"]
  co2exp <- X[, "v_IC_to_gas_CO2_export"]
  redh2 <- X[, "v_redox_to_H2"]
  lacf <- X[, "v_lactate_formation"]
  lacc <- X[, "v_lactate_consumption"]

  pyruvate_source <- 2 * g2p + lacc
  pyruvate_sink <- pfl + por + lacf + obs$VFA_formation_C / 3
  acetylCoA_source <- pfl + por
  acetylCoA_sink <- pta
  acetate_available <- obs$acetate_start + pta + obs$VFA_acetate
  redox_supply <- por + obs$VFA_redox
  redox_sink <- redh2 + co2form

  pred_CH4 <- aceCH4 + h2e + forme
  pred_H2 <- redh2 + fhl - 4 * h2e
  pred_formate <- pfl + co2form - formco2 - 4 * forme - fhl
  pred_acetate <- pta + obs$VFA_acetate - aceCH4
  pred_IC <- por + obs$VFA_IC + aceCH4 + 3 * forme + formco2 + fhl - co2form - h2e - co2exp
  pred_CO2 <- co2exp
  pred_lactate <- lacf - lacc

  COD_used_from_CH4 <- 2 * pred_CH4
  residual_COD_unexplained <- obs$COD_removed - COD_used_from_CH4
  COD_excess <- pmax(COD_used_from_CH4 - obs$COD_removed, 0)

  constraint_pyruvate_excess <- pmax(pyruvate_sink - pyruvate_source, 0)
  constraint_acetylCoA_excess <- pmax(acetylCoA_sink - acetylCoA_source, 0)
  constraint_acetate_excess <- pmax(aceCH4 - acetate_available, 0)
  constraint_redox_excess <- pmax(redox_sink - redox_supply, 0)
  glucose_unutilized <- obs$glucose_input - g2p
  glucose_CH4eq_used <- 3 * g2p

  H2_production <- redh2 + fhl
  H2_consumption <- 4 * h2e
  H2_turnover_min <- pmin(H2_production, H2_consumption)

  formate_production <- pfl + co2form
  formate_consumption <- 4 * forme + formco2 + fhl
  formate_turnover_min <- pmin(formate_production, formate_consumption)

  C1_total_entry <- h2e + forme
  total_CH4_model <- pred_CH4
  acetate_CH4_fraction_total <- ifelse(total_CH4_model > 1e-12, aceCH4 / total_CH4_model, NA_real_)
  H2_CH4_fraction_total <- ifelse(total_CH4_model > 1e-12, h2e / total_CH4_model, NA_real_)
  formate_CH4_fraction_total <- ifelse(total_CH4_model > 1e-12, forme / total_CH4_model, NA_real_)
  C1_CH4_fraction_total <- ifelse(total_CH4_model > 1e-12, C1_total_entry / total_CH4_model, NA_real_)
  formate_fraction_within_C1 <- ifelse(C1_total_entry > 1e-12, forme / C1_total_entry, NA_real_)
  H2_fraction_within_C1 <- ifelse(C1_total_entry > 1e-12, h2e / C1_total_entry, NA_real_)
  acetate_net_consumption <- pmax(-obs$dacetate, 0)

  list(
    X = X,
    pred_CH4 = pred_CH4, pred_CO2 = pred_CO2, pred_H2 = pred_H2, pred_formate = pred_formate,
    pred_acetate = pred_acetate, pred_IC = pred_IC, pred_lactate = pred_lactate,
    COD_used_from_CH4 = COD_used_from_CH4, residual_COD_unexplained = residual_COD_unexplained, COD_excess = COD_excess,
    pyruvate_source = pyruvate_source, pyruvate_sink = pyruvate_sink,
    acetylCoA_source = acetylCoA_source, acetylCoA_sink = acetylCoA_sink,
    acetate_available = acetate_available, redox_supply = redox_supply, redox_sink = redox_sink,
    constraint_pyruvate_excess = constraint_pyruvate_excess,
    constraint_acetylCoA_excess = constraint_acetylCoA_excess,
    constraint_acetate_excess = constraint_acetate_excess,
    constraint_redox_excess = constraint_redox_excess,
    glucose_unutilized = glucose_unutilized, glucose_CH4eq_used = glucose_CH4eq_used,
    H2_production = H2_production, H2_consumption = H2_consumption, H2_turnover_min = H2_turnover_min,
    formate_production = formate_production, formate_consumption = formate_consumption, formate_turnover_min = formate_turnover_min,
    C1_total_entry = C1_total_entry,
    acetate_CH4_fraction_total = acetate_CH4_fraction_total,
    H2_CH4_fraction_total = H2_CH4_fraction_total,
    formate_CH4_fraction_total = formate_CH4_fraction_total,
    C1_CH4_fraction_total = C1_CH4_fraction_total,
    formate_fraction_within_C1 = formate_fraction_within_C1,
    H2_fraction_within_C1 = H2_fraction_within_C1,
    acetate_net_consumption = acetate_net_consumption
  )
}

metabolite_loss_fast <- function(x) {
  p <- predict_numeric(x)
  r_CH4 <- scale_abs_resid(p$pred_CH4, obs$dCH4, tol_CH4_obj)
  r_CO2 <- scale_abs_resid(p$pred_CO2, obs$dCO2, tol_CO2_obj)
  r_acetate <- scale_abs_resid(p$pred_acetate, obs$dacetate, tol_acetate_obj)
  r_IC <- scale_abs_resid(p$pred_IC, obs$dIC, tol_IC_obj)
  r_H2 <- scale_abs_resid(p$pred_H2, obs$dH2, tol_H2_obj)
  r_formate <- scale_abs_resid(p$pred_formate, obs$dformate, tol_formate_obj)
  r_lactate <- scale_abs_resid(p$pred_lactate, obs$dlactate, tol_lactate_obj)

  SSE_primary <- sum(w_CH4 * r_CH4^2, w_acetate * r_acetate^2, w_IC * r_IC^2, w_CO2 * r_CO2^2, na.rm = TRUE)
  SSE_secondary <- sum(w_H2 * r_H2^2, w_formate * r_formate^2, w_lactate * r_lactate^2, na.rm = TRUE)
  constraint_penalty <- sum(
    (p$constraint_pyruvate_excess / 1)^2,
    (p$constraint_acetylCoA_excess / 1)^2,
    (p$constraint_acetate_excess / 1)^2,
    (p$constraint_redox_excess / 1)^2,
    (p$COD_excess / 1)^2,
    na.rm = TRUE
  )
  flux_parsimony <- sum((x / pmax(abs(x), 1))^2, na.rm = TRUE)
  glucose_use_soft <- sum((p$glucose_unutilized / pmax(obs$glucose_input, 1))^2, na.rm = TRUE)

  total <- SSE_primary + SSE_secondary +
    lambda_constraint * constraint_penalty +
    lambda_flux_parsimony * flux_parsimony +
    lambda_glucose_use_soft * glucose_use_soft

  list(SSE_primary = SSE_primary, SSE_secondary = SSE_secondary, constraint_penalty = constraint_penalty,
       flux_parsimony = flux_parsimony, glucose_use_soft = glucose_use_soft, total = total)
}

qc_metrics <- function(x) {
  p <- predict_numeric(x)
  abs_CH4 <- abs(p$pred_CH4 - obs$dCH4)
  abs_CO2 <- abs(p$pred_CO2 - obs$dCO2)
  abs_acetate <- abs(p$pred_acetate - obs$dacetate)
  abs_IC <- abs(p$pred_IC - obs$dIC)
  constraint_max <- pmax(p$constraint_pyruvate_excess, p$constraint_acetylCoA_excess,
                         p$constraint_acetate_excess, p$constraint_redox_excess,
                         p$COD_excess, na.rm = TRUE)
  tibble(
    max_abs_CH4 = max(abs_CH4, na.rm = TRUE),
    max_abs_CO2 = max(abs_CO2, na.rm = TRUE),
    max_abs_acetate = max(abs_acetate, na.rm = TRUE),
    max_abs_IC = max(abs_IC, na.rm = TRUE),
    max_constraint_excess = max(constraint_max, na.rm = TRUE),
    qc_primary_pass = max_abs_CH4 <= qc_tol_CH4 &
      max_abs_CO2 <= qc_tol_CO2 &
      max_abs_acetate <= qc_tol_acetate &
      max_abs_IC <= qc_tol_IC &
      max_constraint_excess <= qc_tol_constraint
  )
}

turnover_scores <- function(x, prm) {
  p <- predict_numeric(x)
  H2_turnover_score <- sum(log1p(p$H2_turnover_min / turnover_scale_H2), na.rm = TRUE)
  formate_turnover_score <- sum(log1p(p$formate_turnover_min / turnover_scale_formate), na.rm = TRUE)
  C1_entry_score <- sum(log1p(p$C1_total_entry / C1_entry_scale), na.rm = TRUE)
  acetate_excess <- pmax(p$X[, "v_acetate_to_CH4"] - p$acetate_net_consumption - prm$acetate_slack_mmol, 0)
  acetate_excess_penalty <- sum((acetate_excess / acetate_excess_scale)^2, na.rm = TRUE)
  list(H2_turnover_score = H2_turnover_score,
       formate_turnover_score = formate_turnover_score,
       C1_entry_score = C1_entry_score,
       acetate_excess_penalty = acetate_excess_penalty,
       acetate_excess = acetate_excess)
}

objective_turnover_aware <- function(x, prm, allowed_obj) {
  ml <- metabolite_loss_fast(x)
  sc <- turnover_scores(x, prm)
  guard <- pmax(ml$total - allowed_obj, 0)^2
  total <- ml$total + lambda_objective_guard * guard +
    prm$w_acetate_excess * sc$acetate_excess_penalty -
    prm$w_turnover * (sc$H2_turnover_score + sc$formate_turnover_score) -
    prm$w_C1_entry * sc$C1_entry_score
  if (!is.finite(total)) total <- 1e100
  total
}

# ============================================================
# 6. Starting-point generation: retain the 05G approach and add highTurnover starts
# ============================================================
make_smart_start <- function(bounds, branch = "mixed") {
  lb_mat <- bounds$lb_mat; ub_mat <- bounds$ub_mat
  X <- matrix(0, nrow = G, ncol = V, dimnames = list(groups_keep, var_names))

  for (i in seq_len(G)) {
    ch4 <- max(obs$dCH4[i], 0)
    ace_consumed <- max(-obs$dacetate[i], 0)
    glucose_avail <- obs$glucose_input[i]

    ace0 <- min(ace_consumed, 0.92 * ch4, ub_mat[i, "v_acetate_to_CH4"])
    rem_ch4 <- max(ch4 - ace0, 0)

    if (branch == "h2") {
      h2e0 <- 0.75 * rem_ch4; forme0 <- 0.25 * rem_ch4
    } else if (branch == "formate") {
      h2e0 <- 0.25 * rem_ch4; forme0 <- 0.75 * rem_ch4
    } else if (branch == "highC1") {
      h2e0 <- 0.45 * rem_ch4; forme0 <- 0.55 * rem_ch4
    } else if (branch == "highTurnover") {
      h2e0 <- 0.35 * rem_ch4; forme0 <- 0.65 * rem_ch4
    } else if (branch == "acetate") {
      h2e0 <- 0.10 * rem_ch4; forme0 <- 0.10 * rem_ch4
      ace0 <- min(ch4 - h2e0 - forme0, ub_mat[i, "v_acetate_to_CH4"])
    } else {
      h2e0 <- 0.50 * rem_ch4; forme0 <- 0.50 * rem_ch4
    }

    h2e0 <- min(h2e0, ub_mat[i, "v_H2_entry"])
    forme0 <- min(forme0, ub_mat[i, "v_formate_entry"])

    g2p0 <- min(glucose_avail, ub_mat[i, "v_glucose_to_pyruvate"])
    pyruv0 <- 2 * g2p0
    formate_need <- max(4 * forme0 + obs$dformate[i], 0)
    h2_need <- max(4 * h2e0 - obs$dH2[i], 0)

    if (branch %in% c("formate", "highC1", "highTurnover")) {
      pfl0 <- min(formate_need * 0.8, 0.65 * pyruv0, ub_mat[i, "v_PFL_formate_supply"])
    } else {
      pfl0 <- min(formate_need * 0.3, 0.35 * pyruv0, ub_mat[i, "v_PFL_formate_supply"])
    }
    por0 <- min(max(h2_need - obs$VFA_redox[i], 0), max(pyruv0 - pfl0, 0), ub_mat[i, "v_POR_reduced_Fd_supply"])
    if (!is.finite(por0) || por0 < 0) por0 <- 0

    pta0 <- min(pfl0 + por0, ub_mat[i, "v_acetylCoA_to_acetate_pta_ack"])
    ace0 <- min(max(ace0, max(-obs$dacetate[i] + pta0 + obs$VFA_acetate[i], 0)), ub_mat[i, "v_acetate_to_CH4"])
    redh20 <- min(max(4 * h2e0 - obs$dH2[i], 0), ub_mat[i, "v_redox_to_H2"])
    formate_gap <- max(4 * forme0 - pfl0, 0)
    co2form0 <- min(formate_gap, ub_mat[i, "v_CO2_to_formate"])

    if (branch == "highTurnover") {
      formco20 <- min(0.15 * ub_mat[i, "v_formate_to_CO2"], ub_mat[i, "v_formate_to_CO2"])
      fhl0 <- min(0.08 * ub_mat[i, "v_FHL_formate_to_H2_CO2"], ub_mat[i, "v_FHL_formate_to_H2_CO2"])
    } else {
      formco20 <- 0; fhl0 <- 0
    }
    co2exp0 <- min(max(obs$dCO2[i], 0), ub_mat[i, "v_IC_to_gas_CO2_export"])

    X[i, "v_glucose_to_pyruvate"] <- g2p0
    X[i, "v_PFL_formate_supply"] <- pfl0
    X[i, "v_POR_reduced_Fd_supply"] <- por0
    X[i, "v_acetylCoA_to_acetate_pta_ack"] <- pta0
    X[i, "v_acetate_to_CH4"] <- ace0
    X[i, "v_H2_entry"] <- h2e0
    X[i, "v_formate_entry"] <- forme0
    X[i, "v_CO2_to_formate"] <- co2form0
    X[i, "v_formate_to_CO2"] <- formco20
    X[i, "v_FHL_formate_to_H2_CO2"] <- fhl0
    X[i, "v_IC_to_gas_CO2_export"] <- co2exp0
    X[i, "v_redox_to_H2"] <- redh20
    X[i, "v_lactate_formation"] <- 0
    X[i, "v_lactate_consumption"] <- min(max(-obs$dlactate[i], 0), ub_mat[i, "v_lactate_consumption"])
  }

  vec <- vec_from_mat(X); names(vec) <- joint_names
  cap_value(vec, bounds$lb, bounds$ub)
}

make_random_start <- function(bounds) {
  u <- runif(length(bounds$lb), 0, 1)
  x <- bounds$lb + u * (bounds$ub - bounds$lb)
  names(x) <- joint_names
  x
}

jitter_start <- function(x0, bounds, strength = 0.30) {
  noise <- rnorm(length(x0), mean = 0, sd = strength)
  x <- x0 * exp(noise)
  zero_hit <- which(x0 <= 1e-9 & bounds$ub > 0)
  if (length(zero_hit) > 0) x[zero_hit] <- runif(length(zero_hit), 0, 0.15) * bounds$ub[zero_hit]
  names(x) <- joint_names
  cap_value(x, bounds$lb, bounds$ub)
}

make_starts <- function(bounds, n_starts) {
  branches <- c("acetate", "mixed", "h2", "formate", "highC1", "highTurnover")
  base_starts <- lapply(branches, function(b) make_smart_start(bounds, b))
  starts <- list()
  for (b in seq_along(base_starts)) starts[[length(starts) + 1]] <- base_starts[[b]]
  while (length(starts) < ceiling(0.80 * n_starts)) {
    b <- sample(seq_along(base_starts), 1)
    starts[[length(starts) + 1]] <- jitter_start(base_starts[[b]], bounds, strength = runif(1, 0.15, 0.70))
  }
  while (length(starts) < n_starts) starts[[length(starts) + 1]] <- make_random_start(bounds)
  starts
}

run_one_optim <- function(start, bounds, maxit, fn) {
  tryCatch({
    optim(par = start, fn = fn, method = "L-BFGS-B", lower = bounds$lb, upper = bounds$ub,
          control = list(maxit = maxit))
  }, error = function(e) {
    list(par = start, value = Inf, convergence = 999, message = paste("ERROR:", e$message))
  })
}

make_candidate_record <- function(fit, start_id, stage_label, param_row = NULL, allowed_obj = NA_real_, scenario = "Full") {
  ml <- metabolite_loss_fast(fit$par)
  qc <- qc_metrics(fit$par)
  p <- predict_numeric(fit$par)

  if (is.null(param_row)) {
    prm <- tibble(param_id = stage_label, turnover_level = NA_character_, C1_level = NA_character_, acetate_level = NA_character_,
                  w_turnover = NA_real_, w_C1_entry = NA_real_, w_acetate_excess = NA_real_, acetate_slack_mmol = NA_real_)
    sc <- NULL
  } else {
    prm <- param_row
    sc <- turnover_scores(fit$par, prm)
  }

  tibble(
    scenario = scenario,
    stage = stage_label,
    param_id = prm$param_id,
    turnover_level = prm$turnover_level,
    C1_level = prm$C1_level,
    acetate_level = prm$acetate_level,
    w_turnover = prm$w_turnover,
    w_C1_entry = prm$w_C1_entry,
    w_acetate_excess = prm$w_acetate_excess,
    acetate_slack_mmol = prm$acetate_slack_mmol,
    allowed_obj = allowed_obj,
    start_id = start_id,
    value = fit$value,
    convergence = fit$convergence,
    message = ifelse(is.null(fit$message), "", as.character(fit$message)),
    total_metabolite_objective = ml$total,
    SSE_primary = ml$SSE_primary,
    SSE_secondary = ml$SSE_secondary,
    constraint_penalty = ml$constraint_penalty,
    flux_parsimony = ml$flux_parsimony,
    glucose_use_soft = ml$glucose_use_soft,
    max_abs_CH4 = qc$max_abs_CH4,
    max_abs_CO2 = qc$max_abs_CO2,
    max_abs_acetate = qc$max_abs_acetate,
    max_abs_IC = qc$max_abs_IC,
    max_constraint_excess = qc$max_constraint_excess,
    qc_primary_pass = qc$qc_primary_pass,
    H2_turnover_min_sum = sum(p$H2_turnover_min, na.rm = TRUE),
    formate_turnover_min_sum = sum(p$formate_turnover_min, na.rm = TRUE),
    C1_total_entry_sum = sum(p$C1_total_entry, na.rm = TRUE),
    H2_entry_sum = sum(p$X[, "v_H2_entry"], na.rm = TRUE),
    formate_entry_sum = sum(p$X[, "v_formate_entry"], na.rm = TRUE),
    acetate_to_CH4_sum = sum(p$X[, "v_acetate_to_CH4"], na.rm = TRUE),
    CH4_pred_sum = sum(p$pred_CH4, na.rm = TRUE),
    C1_CH4_fraction_sum = ifelse(sum(p$pred_CH4, na.rm = TRUE) > 1e-12,
                                 sum(p$C1_total_entry, na.rm = TRUE) / sum(p$pred_CH4, na.rm = TRUE), NA_real_),
    acetate_CH4_fraction_sum = ifelse(sum(p$pred_CH4, na.rm = TRUE) > 1e-12,
                                      sum(p$X[, "v_acetate_to_CH4"], na.rm = TRUE) / sum(p$pred_CH4, na.rm = TRUE), NA_real_),
    EAD3_formate_fraction_within_C1 = p$formate_fraction_within_C1[which(groups_keep == "EAD3")],
    turnover_score = ifelse(is.null(sc), NA_real_, sc$H2_turnover_score + sc$formate_turnover_score),
    C1_entry_score = ifelse(is.null(sc), NA_real_, sc$C1_entry_score),
    acetate_excess_penalty = ifelse(is.null(sc), NA_real_, sc$acetate_excess_penalty),
    par = list(fit$par)
  )
}

run_baseline_search <- function(bounds, scenario = "Full") {
  starts <- make_starts(bounds, n_starts_baseline)
  cat("\nBaseline Stage A: metabolite-only, starts =", length(starts), " cores =", n_cores, "\n")

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(tibble)
    library(purrr)
  })

  parallel::clusterExport(
    cl,
    varlist = setdiff(ls(envir = .GlobalEnv), "cl"),
    envir = .GlobalEnv
  )

  recs <- parallel::parLapply(
    cl,
    seq_along(starts),
    function(k) {
      fit <- run_one_optim(
        starts[[k]],
        bounds,
        maxit_baseline,
        function(z) metabolite_loss_fast(z)$total
      )
      make_candidate_record(fit, k, "baseline_metabolite_only", scenario = scenario)
    }
  )

  out <- bind_rows(recs) %>% arrange(total_metabolite_objective)

  cat("  done baseline: best_obj=", signif(min(out$total_metabolite_objective, na.rm = TRUE), 5),
      " QC_pass=", sum(out$qc_primary_pass, na.rm = TRUE), "/", nrow(out), "\n")
  flush.console()

  out
}

run_turnover_search_one_param <- function(bounds, prm, allowed_obj, scenario = "Full") {
  starts <- make_starts(bounds, n_starts_turnover)
  cat("\nTurnover-aware Stage A:", prm$param_id, " starts =", length(starts), " cores =", n_cores, "\n")

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(tibble)
    library(purrr)
  })

  parallel::clusterExport(
    cl,
    varlist = setdiff(ls(envir = .GlobalEnv), "cl"),
    envir = .GlobalEnv
  )

  recs <- parallel::parLapply(
    cl,
    seq_along(starts),
    function(k) {
      fit <- run_one_optim(
        starts[[k]],
        bounds,
        maxit_turnover,
        function(z) objective_turnover_aware(z, prm, allowed_obj)
      )
      make_candidate_record(fit, k, "turnover_aware", prm, allowed_obj, scenario = scenario)
    }
  )

  out <- bind_rows(recs) %>% arrange(value)

  cat("  done turnover:", prm$param_id,
      " best_turnover_obj=", signif(min(out$value, na.rm = TRUE), 5),
      " best_metab_obj=", signif(min(out$total_metabolite_objective, na.rm = TRUE), 5),
      " QC_pass=", sum(out$qc_primary_pass, na.rm = TRUE), "/", nrow(out), "\n")
  flush.console()

  out
}

# ============================================================
# 7. Flux range analysis: assess the H2/formate flux space permitted by the data, rather than selecting a representative solution
# ============================================================
range_target_value <- function(x, target_id) {
  p <- predict_numeric(x)
  X <- p$X
  if (target_id == "H2_entry_sum") return(sum(X[, "v_H2_entry"], na.rm = TRUE))
  if (target_id == "formate_entry_sum") return(sum(X[, "v_formate_entry"], na.rm = TRUE))
  if (target_id == "C1_total_entry_sum") return(sum(p$C1_total_entry, na.rm = TRUE))
  if (target_id == "H2_turnover_min_sum") return(sum(p$H2_turnover_min, na.rm = TRUE))
  if (target_id == "formate_turnover_min_sum") return(sum(p$formate_turnover_min, na.rm = TRUE))
  if (target_id == "EAD3_formate_fraction_within_C1") {
    val <- p$formate_fraction_within_C1[which(groups_keep == "EAD3")]
    ifelse(is.na(val), 0, val)
  } else {
    stop("未知 range target: ", target_id)
  }
}

objective_range <- function(x, allowed_obj, target_id, direction = "max") {
  ml <- metabolite_loss_fast(x)
  guard <- pmax(ml$total - allowed_obj, 0)^2
  target <- range_target_value(x, target_id)
  val <- lambda_objective_guard * guard + ml$total * 0.001
  if (direction == "max") val <- val - target
  if (direction == "min") val <- val + target
  if (!is.finite(val)) val <- 1e100
  val
}

run_range_one <- function(bounds, allowed_obj, target_id, direction) {
  starts <- make_starts(bounds, n_starts_range)
  recs <- list()
  cat("\nRange optimization:", direction, target_id, "starts =", length(starts), "\n")
  for (k in seq_along(starts)) {
    fit <- run_one_optim(starts[[k]], bounds, maxit_range,
                         function(z) objective_range(z, allowed_obj, target_id, direction))
    recs[[k]] <- make_candidate_record(fit, k, paste0("range_", direction, "_", target_id), NULL, allowed_obj)
    recs[[k]]$range_target_id <- target_id
    recs[[k]]$range_direction <- direction
    recs[[k]]$range_target_value <- range_target_value(fit$par, target_id)
  }
  bind_rows(recs) %>%
    filter(total_metabolite_objective <= allowed_obj | qc_primary_pass) %>%
    arrange(if (direction == "max") desc(range_target_value) else range_target_value) %>%
    slice_head(n = 1)
}

# ============================================================
# 8. Flux extraction and edge tables
# ============================================================
extract_flux_one <- function(x) {
  p <- predict_numeric(x)
  X <- p$X
  tibble(
    group = groups_keep,
    v_glucose_to_pyruvate = X[, "v_glucose_to_pyruvate"],
    v_PFL_formate_supply = X[, "v_PFL_formate_supply"],
    v_POR_reduced_Fd_supply = X[, "v_POR_reduced_Fd_supply"],
    v_acetylCoA_to_acetate_pta_ack = X[, "v_acetylCoA_to_acetate_pta_ack"],
    v_acetate_to_CH4 = X[, "v_acetate_to_CH4"],
    v_H2_entry = X[, "v_H2_entry"],
    v_formate_entry = X[, "v_formate_entry"],
    v_CO2_to_formate = X[, "v_CO2_to_formate"],
    v_formate_to_CO2 = X[, "v_formate_to_CO2"],
    v_FHL_formate_to_H2_CO2 = X[, "v_FHL_formate_to_H2_CO2"],
    v_IC_to_gas_CO2_export = X[, "v_IC_to_gas_CO2_export"],
    v_redox_to_H2 = X[, "v_redox_to_H2"],
    v_lactate_formation = X[, "v_lactate_formation"],
    v_lactate_consumption = X[, "v_lactate_consumption"],
    pred_CH4 = p$pred_CH4,
    pred_CO2 = p$pred_CO2,
    pred_H2 = p$pred_H2,
    pred_formate = p$pred_formate,
    pred_acetate = p$pred_acetate,
    pred_IC = p$pred_IC,
    pred_lactate = p$pred_lactate,
    H2_production = p$H2_production,
    H2_consumption = p$H2_consumption,
    H2_turnover_min = p$H2_turnover_min,
    formate_production = p$formate_production,
    formate_consumption = p$formate_consumption,
    formate_turnover_min = p$formate_turnover_min,
    C1_total_entry = p$C1_total_entry,
    acetate_CH4_fraction_total = p$acetate_CH4_fraction_total,
    C1_CH4_fraction_total = p$C1_CH4_fraction_total,
    H2_fraction_within_C1 = p$H2_fraction_within_C1,
    formate_fraction_within_C1 = p$formate_fraction_within_C1,
    COD_used_from_CH4 = p$COD_used_from_CH4,
    residual_COD_unexplained = p$residual_COD_unexplained,
    pyruvate_source = p$pyruvate_source,
    pyruvate_sink = p$pyruvate_sink,
    acetylCoA_source = p$acetylCoA_source,
    acetylCoA_sink = p$acetylCoA_sink,
    acetate_available = p$acetate_available,
    redox_supply = p$redox_supply,
    redox_sink = p$redox_sink
  )
}

extract_flux_table <- function(selected_df) {
  rows <- list()
  if (nrow(selected_df) == 0) return(tibble())
  for (i in seq_len(nrow(selected_df))) {
    tmp <- extract_flux_one(selected_df$par[[i]]) %>%
      mutate(
        stage = selected_df$stage[i],
        param_id = selected_df$param_id[i],
        start_id = selected_df$start_id[i],
        total_metabolite_objective = selected_df$total_metabolite_objective[i],
        qc_primary_pass = selected_df$qc_primary_pass[i],
        .before = 1
      )
    rows[[i]] <- tmp
  }
  bind_rows(rows)
}

make_flux_edges <- function(flux_df) {
  if (nrow(flux_df) == 0) return(tibble())
  rows <- list()
  for (i in seq_len(nrow(flux_df))) {
    r <- flux_df[i, ]
    meta <- r %>% select(stage, param_id, group, start_id, total_metabolite_objective, qc_primary_pass)
    edges <- tibble(
      from = c("glucose", "pyruvate", "pyruvate", "acetyl-CoA", "acetate", "redox", "H2", "CO2/IC", "formate", "formate", "IC", "pyruvate", "lactate"),
      to = c("pyruvate", "formate", "reduced Fd/acetyl-CoA", "acetate", "CH4", "H2", "CH4", "formate", "CH4", "CO2/H2", "gas CO2", "lactate", "pyruvate"),
      reaction = c("glucose_to_pyruvate", "PFL_formate_supply", "POR_reduced_Fd_supply", "pta_ack_acetate", "acetate_to_CH4",
                   "redox_to_H2", "H2_entry", "CO2_to_formate", "formate_entry", "FHL/formate_to_CO2", "IC_to_gas_CO2_export", "lactate_formation", "lactate_consumption"),
      flux = c(r$v_glucose_to_pyruvate, r$v_PFL_formate_supply, r$v_POR_reduced_Fd_supply,
               r$v_acetylCoA_to_acetate_pta_ack, r$v_acetate_to_CH4, r$v_redox_to_H2, r$v_H2_entry,
               r$v_CO2_to_formate, r$v_formate_entry, r$v_FHL_formate_to_H2_CO2 + r$v_formate_to_CO2,
               r$v_IC_to_gas_CO2_export, r$v_lactate_formation, r$v_lactate_consumption)
    )
    rows[[i]] <- bind_cols(meta[rep(1, nrow(edges)), ], edges)
  }
  bind_rows(rows)
}

# ============================================================
# 9. Main run: one scenario x (baseline + 8 parameter combinations)
# ============================================================
all_candidates_by_scenario <- list()
scenario_summary_list <- list()

for (scenario_i in scenarios_to_run) {
  cat("\n============================================================\n")
  cat("开始 Stage A scenario:", scenario_i, "\n")
  cat("============================================================\n")

  bounds <- make_joint_bounds(scenario_i)

  # Run a separate baseline for each scenario and use its own best_baseline_obj to determine allowed_obj.
  baseline_candidates <- run_baseline_search(bounds, scenario_i)
  best_baseline_obj <- min(baseline_candidates$total_metabolite_objective, na.rm = TRUE)
  allowed_obj <- max(best_baseline_obj * objective_guard_factor, best_baseline_obj + objective_guard_delta)

  cat("\nScenario:", scenario_i, " baseline best objective =", best_baseline_obj, "\n")
  cat("Scenario:", scenario_i, " allowed objective =", allowed_obj, "\n")
  cat("Scenario:", scenario_i, " baseline QC pass =", sum(baseline_candidates$qc_primary_pass, na.rm = TRUE), "/", nrow(baseline_candidates), "\n\n")

  if (sum(baseline_candidates$qc_primary_pass, na.rm = TRUE) == 0) {
    warning("Scenario ", scenario_i, " baseline 没有 QC pass。请后续检查 constraints 和 bounds。")
  }

  turnover_candidates_list <- list()
  for (i in seq_len(nrow(param_grid))) {
    prm <- param_grid[i, ]
    turnover_candidates_list[[prm$param_id]] <- run_turnover_search_one_param(bounds, prm, allowed_obj, scenario_i)

    # Save after each parameter combination to avoid losing all results if an overnight run is interrupted.
    scenario_interim <- bind_rows(baseline_candidates, bind_rows(turnover_candidates_list))
    all_interim <- bind_rows(all_candidates_by_scenario, list(scenario_interim))
    interim_file <- file.path(work_dir, paste0("08A5000_AD2EAD3_ONLY_8PARAM_GLU7p264_interim_", scenario_i, "_candidates.rds"))
    saveRDS(scenario_interim, interim_file)
    saveRDS(all_interim, file.path(work_dir, paste0("08A5000_AD2EAD3_ONLY_8PARAM_GLU7p264_interim_", scenario_i, "_ALL_candidates.rds")))

    cat("  interim saved:", interim_file, "\n")
    flush.console()
  }

  turnover_candidates <- bind_rows(turnover_candidates_list)
  scenario_all <- bind_rows(baseline_candidates, turnover_candidates)
  all_candidates_by_scenario[[scenario_i]] <- scenario_all

  scenario_summary_list[[scenario_i]] <- scenario_all %>%
    filter(stage == "turnover_aware") %>%
    group_by(scenario, param_id, turnover_level, C1_level, acetate_level) %>%
    summarise(
      n = n(),
      qc_pass = sum(qc_primary_pass, na.rm = TRUE),
      min_metab_obj = min(total_metabolite_objective, na.rm = TRUE),
      median_metab_obj = median(total_metabolite_objective, na.rm = TRUE),
      max_metab_obj = max(total_metabolite_objective, na.rm = TRUE),
      median_C1_total_entry_sum = median(C1_total_entry_sum, na.rm = TRUE),
      max_C1_total_entry_sum = max(C1_total_entry_sum, na.rm = TRUE),
      median_formate_entry_sum = median(formate_entry_sum, na.rm = TRUE),
      max_formate_entry_sum = max(formate_entry_sum, na.rm = TRUE),
      median_H2_entry_sum = median(H2_entry_sum, na.rm = TRUE),
      max_H2_entry_sum = max(H2_entry_sum, na.rm = TRUE),
      .groups = "drop"
    )

  # Save the complete RDS after each scenario.
  all_candidates <- bind_rows(all_candidates_by_scenario)
  saveRDS(all_candidates, rds_file)
  cat("Scenario 完成并保存本 scenario RDS:", scenario_i, " -> ", rds_file, "\n")
}

all_candidates <- bind_rows(all_candidates_by_scenario)
saveRDS(all_candidates, rds_file)

scenario_summary <- bind_rows(scenario_summary_list)

# The simple selection table is for diagnostics only, not a formal Stage B result.
selected_by_param_scenario <- all_candidates %>%
  filter(stage == "turnover_aware", qc_primary_pass) %>%
  group_by(scenario, param_id) %>%
  arrange(total_metabolite_objective, value, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

# ============================================================
# 10. Excel output: write summaries only; save complete candidate solutions in RDS
# ============================================================
wb <- createWorkbook()

addWorksheet(wb, "00_parameters")
writeData(wb, "00_parameters", tibble(
  parameter = c("liquid_volume_L", "feed_volume_L", "COD_volume_L", "gas_molar_volume_mL_per_mmol",
                "glucose_feed_g_L", "n_starts_baseline", "n_starts_turnover",
                "n_cores", "scenarios_to_run",
                "RUN_MAIN_ONLY", "RUN_RANGE_ANALYSIS",
                "objective_guard_factor", "objective_guard_delta"),
  value = c(liquid_volume_L, feed_volume_L, COD_volume_L, gas_molar_volume_mL_per_mmol,
            glucose_feed_g_L, n_starts_baseline, n_starts_turnover,
            n_cores, paste(scenarios_to_run, collapse = "; "),
            RUN_MAIN_ONLY, RUN_RANGE_ANALYSIS,
            objective_guard_factor, objective_guard_delta)
))

addWorksheet(wb, "01_converted_constraints")
writeData(wb, "01_converted_constraints", constraints)

addWorksheet(wb, "02_param_grid")
writeData(wb, "02_param_grid", param_grid)

addWorksheet(wb, "03_scenario_param_summary")
writeData(wb, "03_scenario_param_summary", scenario_summary)

addWorksheet(wb, "04_selected_diagnostic")
writeData(wb, "04_selected_diagnostic", selected_by_param_scenario %>% select(-par))

# With many candidates, save only lightweight tables without par in Excel; see RDS for complete par values.
addWorksheet(wb, "05_all_candidates_no_par")
writeData(wb, "05_all_candidates_no_par", all_candidates %>% select(-par))

saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\n============================================================\n")
cat("Stage A 完成：", output_file, "\n")
cat("Stage A RDS：", rds_file, "\n")
cat("Constraints：", converted_constraints_file, "\n")
cat("说明：完整通量 par 保存在本 scenario 的 RDS；本脚本不进入 Stage B。\n")
cat("============================================================\n")
