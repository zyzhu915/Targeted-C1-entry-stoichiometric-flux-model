# ============================================================
# 08B_STAGEB_ON_08A5000_FULL_8PARAM_AD2EAD3ONLY_CUTOFF_GATE25_MAXRANK_PFLMIN_ALLV.R
#
# Purpose:
#   Apply the production v26c-style Stage B logic to the AD2/EAD3-only Stage A Full candidate pool with 8 parameter sets:
#   1) Perform omics-based screening separately for the 8 param_id values;
#   2) Evaluate objective <= 0.50 / 0.40 / 0.30 / 0.20 separately;
#   3) Evaluate capacity mode = multiply / sqrt separately;
#   4) Evaluate penalty mode = mixed / all_bilateral separately;
#   5) Core gate: both H2_entry_penalty and formate_entry_penalty must be in the best 25% of the current candidate pool;
#   6) Retain all candidates passing the core gate, with no fallback;
#   7) Output one sheet per param_id; use only AD2/EAD3 omics capacity shares for Stage B and retain all AD2/EAD3 v fluxes, fractions, penalties, ranks, and final rankings.
#
# Important:
#   EAD3_formate_fraction / EAD3_formate_entry are result annotations only;
#   they do not participate in the core gate or selection ranking.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(openxlsx)
  library(stringr)
  library(purrr)
})

# ---------------- User paths ----------------
# Note: use forward slashes (/) rather than Windows backslashes in R to avoid escape-sequence issues.
# This version uses the AD2/EAD3-only Stage A Full candidate pool; AD1 does not participate in Stage B.
work_dir <- "E:/研究生课题/实验4-产甲烷实验/正式实验2025.5.21/2026.6.24代谢模型/AD2EAD3_ONLY_STAGEA_HPC_scripts_FINAL/AD2EAD3_ONLY_STAGEA_HPC/Full"
output_dir <- work_dir

stageA_rds_file <- file.path(work_dir, "08A5000_AD2EAD3_ONLY_8PARAM_GLU7p264_STAGEA_Full_all_candidates_with_par.rds")

# Converted constraints are used only to re-audit the objective for AD2/EAD3.
# Automatically identify possible filenames to accommodate naming differences between versions.
converted_constraints_candidates <- c(
  file.path(work_dir, "08A5000_AD2EAD3_ONLY_8PARAM_GLU7p264_STAGEA_Full_converted_constraints_from_raw24h.xlsx"),
  file.path(work_dir, "08A5000_AD2EAD3_ONLY_8PARAM_GLU7p264_STAGEA_Full_converted_constraints.xlsx"),
  file.path(work_dir, "08A5000_8PARAM_GLU7p264_STAGEA_Full_converted_constraints_from_raw24h.xlsx")
)
converted_constraints_file <- converted_constraints_candidates[file.exists(converted_constraints_candidates)][1]

capacity_file <- "E:/研究生课题/实验4-产甲烷实验/正式实验2025.5.21/2026.6.24代谢模型/2026.7.21 PFL规则修改/09B_v3_reaction_capacity_from_02G_M00357_PFLmin_for_05A.xlsx"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Prevent openxlsx from failing to create workbookTemp in system temporary paths containing Chinese usernames.
safe_tmp <- file.path(output_dir, "_R_tmp")
if (!dir.exists(safe_tmp)) dir.create(safe_tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = safe_tmp, TMP = safe_tmp, TEMP = safe_tmp)

output_file <- file.path(output_dir, "08B_STAGEB_ON_AD2EAD3_ONLY_STAGEA_FULL_8PARAM_CUTOFF_GATE25_MAXRANK_PFLMIN_results.xlsx")
output_rds_file <- file.path(output_dir, "08B_STAGEB_ON_AD2EAD3_ONLY_STAGEA_FULL_8PARAM_CUTOFF_GATE25_MAXRANK_PFLMIN_all_results.rds")
cat("\n路径检查：\n")
cat("work_dir =", work_dir, "\n")
cat("output_dir =", output_dir, "\n")
cat("Stage A RDS =", stageA_rds_file, "\n")
cat("converted constraints =", converted_constraints_file, "\n")
cat("09B capacity =", capacity_file, "\n\n")

# ---------------- Stage B settings ----------------
objective_cutoffs <- c(0.50, 0.40, 0.30, 0.20)
capacity_modes_to_run <- c("multiply", "sqrt")
penalty_modes_to_run <- c("mixed", "all_bilateral")
core_gate_frac <- 0.25
trend_tolerance <- 0.03

# AD2/EAD3-specific metabolite objective, recalculated from saved Stage A flux candidates.
# These are the same objective weights/tolerances used in Stage A.
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

core_reactions <- c("v_H2_entry", "v_formate_entry")
aux_reactions_for_omics <- c(
  "v_CO2_to_formate",
  "v_formate_to_CO2",
  "v_FHL_formate_to_H2_CO2",
  "v_PFL_formate_supply",
  "v_redox_to_H2",
  "v_acetate_to_CH4",
  "v_acetylCoA_to_acetate_pta_ack",
  "v_IC_to_gas_CO2_export"
)

all_omics_reactions <- unique(c(core_reactions, aux_reactions_for_omics))

groups_keep <- c("AD2", "EAD3")

# ---------------- helper ----------------
objective_cutoff_label <- function(x) {
  paste0("objle", sprintf("%.2f", x) %>% str_replace("\\.", "p"))
}

safe_div <- function(a, b) {
  ifelse(is.na(b) | abs(b) < 1e-12, NA_real_, a / b)
}

rank_min_safe <- function(x, desc = FALSE) {
  if (all(is.na(x))) return(rep(NA_integer_, length(x)))
  if (desc) {
    min_rank(desc(x))
  } else {
    min_rank(x)
  }
}

# conservative ranking: smaller = more conservative
# Support two param_id naming conventions:
#   1) turnover_weak__C1_weak__acetate_weak
#   2) weak_C1_weak_acetate_weak
parse_param_levels <- function(param_id) {
  tibble(param_id = param_id) %>%
    mutate(
      turnover_level = case_when(
        str_detect(param_id, "turnover_weak") | str_detect(param_id, "^weak_") ~ "weak",
        str_detect(param_id, "turnover_medium") | str_detect(param_id, "^medium_") ~ "medium",
        str_detect(param_id, "turnover_strong") | str_detect(param_id, "^strong_") ~ "strong",
        TRUE ~ NA_character_
      ),
      C1_level = case_when(
        str_detect(param_id, "C1_weak") ~ "weak",
        str_detect(param_id, "C1_medium") ~ "medium",
        TRUE ~ NA_character_
      ),
      acetate_level = case_when(
        str_detect(param_id, "acetate_none") ~ "none",
        str_detect(param_id, "acetate_weak") ~ "weak",
        TRUE ~ NA_character_
      ),
      turnover_rank = recode(turnover_level, weak = 1L, medium = 2L, strong = 3L, .default = 99L),
      C1_rank = recode(C1_level, weak = 1L, medium = 2L, .default = 99L),
      acetate_rank = recode(acetate_level, none = 1L, weak = 2L, .default = 99L),
      conservative_rank = turnover_rank * 100 + C1_rank * 10 + acetate_rank
    )
}

short_sheet_name <- function(pid) {
  nm <- pid
  nm <- str_replace(nm, "turnover_", "t")
  nm <- str_replace(nm, "__C1_", "_c")
  nm <- str_replace(nm, "__acetate_", "_a")
  nm <- paste0("P_", nm)
  substr(nm, 1, 31)
}

safe_add_sheet <- function(wb, sheet_name) {
  if (sheet_name %in% names(wb)) removeWorksheet(wb, sheet_name)
  addWorksheet(wb, sheet_name)
}

# ---------------- Read Stage A candidates ----------------
if (!file.exists(stageA_rds_file)) stop("Stage A RDS 不存在：", stageA_rds_file)
all_candidates <- readRDS(stageA_rds_file)

# Read the 08A v1000 Full-only Stage A candidate pool.
# If the RDS lacks a scenario column, set it to Full for all rows.
if (!"scenario" %in% names(all_candidates)) {
  all_candidates$scenario <- "Full"
}

stageA_candidates <- all_candidates %>%
  filter(stage == "turnover_aware") %>%
  filter(scenario == "Full" | is.na(scenario))

# Note: this Stage B step treats AD1 as absent.
# Therefore, no longer use Stage A qc_primary_pass to filter the candidate pool here;
# subsequent objective cutoffs, ranking, and omics penalties are based only on AD2/EAD3.
# Retain qc_primary_pass and total_metabolite_objective only as global audit columns for all three groups.
if (nrow(stageA_candidates) == 0) stop("没有 turnover_aware + Full 的候选解。")

param_info <- bind_rows(lapply(sort(unique(stageA_candidates$param_id)), parse_param_levels))

# ---------------- Read 09B capacity ----------------
if (!file.exists(capacity_file)) stop("09B capacity 文件不存在：", capacity_file)
cap <- readxl::read_excel(capacity_file, sheet = "03_capacity_for_05A") %>%
  as_tibble()

required_cols <- c("capacity_mode", "reaction_id", "group", "reaction_capacity", "capacity_share")
missing_cols <- setdiff(required_cols, names(cap))
if (length(missing_cols) > 0) {
  stop("09B sheet 03_capacity_for_05A 缺少列：", paste(missing_cols, collapse = ", "))
}

cap_long <- cap %>%
  transmute(
    capacity_mode = as.character(capacity_mode),
    reaction_id = as.character(reaction_id),
    group = as.character(group),
    reaction_capacity = suppressWarnings(as.numeric(reaction_capacity)),
    capacity_share_original_all3 = suppressWarnings(as.numeric(capacity_share)),
    used_in_05A_main = ifelse("used_in_05A_main" %in% names(cap), as.character(used_in_05A_main), NA_character_),
    used_in_05A_trend = ifelse("used_in_05A_trend" %in% names(cap), as.character(used_in_05A_trend), NA_character_)
  ) %>%
  filter(
    capacity_mode %in% capacity_modes_to_run,
    group %in% groups_keep,
    reaction_id %in% all_omics_reactions,
    is.finite(reaction_capacity)
  ) %>%
  group_by(capacity_mode, reaction_id) %>%
  mutate(
    # AD2-EAD3 ONLY: treat AD1 as absent.
    # Renormalize omics capacity shares using only AD2/EAD3 reaction_capacity values.
    capacity_sum_pair = sum(reaction_capacity, na.rm = TRUE),
    capacity_share = ifelse(capacity_sum_pair > 1e-12, reaction_capacity / capacity_sum_pair, 0)
  ) %>%
  ungroup()

if (nrow(cap_long) == 0) stop("09B capacity 没有成功映射到目标 reactions。")

# ---------------- Extract fluxes from par ----------------
extract_flux_from_par_only <- function(selected_df) {
  rows <- vector("list", nrow(selected_df))
  for (i in seq_len(nrow(selected_df))) {
    par_vec <- selected_df$par[[i]]
    tmp <- tibble(full_name = names(par_vec), flux = as.numeric(par_vec)) %>%
      separate(full_name, into = c("group", "reaction"), sep = "__", remove = TRUE) %>%
      filter(group %in% groups_keep) %>%
      pivot_wider(names_from = reaction, values_from = flux) %>%
      mutate(
        scenario = selected_df$scenario[i],
        stage = selected_df$stage[i],
        param_id = selected_df$param_id[i],
        start_id = selected_df$start_id[i],
        total_metabolite_objective = selected_df$total_metabolite_objective[i],
        value = selected_df$value[i],
        qc_primary_pass = selected_df$qc_primary_pass[i],
        C1_total_entry = v_H2_entry + v_formate_entry,
        H2_fraction_within_C1 = safe_div(v_H2_entry, C1_total_entry),
        formate_fraction_within_C1 = safe_div(v_formate_entry, C1_total_entry),
        H2_production = v_redox_to_H2 + v_FHL_formate_to_H2_CO2,
        H2_consumption = 4 * v_H2_entry,
        H2_turnover_min = pmin(H2_production, H2_consumption),
        formate_production = v_PFL_formate_supply + v_CO2_to_formate,
        formate_consumption = 4 * v_formate_entry + v_formate_to_CO2 + v_FHL_formate_to_H2_CO2,
        formate_turnover_min = pmin(formate_production, formate_consumption),
        acetate_CH4_fraction_total = v_acetate_to_CH4 / (v_acetate_to_CH4 + v_H2_entry + v_formate_entry),
        C1_CH4_fraction_total = C1_total_entry / (v_acetate_to_CH4 + v_H2_entry + v_formate_entry),
        .before = 1
      )
    rows[[i]] <- tmp
  }
  bind_rows(rows)
}

cat("Extracting flux table from Stage A candidates...\n")
flux_group <- extract_flux_from_par_only(stageA_candidates)

# ============================================================
# AD2/EAD3-specific metabolite objective
# ------------------------------------------------------------
# Stage A in this version is already AD2/EAD3-only.
# Still recalculate the metabolite objective for AD2/EAD3 from the saved candidate flux solutions
# for objective cutoffs and final ranking; do not rerun Stage A.
# Retain the original total_metabolite_objective only as an original Stage A audit column.
# ============================================================
if (length(converted_constraints_file) == 0 || is.na(converted_constraints_file) || !file.exists(converted_constraints_file)) {
  stop("没有找到 converted constraints xlsx。请确认 work_dir 下存在 AD2/EAD3-only 的 converted_constraints_from_raw24h.xlsx，或把 converted_constraints_candidates 中的文件名改成你的实际文件名。")
}
constraints_pair <- readxl::read_excel(converted_constraints_file, sheet = "02_converted_constraints") %>%
  as_tibble() %>%
  filter(group %in% groups_keep) %>%
  transmute(
    group = as.character(group),
    dCH4 = as.numeric(dCH4_mmol),
    dCO2 = as.numeric(dCO2_mmol),
    dH2 = as.numeric(dH2_mmol),
    dformate = as.numeric(dformate_mmol),
    dacetate = as.numeric(dacetate_mmol),
    dIC = as.numeric(dIC_mmol),
    dlactate = as.numeric(dlactate_mmol),
    glucose_input = as.numeric(glucose_input_mmol),
    acetate_start = as.numeric(acetate_start_mmol),
    COD_removed = as.numeric(COD_removed_mmolO2),
    VFA_acetate = as.numeric(VFA_acetate_supply_fixed),
    VFA_IC = as.numeric(VFA_IC_supply_fixed),
    VFA_redox = as.numeric(VFA_redox_supply_fixed),
    VFA_formation_C = as.numeric(VFA_formation_Cmmol_fixed)
  )

rxn_cols_for_parsimony <- c(
  "v_glucose_to_pyruvate", "v_PFL_formate_supply", "v_POR_reduced_Fd_supply",
  "v_acetylCoA_to_acetate_pta_ack", "v_acetate_to_CH4",
  "v_H2_entry", "v_formate_entry", "v_CO2_to_formate", "v_formate_to_CO2",
  "v_FHL_formate_to_H2_CO2", "v_IC_to_gas_CO2_export", "v_redox_to_H2",
  "v_lactate_formation", "v_lactate_consumption"
)

flux_group <- flux_group %>%
  left_join(constraints_pair, by = "group") %>%
  mutate(
    pyruvate_source = 2 * v_glucose_to_pyruvate + v_lactate_consumption,
    pyruvate_sink = v_PFL_formate_supply + v_POR_reduced_Fd_supply + v_lactate_formation + VFA_formation_C / 3,
    acetylCoA_source = v_PFL_formate_supply + v_POR_reduced_Fd_supply,
    acetylCoA_sink = v_acetylCoA_to_acetate_pta_ack,
    acetate_available = acetate_start + v_acetylCoA_to_acetate_pta_ack + VFA_acetate,
    redox_supply = v_POR_reduced_Fd_supply + VFA_redox,
    redox_sink = v_redox_to_H2 + v_CO2_to_formate,

    pred_CH4 = v_acetate_to_CH4 + v_H2_entry + v_formate_entry,
    pred_H2 = v_redox_to_H2 + v_FHL_formate_to_H2_CO2 - 4 * v_H2_entry,
    pred_formate = v_PFL_formate_supply + v_CO2_to_formate - v_formate_to_CO2 - 4 * v_formate_entry - v_FHL_formate_to_H2_CO2,
    pred_acetate = v_acetylCoA_to_acetate_pta_ack + VFA_acetate - v_acetate_to_CH4,
    pred_IC = v_POR_reduced_Fd_supply + VFA_IC + v_acetate_to_CH4 + 3 * v_formate_entry +
      v_formate_to_CO2 + v_FHL_formate_to_H2_CO2 - v_CO2_to_formate - v_H2_entry - v_IC_to_gas_CO2_export,
    pred_CO2 = v_IC_to_gas_CO2_export,
    pred_lactate = v_lactate_formation - v_lactate_consumption,

    COD_used_from_CH4 = 2 * pred_CH4,
    COD_excess = pmax(COD_used_from_CH4 - COD_removed, 0),
    constraint_pyruvate_excess = pmax(pyruvate_sink - pyruvate_source, 0),
    constraint_acetylCoA_excess = pmax(acetylCoA_sink - acetylCoA_source, 0),
    constraint_acetate_excess = pmax(v_acetate_to_CH4 - acetate_available, 0),
    constraint_redox_excess = pmax(redox_sink - redox_supply, 0),
    glucose_unutilized = glucose_input - v_glucose_to_pyruvate,

    r_CH4 = (pred_CH4 - dCH4) / tol_CH4_obj,
    r_CO2 = (pred_CO2 - dCO2) / tol_CO2_obj,
    r_acetate = (pred_acetate - dacetate) / tol_acetate_obj,
    r_IC = (pred_IC - dIC) / tol_IC_obj,
    r_H2 = (pred_H2 - dH2) / tol_H2_obj,
    r_formate = (pred_formate - dformate) / tol_formate_obj,
    r_lactate = (pred_lactate - dlactate) / tol_lactate_obj,

    group_SSE_primary = w_CH4 * r_CH4^2 + w_acetate * r_acetate^2 + w_IC * r_IC^2 + w_CO2 * r_CO2^2,
    group_SSE_secondary = w_H2 * r_H2^2 + w_formate * r_formate^2 + w_lactate * r_lactate^2,
    group_constraint_penalty = (constraint_pyruvate_excess / 1)^2 +
      (constraint_acetylCoA_excess / 1)^2 +
      (constraint_acetate_excess / 1)^2 +
      (constraint_redox_excess / 1)^2 +
      (COD_excess / 1)^2,
    group_glucose_use_soft = (glucose_unutilized / pmax(glucose_input, 1))^2
  ) %>%
  rowwise() %>%
  mutate(
    group_flux_parsimony = sum((c_across(all_of(rxn_cols_for_parsimony)) /
                                 pmax(abs(c_across(all_of(rxn_cols_for_parsimony))), 1))^2, na.rm = TRUE),
    group_metabolite_objective = group_SSE_primary + group_SSE_secondary +
      lambda_constraint * group_constraint_penalty +
      lambda_flux_parsimony * group_flux_parsimony +
      lambda_glucose_use_soft * group_glucose_use_soft
  ) %>%
  ungroup()

# candidate summary (one row per candidate)
candidate_summary <- flux_group %>%
  group_by(scenario, param_id, start_id) %>%
  summarise(
    total_metabolite_objective = first(total_metabolite_objective),
    value = first(value),
    qc_primary_pass = first(qc_primary_pass),
    C1_total_entry_sum = sum(C1_total_entry, na.rm = TRUE),
    H2_entry_sum = sum(v_H2_entry, na.rm = TRUE),
    formate_entry_sum = sum(v_formate_entry, na.rm = TRUE),
    acetate_to_CH4_sum = sum(v_acetate_to_CH4, na.rm = TRUE),
    H2_turnover_sum = sum(H2_turnover_min, na.rm = TRUE),
    formate_turnover_sum = sum(formate_turnover_min, na.rm = TRUE),
    acetate_CH4_fraction_sum = acetate_to_CH4_sum / (acetate_to_CH4_sum + H2_entry_sum + formate_entry_sum),
    C1_CH4_fraction_sum = C1_total_entry_sum / (acetate_to_CH4_sum + H2_entry_sum + formate_entry_sum),
    AD2EAD3_metabolite_objective = sum(group_metabolite_objective, na.rm = TRUE),
    AD2EAD3_SSE_primary = sum(group_SSE_primary, na.rm = TRUE),
    AD2EAD3_SSE_secondary = sum(group_SSE_secondary, na.rm = TRUE),
    AD2EAD3_constraint_penalty = sum(group_constraint_penalty, na.rm = TRUE),
    AD2EAD3_flux_parsimony = sum(group_flux_parsimony, na.rm = TRUE),
    AD2EAD3_glucose_use_soft = sum(group_glucose_use_soft, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(param_info, by = "param_id")

# wide flux columns for AD2/EAD3 key reactions
key_flux_cols <- c(
  # All Stage A reaction fluxes v: retain these AD2/EAD3 columns in every param_id sheet.
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
  "v_lactate_consumption",

  # Diagnostic columns for assessing H2/formate entry, turnover, and methane sources
  "C1_total_entry",
  "H2_fraction_within_C1",
  "formate_fraction_within_C1",
  "H2_production",
  "H2_consumption",
  "H2_turnover_min",
  "formate_production",
  "formate_consumption",
  "formate_turnover_min",
  "acetate_CH4_fraction_total",
  "C1_CH4_fraction_total"
)

flux_wide <- flux_group %>%
  select(scenario, param_id, start_id, group, all_of(key_flux_cols)) %>%
  pivot_wider(
    id_cols = c(scenario, param_id, start_id),
    names_from = group,
    values_from = all_of(key_flux_cols),
    names_glue = "{group}_{.value}"
  )

candidate_base <- candidate_summary %>%
  left_join(flux_wide, by = c("scenario", "param_id", "start_id")) %>%
  mutate(
    # The following are AD2-EAD3 result annotations only; they do not participate in the core gate or selection ranking.
    EAD3_formate_dominant = EAD3_formate_fraction_within_C1 >= 0.70,
    EAD3_formate_entry_higher_than_AD2 = EAD3_v_formate_entry >= AD2_v_formate_entry,
    EAD3_formate_fraction_higher_than_AD2 = EAD3_formate_fraction_within_C1 >= AD2_formate_fraction_within_C1,
    EAD3_H2_entry_lower_than_AD2 = EAD3_v_H2_entry <= AD2_v_H2_entry
  )

# ---------------- penalty calculation ----------------
calc_penalties_for_pool <- function(pool_df, capacity_mode_i, penalty_mode_i) {
  cap_i <- cap_long %>% filter(capacity_mode == capacity_mode_i)
  if (nrow(cap_i) == 0) stop("No capacity rows for mode: ", capacity_mode_i)

  # flux long for candidates in pool and target omics reactions
  flux_long <- flux_group %>%
    semi_join(pool_df %>% select(scenario, param_id, start_id), by = c("scenario", "param_id", "start_id")) %>%
    select(scenario, param_id, start_id, group, all_of(all_omics_reactions)) %>%
    pivot_longer(cols = all_of(all_omics_reactions), names_to = "reaction_id", values_to = "flux") %>%
    mutate(abs_flux = abs(flux)) %>%
    group_by(scenario, param_id, start_id, reaction_id) %>%
    mutate(total_abs_flux = sum(abs_flux, na.rm = TRUE), flux_share = ifelse(total_abs_flux > 1e-12, abs_flux / total_abs_flux, 0)) %>%
    ungroup()

  joined <- flux_long %>%
    left_join(cap_i %>% select(reaction_id, group, capacity_share, reaction_capacity), by = c("reaction_id", "group")) %>%
    mutate(
      capacity_share = ifelse(is.na(capacity_share), 0, capacity_share),
      diff = flux_share - capacity_share,
      # 05G/06B-style: trend_tolerance = 0.03
      # core reactions always bilateral: penalize deviations in both directions beyond tolerance
      # auxiliary reactions in mixed mode use GIMME-like one-sided penalty
      penalty_component = case_when(
        reaction_id %in% core_reactions ~ pmax(abs(diff) - trend_tolerance, 0)^2,
        penalty_mode_i == "all_bilateral" ~ pmax(abs(diff) - trend_tolerance, 0)^2,
        penalty_mode_i == "mixed" ~ pmax(diff - trend_tolerance, 0)^2,
        TRUE ~ pmax(abs(diff) - trend_tolerance, 0)^2
      )
    )

  reaction_penalty <- joined %>%
    group_by(scenario, param_id, start_id, reaction_id) %>%
    summarise(
      reaction_penalty = sum(penalty_component, na.rm = TRUE),
      flux_cap_support = sum(pmin(flux_share, capacity_share), na.rm = TRUE),
      .groups = "drop"
    )

  out <- reaction_penalty %>%
    group_by(scenario, param_id, start_id) %>%
    summarise(
      H2_entry_penalty = reaction_penalty[match("v_H2_entry", reaction_id)] %>% {ifelse(length(.) == 0 | is.na(.), NA_real_, .)},
      formate_entry_penalty = reaction_penalty[match("v_formate_entry", reaction_id)] %>% {ifelse(length(.) == 0 | is.na(.), NA_real_, .)},
      aux_omics_penalty = mean(reaction_penalty[!(reaction_id %in% core_reactions)], na.rm = TRUE),
      total_omics_penalty = mean(reaction_penalty, na.rm = TRUE),
      n_omics_reactions = n_distinct(reaction_id[is.finite(reaction_penalty)]),
      .groups = "drop"
    )

  out
}

# ---------------- score one pool ----------------
score_one_pool <- function(param_i, cutoff_i, cap_i, pen_i) {
  pool <- candidate_base %>%
    filter(param_id == param_i, AD2EAD3_metabolite_objective <= cutoff_i)

  if (nrow(pool) == 0) return(tibble())

  pen <- calc_penalties_for_pool(pool, cap_i, pen_i)

  out <- pool %>%
    left_join(pen, by = c("scenario", "param_id", "start_id")) %>%
    mutate(
      objective_cutoff = cutoff_i,
      objective_cutoff_label = objective_cutoff_label(cutoff_i),
      objective_cutoff_basis = "AD2EAD3_metabolite_objective",
      capacity_mode = cap_i,
      penalty_mode = pen_i,
      near_pool_n = nrow(pool),
      core_gate_threshold_n = ceiling(core_gate_frac * near_pool_n),
      rank_global_metabolite_objective = rank_min_safe(total_metabolite_objective, desc = FALSE),
      rank_AD2EAD3_metabolite_objective = rank_min_safe(AD2EAD3_metabolite_objective, desc = FALSE),
      # Final Stage B ranking uses AD2/EAD3-specific metabolite objective.
      # total_metabolite_objective is retained as global three-reactor feasibility audit.
      rank_metabolite_objective = rank_AD2EAD3_metabolite_objective,
      rank_H2_entry_penalty = rank_min_safe(H2_entry_penalty, desc = FALSE),
      rank_formate_entry_penalty = rank_min_safe(formate_entry_penalty, desc = FALSE),
      worst_core_entry_rank = pmax(rank_H2_entry_penalty, rank_formate_entry_penalty, na.rm = TRUE),
      rank_aux_omics_penalty = rank_min_safe(aux_omics_penalty, desc = FALSE),
      core_gate_pass = rank_H2_entry_penalty <= core_gate_threshold_n & rank_formate_entry_penalty <= core_gate_threshold_n,
      selection_score = rank_metabolite_objective + worst_core_entry_rank + rank_aux_omics_penalty,
      final_rank_within_pool = rank_min_safe(selection_score, desc = FALSE)
    )
  out
}

# ---------------- all settings ----------------
param_ids <- sort(unique(stageA_candidates$param_id))
cat("param_id count =", length(param_ids), "\n")
cat("param_id list =", paste(param_ids, collapse = "; "), "\n\n")
all_results_list <- list()

cat("Scoring param_id × obj(0.50/0.40/0.30/0.20) × capacity × penalty pools...\n")
for (param_i in param_ids) {
  cat("  param:", param_i, "\n")
  for (cutoff_i in objective_cutoffs) {
    for (cap_i in capacity_modes_to_run) {
      for (pen_i in penalty_modes_to_run) {
        key <- paste(param_i, objective_cutoff_label(cutoff_i), cap_i, pen_i, sep = "__")
        all_results_list[[key]] <- score_one_pool(param_i, cutoff_i, cap_i, pen_i)
      }
    }
  }
}

all_scored <- bind_rows(all_results_list)
if (nrow(all_scored) == 0) stop("所有 pool 都为空。")

# candidate audit: one row per param_id before Stage B settings
candidate_audit <- candidate_base %>%
  group_by(scenario, param_id, turnover_level, C1_level, acetate_level, conservative_rank) %>%
  summarise(
    stageA_candidate_n = n(),
    stageA_global_qc_pass_n = sum(qc_primary_pass, na.rm = TRUE),
    n_AD2EAD3_obj_le_0p50 = sum(AD2EAD3_metabolite_objective <= 0.50, na.rm = TRUE),
    n_AD2EAD3_obj_le_0p40 = sum(AD2EAD3_metabolite_objective <= 0.40, na.rm = TRUE),
    n_AD2EAD3_obj_le_0p30 = sum(AD2EAD3_metabolite_objective <= 0.30, na.rm = TRUE),
    n_AD2EAD3_obj_le_0p20 = sum(AD2EAD3_metabolite_objective <= 0.20, na.rm = TRUE),
    min_AD2EAD3_metabolite_objective = min(AD2EAD3_metabolite_objective, na.rm = TRUE),
    median_AD2EAD3_metabolite_objective = median(AD2EAD3_metabolite_objective, na.rm = TRUE),
    min_stageA_total_metabolite_objective_audit = min(total_metabolite_objective, na.rm = TRUE),
    median_stageA_total_metabolite_objective_audit = median(total_metabolite_objective, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(conservative_rank)

gated_all <- all_scored %>% filter(core_gate_pass)

selected_each_setting <- gated_all %>%
  group_by(scenario, param_id, objective_cutoff, objective_cutoff_label, capacity_mode, penalty_mode) %>%
  arrange(selection_score, worst_core_entry_rank, rank_aux_omics_penalty, rank_metabolite_objective,
          total_metabolite_objective, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

setting_summary <- all_scored %>%
  group_by(scenario, param_id, turnover_level, C1_level, acetate_level, conservative_rank,
           objective_cutoff, objective_cutoff_label, capacity_mode, penalty_mode) %>%
  summarise(
    near_pool_n = first(near_pool_n),
    core_gate_threshold_n = first(core_gate_threshold_n),
    core_gate_n = sum(core_gate_pass, na.rm = TRUE),
    core_gate_frac_real = core_gate_n / near_pool_n,
    best_selection_score_gate = ifelse(core_gate_n > 0, min(selection_score[core_gate_pass], na.rm = TRUE), NA_real_),
    median_selection_score_gate = ifelse(core_gate_n > 0, median(selection_score[core_gate_pass], na.rm = TRUE), NA_real_),
    median_H2_entry_penalty_gate = ifelse(core_gate_n > 0, median(H2_entry_penalty[core_gate_pass], na.rm = TRUE), NA_real_),
    median_formate_entry_penalty_gate = ifelse(core_gate_n > 0, median(formate_entry_penalty[core_gate_pass], na.rm = TRUE), NA_real_),
    median_aux_omics_penalty_gate = ifelse(core_gate_n > 0, median(aux_omics_penalty[core_gate_pass], na.rm = TRUE), NA_real_),
    median_C1_entry_gate = ifelse(core_gate_n > 0, median(C1_total_entry_sum[core_gate_pass], na.rm = TRUE), NA_real_),
    median_H2_entry_gate = ifelse(core_gate_n > 0, median(H2_entry_sum[core_gate_pass], na.rm = TRUE), NA_real_),
    median_formate_entry_gate = ifelse(core_gate_n > 0, median(formate_entry_sum[core_gate_pass], na.rm = TRUE), NA_real_),
    n_EAD3_formate_dominant_in_gate = sum(core_gate_pass & EAD3_formate_dominant, na.rm = TRUE),
    n_EAD3_formate_entry_higher_than_AD2_in_gate = sum(core_gate_pass & EAD3_formate_entry_higher_than_AD2, na.rm = TRUE),
    n_EAD3_formate_fraction_higher_than_AD2_in_gate = sum(core_gate_pass & EAD3_formate_fraction_higher_than_AD2, na.rm = TRUE),
    median_EAD3_formate_fraction_gate = ifelse(core_gate_n > 0, median(EAD3_formate_fraction_within_C1[core_gate_pass], na.rm = TRUE), NA_real_),
    max_EAD3_formate_fraction_gate = ifelse(core_gate_n > 0, max(EAD3_formate_fraction_within_C1[core_gate_pass], na.rm = TRUE), NA_real_),
    median_EAD3_formate_entry_gate = ifelse(core_gate_n > 0, median(EAD3_v_formate_entry[core_gate_pass], na.rm = TRUE), NA_real_),
    max_EAD3_formate_entry_gate = ifelse(core_gate_n > 0, max(EAD3_v_formate_entry[core_gate_pass], na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>% arrange(conservative_rank, objective_cutoff, capacity_mode, penalty_mode)

# flux group for selected candidates
selected_flux_group <- flux_group %>%
  inner_join(
    selected_each_setting %>%
      select(scenario, param_id, start_id, objective_cutoff, objective_cutoff_label, objective_cutoff_basis, capacity_mode, penalty_mode,
             selection_score, final_rank_within_pool, rank_metabolite_objective,
             rank_H2_entry_penalty, rank_formate_entry_penalty, worst_core_entry_rank,
             rank_aux_omics_penalty, H2_entry_penalty, formate_entry_penalty,
             aux_omics_penalty, total_omics_penalty, n_omics_reactions),
    by = c("scenario", "param_id", "start_id"),
    relationship = "many-to-many"
  ) %>%
  arrange(scenario, param_id, objective_cutoff, capacity_mode, penalty_mode, start_id, group)

# flux group for all gated candidates (optional long table, useful for plotting)
gated_flux_group <- flux_group %>%
  inner_join(
    gated_all %>%
      select(scenario, param_id, start_id, objective_cutoff, objective_cutoff_label, objective_cutoff_basis, capacity_mode, penalty_mode,
             selection_score, final_rank_within_pool, rank_metabolite_objective,
             rank_H2_entry_penalty, rank_formate_entry_penalty, worst_core_entry_rank,
             rank_aux_omics_penalty, H2_entry_penalty, formate_entry_penalty,
             aux_omics_penalty, total_omics_penalty, n_omics_reactions),
    by = c("scenario", "param_id", "start_id"),
    relationship = "many-to-many"
  ) %>%
  arrange(scenario, param_id, objective_cutoff, capacity_mode, penalty_mode, selection_score, start_id, group)


# Save the complete RDS to avoid Excel row-limit issues or subsequent recomputation.
saveRDS(list(
  all_scored = all_scored,
  gated_all = gated_all,
  selected_each_setting = selected_each_setting,
  setting_summary = setting_summary,
  candidate_audit = candidate_audit,
  selected_flux_group = selected_flux_group,
  gated_flux_group = gated_flux_group,
  candidate_base = candidate_base,
  flux_group = flux_group,
  cap_long = cap_long
), output_rds_file)
cat("完整 Stage B RDS 已保存：", output_rds_file, "\n")

# ---------------- Excel output ----------------
wb <- createWorkbook()

safe_add_sheet(wb, "00_note")
writeData(wb, "00_note", tibble(
  item = c("说明", "Stage A RDS", "Capacity file", "Objective cutoffs", "Capacity modes", "Penalty modes", "Core gate", "Ranking", "EAD3 annotation"),
  value = c(
    "08B: 读取 AD2/EAD3-only Stage A Full scenario RDS；8组 Stage A 参数分别进行 v26c/06B-style Stage B。AD1 在 Stage A 和 Stage B 均不参与；objective cutoff、metabolite ranking、capacity share、flux share、penalty、core gate 和 selection ranking 均只基于 AD2/EAD3。objective cutoff 基于 AD2EAD3_metabolite_objective <= 0.50/0.40/0.30/0.20；所有通过 core gate 的候选解保留；每个参数 sheet 保留 AD2/EAD3 全部 v 通量列。",
    stageA_rds_file,
    capacity_file,
    paste(objective_cutoffs, collapse = ", "),
    paste(capacity_modes_to_run, collapse = ", "),
    paste(penalty_modes_to_run, collapse = ", "),
    "H2_entry_penalty 与 formate_entry_penalty 均位于当前 pool 前 25%。不 fallback。trend_tolerance = 0.03。",
    "selection_score = rank_AD2EAD3_metabolite_objective + worst_core_entry_rank + rank_aux_omics_penalty；其中 worst_core_entry_rank = max(rank_H2_entry_penalty, rank_formate_entry_penalty)，恢复原始 max-rank 规则；objective cutoff 与 ranking 均使用 AD2EAD3_metabolite_objective；total_metabolite_objective 保留为 Stage A 原始可行性审计，不参与 Stage B 筛选。",
    "EAD3-vs-AD2 formate fraction/entry 只作为 annotation，不参与 gate 或 ranking；AD1 不参与本版 Stage A/Stage B 的 objective cutoff、ranking 或 omics penalty。"
  )
))

safe_add_sheet(wb, "01_candidate_audit")
writeData(wb, "01_candidate_audit", candidate_audit)

safe_add_sheet(wb, "02_setting_summary")
writeData(wb, "02_setting_summary", setting_summary)

safe_add_sheet(wb, "03_selected_each_setting")
writeData(wb, "03_selected_each_setting", selected_each_setting)

safe_add_sheet(wb, "04_selected_flux_group")
writeData(wb, "04_selected_flux_group", selected_flux_group)

safe_add_sheet(wb, "05_all_gated_candidates")
writeData(wb, "05_all_gated_candidates", gated_all)

safe_add_sheet(wb, "06_all_gated_flux_group")
writeData(wb, "06_all_gated_flux_group", gated_flux_group)

safe_add_sheet(wb, "07_capacity_mapped")
writeData(wb, "07_capacity_mapped", cap_long)

# One sheet per param_id: retain v results and rankings for all gated candidates under that parameter combination.
for (pid in param_ids) {
  sheet_name <- short_sheet_name(pid)
  tmp <- gated_all %>%
    filter(param_id == pid) %>%
    arrange(objective_cutoff, capacity_mode, penalty_mode, selection_score,
            worst_core_entry_rank, rank_aux_omics_penalty, rank_metabolite_objective,
            total_metabolite_objective)
  safe_add_sheet(wb, sheet_name)
  writeData(wb, sheet_name, tmp)
}

# Prevent Excel save failures caused by occasional Windows/RStudio temporary-directory failures.
.td <- tempdir()
dir.create(.td, recursive = TRUE, showWarnings = FALSE)
# Save to a safe temporary directory first, then copy to the final file to avoid occasional openxlsx temporary-directory errors.
tmp_output_file <- file.path(safe_tmp, basename(output_file))
saveWorkbook(wb, tmp_output_file, overwrite = TRUE)
file.copy(tmp_output_file, output_file, overwrite = TRUE)

cat("\n============================================================\n")
cat("Done Excel:", output_file, "\n")
cat("Done RDS:", output_rds_file, "\n")
cat("重点看：01_candidate_audit, 02_setting_summary, 03_selected_each_setting, 05_all_gated_candidates, 以及每个 P_ 参数 sheet。\n")
cat("============================================================\n")
