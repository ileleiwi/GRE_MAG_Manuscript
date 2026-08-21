## ---------------------------
##
## Script: calculate_redox_ranking_index 
##
## Purpose: construct literature informed ordinal redox ranking for Fig. 5
##          
## Author: Ikaia Leleiwi
##
## Date Created: August 21st, 2026
##
## Copyright (c) Ikaia Leleiwi, 2026
## Email: leleiwi1@llnl.gov
##
## ---------------------------
##
## Notes:
# more negative redox_index values correspond to larger positive reference reduction potentials
# redox_index = -E_ref
##   
##standard electrode potentials
## https://www.nist.gov/system/files/documents/2019/04/02/jpcrd355.pdf
## ---------------------------

## ---------------------------

## Libraries ##
library(tidyverse)

## Functions ##
redox_index <- function(E_ref){
  -E_ref
}

#electron-weighted mean to reconstruct reference potential from component half-reactions
combine_E <- function(E_ref, n_e) {
  stopifnot(length(E_ref) == length(n_e))
  sum(E_ref * n_e) / sum(n_e)
}


reaction_tbl <- tribble(
  ~short_name, ~long_name, ~E_ref, ~n_e, ~source_note,
  
  # Fe3+ + e- -> Fe2+
  "Fe3_Fe2",
  "Iron Reduction Fe(III) => Fe(II)",
  0.771, 1,
  "Bratsch",
  
  "Fe2_Fe3",
  "Iron Oxidation Fe(II) => Fe(III)",
  -0.771, 1,
  "Reverse of Fe3+/Fe2+",
  
  # NO3- + H2O + 2e- -> NO2- + 2OH-
  "NO2_NO3",
  "nitrite => nitrate",
  -0.017, 2,
  "Reverse of nitrate => nitrite",
  
  "NO3_NO2",
  "nitrate => nitrite",
  0.017, 2,
  "Bratsch",
  
  # N2 + 6H2O + 6e- -> 2NH3 + 6OH-
  "N2_NH3",
  "nitrogen => ammonia",
  -0.736, 6,
  "Bratsch",
  
  "NH3_N2",
  "ammonia => nitrogen",
  0.736, 6,
  "Reverse of nitrogen => ammonia",
  
  # NO2- / NH3, OH- couple
  "NH3_NO2",
  "Bacterial/Archaeal ammonia oxidation",
  0.165, 6,
  "Reverse of Bratsch NO2-/NH3,OH- couple; verify exact source value",
  
  "NO2_NH3",
  "nitrite => ammonia",
  -0.165, 6,
  "Bratsch NO2-/NH3,OH- couple; verify exact source value",
  
  # 2SO3^2- + 3H2O + 4e- -> S2O3^2- + 6OH-
  "S2O3_SO3",
  "thiosulfate => sulfite",
  0.566, 4,
  "Reverse of Bratsch SO3^2-/S2O3^2-,OH- couple",
  
  "SO3_S2O3",
  "sulfite => thiosulfate",
  -0.566, 4,
  "Bratsch SO3^2-/S2O3^2-,OH- couple",
  
  # S4O6^2- + 2e- -> 2S2O3^2-
  "S4O6_S2O3",
  "tetrathionate => thiosulfate",
  0.024, 2,
  "Bratsch",
  
  # NO2- + H2O + e- -> NO + 2OH-
  "NO2_NO",
  "nitrite => nitric oxide",
  -0.481, 1,
  "Bratsch",
  
  # 2NO + H2O + 2e- -> N2O + 2OH-
  "NO_N2O",
  "nitric oxide => nitrous oxide",
  0.759, 2,
  "Bratsch",
  
  # NO3- + 10H+ + 8e- -> NH4+ + 3H2O
  "DNRA",
  "Dissimilatory nitrate reduction to ammonium (DNRA)",
  0.880, 8,
  "Bratsch direct NO3-/NH4+ couple",
  
  # N2O + 2H+ + 2e- -> N2 + H2O
  "N2O_N2",
  "nitrous oxide => nitrogen",
  1.769, 2,
  "Bratsch",
  
  # O2 + 4H+ + 4e- -> 2H2O
  "O2_H2O",
  "oxygen => water",
  1.2291, 4,
  "Bratsch oxygen/water couple"
) |>
  mutate(
    redox_index = redox_index(E_ref)
  )

#composite: sulfate -> H2S
#Component Bratsch reference potentials:
#
# SO4^2- -> S0    E_ref = 0.353 V, n = 6 e-
# S0 -> H2S       E_ref = 0.144 V, n = 2 e-

sulfate_steps <- tribble(
  ~step,       ~E_ref, ~n_e,
  "SO4_to_S0",  0.353, 6,
  "S0_to_H2S",  0.144, 2
)

SO4_H2S_E_ref <- combine_E(
  E_ref = sulfate_steps$E_ref,
  n_e = sulfate_steps$n_e
)

SO4_H2S_result <- tibble(
  short_name = "SO4_H2S",
  long_name =
    "dissimilatory sulfate reduction (and oxidation) sulfate => sulfide",
  E_ref = SO4_H2S_E_ref,
  n_e = sum(sulfate_steps$n_e),
  source_note =
    "Composite SO4->S0 and S0->H2S reference potentials; electron-weighted",
  redox_index = redox_index(SO4_H2S_E_ref)
)

#composite thiosulfate -> sulfate
#step1: SO4^2- + 4H+ + 2e- -> SO2 + 2H2O
#step2: 2SO2 + 3H+ + 4e- -> HS2O3- + H2O
#reaction one occurs twice
#2SO4^2- + 8H+ + 4e- -> 2SO2 + 4H2O
#then add rection 2
#2SO4^2- + 11H+ + 8e- -> HS2O3- + 5H2O (sulfate -> thiosulfate)
#reverse process for thiosulfate -> sulfate
#HS2O3- + 5H2O -> 2SO4^2- + 11H+ + 8e-

S2O3_SO4_steps <- tribble(
  ~step,         ~E_ref, ~n_e_per_rxn, ~multiplier,
  "SO4_to_SO2",   0.158, 2,             2,
  "SO2_to_HS2O3", 0.430, 4,             1
) |>
  mutate(
    n_e_combined = n_e_per_rxn * multiplier
  )

S2O3_SO4_E_reduction <- combine_E(
  E_ref = S2O3_SO4_steps$E_ref,
  n_e = S2O3_SO4_steps$n_e_combined
)

S2O3_SO4_E_oxidation <- -S2O3_SO4_E_reduction

S2O3_SO4_result <- tibble(
  short_name = "S2O3_SO4",
  long_name =
    "Thiosulfate oxidation by SOX complex, thiosulfate => sulfate",
  E_ref = S2O3_SO4_E_oxidation,
  n_e = sum(S2O3_SO4_steps$n_e_combined),
  source_note = paste(
    "Composite Bratsch couples:",
    "SO4^2-/SO2 (0.158 V) and SO2/HS2O3- (0.430 V);",
    "electron-weighted and reversed for thiosulfate oxidation"
  ),
  redox_index = redox_index(S2O3_SO4_E_oxidation)
)

#cytochrome ranking heuristic
#all cytochromes use O2 as terminal e- acceptor so use O2/H2O reference potential
cyt_bd <- mean(c(176, 168, 258, 172, 182, 256))

cyt_cbb3 <- mean(
  c(234, 320, 351, -59, 418, 185, 245, 225, 310, 215, -38, 215)
)

cyt_o <- 280

cyt_c <- mean(c(250, 290, 230))

cyt_redox_potentials_mV <- c(
  cyt_bd = cyt_bd,
  cyt_cbb3 = cyt_cbb3,
  cyt_o = cyt_o,
  cyt_c = cyt_c
)


# Min-max scale observed heme midpoint potentials from 0 to 1.
cyt_proportions <- (
  cyt_redox_potentials_mV - min(cyt_redox_potentials_mV)
) / (
  max(cyt_redox_potentials_mV) - min(cyt_redox_potentials_mV)
)


# Anchor all cytochromes to the O2/H2O literature potential.
O2_E_ref <- reaction_tbl |>
  filter(short_name == "O2_H2O") |>
  pull(E_ref)

O2_redox_index <- redox_index(O2_E_ref)


# Separate cytochromes over a 1% range around the oxygen ranking index.
#
cyt_deviation_fraction <- 0.01
cyt_deviation_range <-
  cyt_deviation_fraction * abs(O2_redox_index)

cyt_redox_index_values <-
  O2_redox_index - (cyt_proportions * cyt_deviation_range)

cyt_results <- tibble(
  short_name = names(cyt_redox_index_values),
  
  long_name = c(
    "High affinity: Cytochrome bd ubiquinol oxidase",
    "High affinity: Cytochrome c oxidase, cbb3-type",
    "Low affinity: Cytochrome o ubiquinol oxidase",
    "Low affinity: Cytochrome c oxidase"
  ),
  
  E_ref = O2_E_ref,
  
  n_e = NA_real_,
  
  source_note =
    "Heuristic ranking anchored to O2/H2O E_ref and separated using relative heme midpoint potentials",
  
  redox_index = as.numeric(cyt_redox_index_values)
)


cyt_cp_result <- cyt_results |>
  filter(short_name == "cyt_c") |>
  mutate(
    short_name = "cyt_cp",
    long_name = "Low affinity: Cytochrome c oxidase, prokaryotes"
  )

ranking_df <- bind_rows(
  reaction_tbl,
  SO4_H2S_result,
  S2O3_SO4_result,
  cyt_results,
  cyt_cp_result
) |>
  # O2_H2O is used only as the cytochrome anchor, not as a separate
  # functional category in the manuscript figure.
  filter(short_name != "O2_H2O") |>
  arrange(redox_index) |>
  mutate(
    rank = row_number(),
    
   #dimensionless 0-1 version for plotting/annotation.
    # 0 = smallest redox_index; 1 = largest redox_index.
    redox_index_scaled = (
      redox_index - min(redox_index, na.rm = TRUE)
    ) / (
      max(redox_index, na.rm = TRUE) -
        min(redox_index, na.rm = TRUE)
    )
  ) |>
  select(
    rank,
    short_name,
    long_name,
    E_ref,
    redox_index,
    redox_index_scaled,
    n_e,
    source_note
  )

aerobic_ammonia_oxidation <- ranking_df |>
  filter(short_name == "NH3_NO2") |>
  mutate(
    short_name = "NH3_NO2_aerobic",
    long_name = "Bacterial (aerobic-specific) ammonia oxidation"
  )

ranking_df <- bind_rows(
  ranking_df,
  aerobic_ammonia_oxidation
) |>
  arrange(redox_index) |>
  mutate(rank = row_number())

print(
  ranking_df |>
    select(
      rank,
      short_name,
      long_name,
      E_ref,
      redox_index,
      redox_index_scaled
    ),
  n = Inf
)

dir.create("data", showWarnings = FALSE, recursive = TRUE)

write_tsv(
  ranking_df,
  here("data", "redox_favorability_ranking.tsv")
)
