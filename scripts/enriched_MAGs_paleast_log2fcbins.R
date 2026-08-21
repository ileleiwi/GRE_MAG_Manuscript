## ---------------------------
##
## Script: enriched_MAGs_paleast_log2fcbins
##
## Purpose: find set of MAG that are enriched in each treatment
##          
## Author: Ikaia Leleiwi
##
## Date Created: February 26th, 2024
##
## Copyright (c) Ikaia Leleiwi, 2024
## Email: leleiwi1@llnl.gov
##
## ---------------------------
##
## Notes: Method used is the present absent least conservative method (pa_least) and Log2fc
## Least conservative (pa_least)
## Genome has to be found in a single 13C bin and not in the equivalent 12C bin and found in both reps, 
## 12C can only be found in bins lighter than the 13C bin
## 
## Log2fc
## Multiple instances from a single sample in a fraction bin were averaged before being input into DeSeq2
## Genome is present if adjpval < 0.05 and log2fc > 0 comparing 13C vs 12C
## AN has no instances of 12C in density bin 4, so that comparison could not be made
## 
##   
##
## ---------------------------

## Libraries ##
library(tidyverse)
library(janitor)
library(ggridges)
library(readxl)
library(here)
## Data ##

meta <- read_xlsx(here("data", "metadata", "SIP_BULK_MG_sampleData_short.xlsx")) |>
  clean_names() |>
  filter(density != "Bulk") |>
  mutate(fraction = ifelse(fraction == "NA", NA, fraction),
         density = as.double(density)) 

rename_cols <- function(x){
  
  # x <- str_remove(x, "GRE\\.bulkMG\\.")
  # x <- str_remove(x, "GRE\\.SIPMG\\.")
  x <- str_remove(x, "_Mean")
  
  return(x)
}

rename_cols_getmm <- function(x){
  
  x <- str_remove(x, "GRE\\.bulkMG\\.")
  x <- str_remove(x, "GRE\\.SIPMG\\.")
  
  return(x)
}


rabund <- read_tsv(here("data", "drep_mags_relabund_10cov_95id.tsv")) |>
  rename_all(rename_cols) |>
  filter(if_any(!matches("bins"), ~ . != 0)) |>
  dplyr::rename("genome" = "bins") |>
  pivot_longer(cols = -genome,
               names_to = "jgi_sample_id",
               values_to = "relative_abundance") 

getmm <- read_tsv(here("data", "drep_mags_GeTMM_10cov_95id.tsv")) |>
  filter(bins %in% rabund$genome) |>
  rename_all(rename_cols) |>
  rename_all(rename_cols_getmm)

gc <- read_tsv(here("data", "all_mags_gc.tsv"))

deseq <- read_tsv(here("data", "DESeq2_results_MAGs.tsv"))


meta |>
  ggplot(aes(x=density)) +
  geom_histogram(aes(color = isotope),
                 bins = 300) +
  facet_wrap(~treatment)

split_df <- function(df, n) {
  ## function that breaks a dataframe into multiple dataframes where nrow = n
  
  if(n > nrow(df)) { break }
  
  rows_per_df <- ceiling(nrow(df)/n)
  
  df_list <- split(df, rep(1:n, each = rows_per_df, length.out = nrow(df)))
  
  return(df_list)
  
}

count_splits <- function(data) {
  ## takes a dataframe, subsets it into a list of dataframes on nrow = n
  ## checks that each subset density range contains at least one instance of 'present' from 12C and 13C
  ## returns vector of all the nrow values where there is a 12C and 13C 'present' instance in all of the density bins
  n_valid_split_vec <- c() 
  for(i in 1:nrow(data)){
    df_list <- split_df(data, i) #makes a list of dataframes 
    df_bool_vec_12 <- c()
    df_bool_vec_13 <- c()
    for(l in 1:length(df_list)){ #checks that each dataframe in the list has 12C and 13C observations
      df_bool_vec_12 <- c(df_bool_vec_12, sum(df_list[[l]]$`12C`) > 0) 
      df_bool_vec_13 <- c(df_bool_vec_13, sum(df_list[[l]]$`13C`) > 0) 
    }
    if(sum(df_bool_vec_12) == length(df_bool_vec_12) &
       sum(df_bool_vec_13) == length(df_bool_vec_13)) {
      n_valid_split_vec <- c(n_valid_split_vec, i) #save all valid splits in vector, valid split = # of df's in df_list = potential number of bins for density
    }
  }
  return(n_valid_split_vec)
}


meta_by_trt <- meta |> 
  group_by(treatment) |>
  nest() |>
  mutate(data = map(data, ~ dplyr::count(.x, density, isotope) |>
                      arrange(density) |>
                      pivot_wider(names_from = "isotope",
                                  values_from = "n",
                                  values_fill = 0) |>
                      count_splits())) |>
  unnest(cols = c(data)) |>
  rownames_to_column(var = "idx") |>
  pivot_wider(names_from = "treatment",
              values_from = "data")

#AN treatment has lowest number of valid bins n=4

an_bins <- meta |> 
  filter(treatment == "AN",
         !is.na(density)) |>
  arrange(density) |>
  split_df(4)

lf_bins <- meta |> 
  filter(treatment == "LF",
         !is.na(density)) |>
  arrange(density) |>
  split_df(6)

for(i in an_bins){ print(range(i$density))}

for(i in lf_bins){ print(range(i$density))}

#these are the bins!
dens_bins_df_4 <- meta |> 
  arrange(density) |>
  mutate(dens_bin = case_when(density <= 1.7340790790 ~ 1,
                              density > 1.7340790790 & density < 1.736811 ~ 2,
                              density >= 1.736811 & density <= 1.7396794740 ~ 3,
                              density > 1.7396794740 ~ 4,
                              T ~ 0),
         dens_bin = factor(dens_bin, levels = seq(1,4)),
         rep = factor(rep, levels = seq(1,3)))


dens_bins_df_4 |>
  ggplot(aes(x = density)) +
  geom_histogram(aes(color = rep),
                 bins = 100) +
  facet_wrap(~dens_bin+treatment)

dens_bins_df_4 |>
  ggplot(aes(x= density, y = treatment, fill = rep, color = rep)) +
  geom_density_ridges(stat = "binline", bins = 300) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.0001), breaks = seq(1.731511, 1.741564, by = 0.0005)) +
  theme_ridges() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  facet_wrap(~dens_bin + isotope)



#present absent method
convert_to_binary <- function(df) {
  df %>%
    mutate_if(is.numeric, ~ if_else(. > 0, 1, 0))
}


pa_df_trt <- rabund |>
  left_join(dens_bins_df_4, by = "jgi_sample_id") |>
  filter(!is.na(dens_bin)) |>
  mutate(present_absence = ifelse(relative_abundance > 0, 1, 0),
         rep = paste0("rep_", as.character(rep)),
         rep = factor(rep, levels = paste0("rep_",seq(1,3)))) |>
  group_by(genome, isotope, dens_bin, rep, treatment) |>
  summarise(total = sum(present_absence), .groups = "drop") |>
  pivot_wider(names_from = "dens_bin",
              values_from = "total",
              values_fill = 0) |>
  group_by(isotope, rep, treatment) |>
  nest() |>
  mutate(data = map(data, ~ convert_to_binary(.x) |>
                      mutate(total = rowSums(across(where(is.numeric)))))) |>
  unnest(cols = c(data)) 



remove_commas <- function(x) {
  # Remove leading and trailing commas and spaces
  x <- gsub("^,\\s*|,\\s*$", "", x)
  
  # Remove adjacent commas separated by space
  x <- gsub("\\s*,\\s*", ",", x)
  
  # Remove adjacent commas with empty entries
  x <- gsub(",(?=,)", "", x, perl = TRUE)
  
  # Remove any trailing commas
  x <- gsub(",$", "", x)
  
  # Remove any leading commas
  x <- gsub("^,", "", x)
  
  return(x)
}

pa_isotope_genome_rep_trt <- pa_df_trt |>
  pivot_longer(cols = c(`1`, `2`, `3`, `4`),
               names_to = "bin",
               values_to = "present_in_bin") |>
  mutate(bin = as.character(bin)) |>
  filter(present_in_bin == 1) |>
  mutate(bin = ifelse(isotope == "12C", paste0(bin, "_l"), paste0(bin, "_h"))) |> #light, heavy
  group_by(genome, isotope, rep, treatment) |>
  summarise(all_bins_present_in = paste(bin, collapse = ", "),
            .groups = "drop") |>
  pivot_wider(names_from = "rep",
              values_from = "all_bins_present_in",
              values_fill = "") |>
  pivot_longer(cols = starts_with("rep"),
               names_to = "rep",
               values_to = "all_bins_present_in") |>
  pivot_wider(names_from = "treatment",
              values_from = "all_bins_present_in",
              values_fill = "") |>
  pivot_longer(cols = c("AN", "LF", "HF", "OX"),
               names_to = "treatment",
               values_to = "all_bins_present_in") |>
  rowwise() |>
  mutate(total_bins_in = nchar(all_bins_present_in),
         total_bins_in = case_when(total_bins_in == 0 ~ 0,
                                   total_bins_in == 3 ~ 1,
                                   total_bins_in == 8 ~ 2,
                                   total_bins_in == 13 ~ 3,
                                   total_bins_in == 18 ~ 4)) |>
  pivot_wider(names_from = "isotope",
              values_from = "total_bins_in",
              values_fill = 0) |>
  pivot_wider(names_from = "rep",
              values_from = "all_bins_present_in",
              values_fill = "") |>
  group_by(genome, treatment) |>
  summarise(`total_bins_13C` = sum(`13C`),
            `total_bins_12C` = sum(`12C`),
            rep_1 = paste(rep_1, collapse = ", "),
            rep_2 = paste(rep_2, collapse = ", "),
            rep_3 = paste(rep_3, collapse = ", "),
            .groups = "drop") |>
  filter(total_bins_13C > 0) |>
  mutate_at(vars(starts_with("rep_")), ~remove_commas(.)) |>
  mutate_at(vars(starts_with("rep_")), ~str_replace_all(., ",", ", ")) |>
  mutate(total_reps_13C = case_when(str_detect(rep_1, "h") & str_detect(rep_2, "h") & str_detect(rep_3, "h") ~ 3,
                                    str_detect(rep_1, "h") & str_detect(rep_2, "h") & !str_detect(rep_3, "h") ~ 2,
                                    str_detect(rep_1, "h") & !str_detect(rep_2, "h") & str_detect(rep_3, "h") ~ 2,
                                    str_detect(rep_1, "h") & !str_detect(rep_2, "h") & !str_detect(rep_3, "h") ~ 1,
                                    !str_detect(rep_1, "h") & str_detect(rep_2, "h") & !str_detect(rep_3, "h") ~ 1,
                                    !str_detect(rep_1, "h") & str_detect(rep_2, "h") & str_detect(rep_3, "h") ~ 2,
                                    !str_detect(rep_1, "h") & !str_detect(rep_2, "h") & str_detect(rep_3, "h") ~ 1,
                                    !str_detect(rep_1, "h") & !str_detect(rep_2, "h") & str_detect(rep_3, "h") ~ 0),
         total_reps_12C = case_when(str_detect(rep_1, "l") & str_detect(rep_2, "l") & str_detect(rep_3, "l") ~ 3,
                                    str_detect(rep_1, "l") & str_detect(rep_2, "l") & !str_detect(rep_3, "l") ~ 2,
                                    str_detect(rep_1, "l") & !str_detect(rep_2, "l") & !str_detect(rep_3, "l") ~ 1,
                                    !str_detect(rep_1, "l") & str_detect(rep_2, "l") & !str_detect(rep_3, "l") ~ 1,
                                    !str_detect(rep_1, "l") & str_detect(rep_2, "l") & str_detect(rep_3, "l") ~ 2,
                                    !str_detect(rep_1, "l") & !str_detect(rep_2, "l") & str_detect(rep_3, "l") ~ 1,
                                    !str_detect(rep_1, "l") & !str_detect(rep_2, "l") & !str_detect(rep_3, "l") ~ 0)) |>
  left_join(gc, by = c("genome" = "fasta")) |>
  dplyr::select(genome, starts_with("total"), everything())

count_bins <- function(x){
  
  h <- c(str_count(x, "1_h"),
         str_count(x, "2_h"),
         str_count(x, "3_h"),
         str_count(x, "4_h"))
  
  l <- c(str_count(x, "1_l"),
         str_count(x, "2_l"),
         str_count(x, "3_l"),
         str_count(x, "4_l"))
  
  out <- list(h, l)
  names(out) <- c("h", "l")
  return(out)
}

build_lists <- function(r1, r2, r3){
  list1 <- count_bins(r1)
  list2 <- count_bins(r2)
  list3 <- count_bins(r3)
  
  return(list(list1 = list1, list2 = list2, list3 = list3))
}

condition_pass_fail <- function(p1, p2, p3){
  
  condition_pass <- case_when(p1 & p2 ~ 1,
                              p1 & p3 ~ 1,
                              p2 & p3 ~ 1,
                              T ~ 0)
  
  return(condition_pass)
  
}


least_conservative <- function(r1, r2, r3){
  # Genome has to be found in a single 13C bin and not in the equivalent 12C bin and found in both reps
  # 12C can only be found in bins lighter than the 13C bin
  list_all <- build_lists(r1, r2, r3)
  
  
  l1_h <- which(list_all$list1[["h"]] == 1)
  l1_l <- which(list_all$list1[["l"]] == 1)
  l2_h <- which(list_all$list2[["h"]] == 1)
  l2_l <- which(list_all$list2[["l"]] == 1)
  l3_h <- which(list_all$list3[["h"]] == 1)
  l3_l <- which(list_all$list3[["l"]] == 1)
  
  if(length(l1_l) > 0 & length(l1_h) > 0) {
    l1_pass <- ifelse(length(l1_h) > 0 & max(l1_h) > max(l1_l), T, F)
  }else if(length(l1_l) > 0 & length(l1_h) == 0){
    l1_pass <- F
  }else{
    l1_pass <- ifelse(length(l1_h) > 0, T, F)
  }
  
  if(length(l2_l) > 0 & length(l2_h) > 0) {
    l2_pass <- ifelse(length(l2_h) > 0 & max(l2_h) > max(l2_l), T, F)
  }else if(length(l2_l) > 0 & length(l2_h) == 0){
    l2_pass <- F
  }else{
    l2_pass <- ifelse(length(l2_h) > 0, T, F)
  }
  
  if(length(l3_l) > 0 & length(l3_h) > 0) {
    l3_pass <- ifelse(length(l3_h) > 0 & max(l3_h) > max(l3_l), T, F)
  }else if(length(l3_l) > 0 & length(l3_h) == 0){
    l3_pass <- F
  }else{
    l3_pass <- ifelse(length(l3_h) > 0, T, F)
  }
  
  return(condition_pass_fail(l1_pass, l2_pass, l3_pass))
  
}


pa_final_df_trt <- pa_isotope_genome_rep_trt |>
  group_by(treatment) |>
  nest() |>
  mutate(data = map(data, ~rowwise(.x) |>
                      mutate(least = least_conservative(rep_1, rep_2, rep_3)))) |>
  unnest(cols = "data")

pa_final_df_trt |>
  filter(least == 1,
         treatment == "AN") 


pa_final_df_trt |>
  filter(least == 1) |>
  group_by(treatment) |>
  dplyr::count(least) |>
  dplyr::select(-least) 

pa_final_df_trt2 <- pa_final_df_trt |>
  mutate(treatment = paste0("pa_", treatment)) |>
  group_by(genome, treatment) |>
  summarise(total = sum(least), .groups = "drop") |>
  pivot_wider(names_from = "treatment",
              values_from = "total",
              values_fill = 0)

#log2fc
rename_cols_deseq <- function(x){
  trt <- str_extract(x, "^[A-Z]{2}")
  bin <- str_extract(x, "\\d")
  paste0("log2fc_", bin, "_", trt)
}


deseq_df <- deseq |>
  mutate(present_absent = ifelse(padj < 0.05 &
                                   log2FoldChange > 0, 1, 0)) |>
  group_by(row, id) |>
  summarise(total = sum(present_absent), .groups = "drop") |>
  pivot_wider(names_from = "id",
              values_from = "total",
              values_fill = 0) |>
  dplyr::select(genome = row,
                contains("_")) %>%
  rename_with(.fn = ~rename_cols_deseq(.), .cols = c(starts_with("AN"), 
                                                     starts_with("LF"),
                                                     starts_with("HF"), 
                                                     starts_with("OX"))) |>
  pivot_longer(cols = -genome,
               names_to = "method_bin_trt",
               values_to = "pa") |>
  separate(col = "method_bin_trt",
           into = c("method", "bin", "treatment"),
           sep = "_",
           remove = TRUE) |>
  group_by(genome, treatment) |>
  summarise(total = sum(pa), .groups = "drop") |>
  mutate(treatment = paste0("log2fc_", treatment)) |>
  pivot_wider(names_from = "treatment",
              values_from = "total",
              values_fill = 0)

pa_log2fc_df <- pa_final_df_trt2 |>
  left_join(deseq_df) |>
  pivot_longer(cols = -genome,
               names_to = "method_trt",
               values_to = "pa") |>
  separate(col = "method_trt",
           into = c("method", "treatment"),
           sep = "_",
           remove = TRUE) |>
  group_by(treatment) |>
  nest() |>
  mutate(data = map(data, ~group_by(.x, genome) |>
                      summarise(total = sum(pa)) |>
                      mutate(total = ifelse(total > 0, 1, 0)))) |>
  unnest(cols = "data") |>
  pivot_wider(names_from = "treatment",
              values_from = "total")

write_tsv(
  pa_log2fc_df,
  here("data", "enriched_mags.tsv")
)





pa_final_df_trt2 |>
  left_join(deseq_df) |>
  pivot_longer(cols = starts_with("log"),
               names_to = "ds2",
               values_to = "num") |>
  mutate(num = ifelse(num > 0, 1, num)) |>
  pivot_wider(names_from = "ds2",
              values_from = "num",
              values_fill = 0) |>
  pivot_longer(cols = -genome,
               names_to = "name",
               values_to = "yes_no") |>
  separate(name, into = c("method", "treatment"), sep = "_") |>
  filter(yes_no == 1) |>
  group_by(treatment) |>
  nest() |>
  mutate(
    method_sets = map(data, ~ .x %>% group_by(method) %>% summarise(genomes = list(unique(genome)), .groups = "drop")),
    overlap_stats = map(method_sets, function(df) {
      methods <- df$method
      if (length(methods) < 2) return(tibble(
        shared = 0,
        only_1 = length(df$genomes[[1]]),
        only_2 = 0,
        method1 = methods[1],
        method2 = NA
      ))
      
      genomes1 <- df$genomes[[1]]
      genomes2 <- df$genomes[[2]]
      
      shared <- intersect(genomes1, genomes2)
      only_1 <- setdiff(genomes1, genomes2)
      only_2 <- setdiff(genomes2, genomes1)
      
      tibble(
        shared = length(shared),
        only_1 = length(only_1),
        only_2 = length(only_2),
        method1 = methods[1],
        method2 = methods[2]
      )
    })
    
  ) |>
  select(treatment, overlap_stats) |>
  unnest(overlap_stats)

pa_final_df_trt2 |>
  left_join(deseq_df) |>
  pivot_longer(cols = starts_with("log"),
               names_to = "ds2",
               values_to = "num") |>
  mutate(num = ifelse(num > 0, 1, num)) |>
  pivot_wider(names_from = "ds2",
              values_from = "num",
              values_fill = 0) |>
  pivot_longer(cols = -genome,
               names_to = "name",
               values_to = "yes_no") |>
  separate(name, into = c("method", "treatment"), sep = "_") |>
  filter(yes_no == 1) |>
  group_by(treatment, method) |>
  summarise(total = sum( yes_no)) |>
  pivot_wider(names_from = "treatment",
              values_from = "total")
