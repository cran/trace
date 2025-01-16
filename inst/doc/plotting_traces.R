## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=7.2, fig.height=6) 

## ----setup,warning=FALSE, message=FALSE---------------------------------------
library(trace)
library(dplyr)
library(ggplot2)

## ----fsa_import---------------------------------------------------------------
fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())


## ----ladders,warning=FALSE, message=FALSE-------------------------------------
find_ladders(
  fsa_list,
  show_progress_bar = FALSE
)

## ----find_fragments, warning=FALSE--------------------------------------------
fragments_list <- find_fragments(
  fsa_list,
  min_bp_size = 300
)

## ----add_metadata-------------------------------------------------------------

add_metadata(
  fragments_list = fragments_list,
  metadata_data.frame = metadata,
  unique_id = "unique_id",
  metrics_group_id = "metrics_group_id",
  metrics_baseline_control = "metrics_baseline_control",
  batch_run_id = "batch_run_id",
  batch_sample_id = "batch_sample_id"
)

## ----find_alleles_and_repeats, warning=FALSE, message=FALSE-------------------
find_alleles(
  fragments_list = fragments_list
)


call_repeats(
  fragments_list = fragments_list
)

## -----------------------------------------------------------------------------
plot_traces(fragments_list[1], xlim = c(100, 150), show_peaks = FALSE)

## -----------------------------------------------------------------------------
extract_trace_table(fragments_list[1]) |>
  filter(between(calculated_repeats, 100, 150)) |>
  ggplot(aes(calculated_repeats, signal)) +
  geom_line() 

## -----------------------------------------------------------------------------

traces_df <- extract_trace_table(fragments_list)


## -----------------------------------------------------------------------------

samples_for_plotting <- metadata[which(metadata$metrics_group_id == "CC6"), "unique_id"]

traces_for_plotting_df <- extract_trace_table(fragments_list[samples_for_plotting])


## -----------------------------------------------------------------------------
traces_for_plotting_with_metadata_df <- dplyr::left_join(
  traces_for_plotting_df,
  metadata,
  by = join_by(unique_id)
)

## -----------------------------------------------------------------------------
traces_for_plotting_with_metadata_df |>
  dplyr::filter(between(calculated_repeats, 100, 150)) |>
  ggplot(aes(x = calculated_repeats,
             y = signal,
             colour = as.factor(treatment))) +
  geom_line(aes(group = unique_id)) +
  facet_wrap(vars(paste("Day", day, treatment, "nM Branaplam")), ncol = 1)




## -----------------------------------------------------------------------------

alleles_df <- extract_alleles(fragments_list[samples_for_plotting])

sample_to_exclude <- alleles_df[which(is.na(alleles_df$allele_repeat)), "unique_id"]

traces_for_plotting_with_metadata_df <- traces_for_plotting_with_metadata_df |>
  dplyr::filter(between(calculated_repeats, 100, 150),
                !unique_id %in% sample_to_exclude) |>
  group_by(unique_id) |>
  mutate(relative_signal = signal / max(signal),
         day_treatment = paste("Day", day, treatment, "nM Branaplam")) |>
  ungroup()


traces_for_plotting_with_metadata_df |>
  ggplot(aes(x = calculated_repeats,
             relative_signal,
             colour = as.factor(treatment))) +
  geom_line(aes(group = unique_id)) +
  facet_wrap(vars(day_treatment), ncol = 1)



## -----------------------------------------------------------------------------


d0_trace <- traces_for_plotting_with_metadata_df |>
  filter(day == 0) |>
  filter(unique_id == unique(unique_id)[1]) |>
  select(-day_treatment)

traces_for_plotting_with_metadata_df |>
  group_by(day, treatment) |>
  filter(unique_id == unique(unique_id)[1]) |>
  ggplot(aes(x = calculated_repeats,
             y = relative_signal,
             colour = paste("Day", day, treatment, "nM Branaplam"))) +
  geom_line(data = d0_trace ,
            aes(group = unique_id),
            colour = "gray40") +
  geom_line(aes(group = unique_id)) +
  facet_wrap(vars(day_treatment), ncol = 1) +
  labs(colour = "") +
  theme_minimal()
  



