library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr, warn.conflicts = FALSE)
library(tidyr)

# Load the arguments (argv)
args = commandArgs(trailingOnly=TRUE)

input_file = args[1]
output_file = sprintf("%s.png", input_file)

df = read_delim(input_file, delim=",", show_col_types = FALSE) %>%
  # Convert string columns to factors
  mutate(bench_name = as.factor(bench_name)) %>%
  mutate(label = as.factor(label)) %>%
  mutate(machine = as.factor(machine)) %>%
  # Make the dataset from:
  #  label,bench_time
  #  baseline,123
  #  current,234
  # Into
  #  time_baseline,time_current
  #  123,234
  pivot_wider(names_from=label, values_from=bench_time, names_prefix="time_") %>%
  # Group by machine and benchmark to compute median
  group_by(machine, bench_name) %>%
    mutate(time_baseline_median = median(time_baseline)) %>%
    mutate(time_norm_baseline = time_baseline / time_baseline_median) %>%
    mutate(time_norm_current  = time_current  / time_baseline_median) %>%
  ungroup() %>%
  # Return time_norm_* columns to a single time column tagged by type
  pivot_longer(
    cols = starts_with("time_norm_"),
    names_to = "type",
    names_prefix = "time_norm_",
    values_to = "time",
    values_drop_na = TRUE
  )

#print(df)

dpi = 200
h = 6
w = 8

# ---------------------------------------------------------------------

p = ggplot(df, aes(x=bench_name, y=time, color=type)) +
  geom_boxplot() +
  scale_color_manual(values = alpha(c("black", "red"))) +
  facet_wrap("machine", ncol=1) +
  theme(axis.text.x = element_text(angle = -20, hjust = 0)) +
  labs(x = "Benchmark", y = "Normalized time", title="Bench6 miniapps")

# TODO: Add commit hash to plot

ggsave(output_file, plot=p, width=w, height=h, dpi=dpi)

print(sprintf("written to %s", output_file))
