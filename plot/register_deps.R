library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)

# Load the arguments (argv)
args = commandArgs(trailingOnly=TRUE)

input_file = args[1]

df = read_delim(input_file, delim=",", comment="#") %>%
  mutate(ndeps = as.factor(ndeps)) %>%
  mutate(time_per_task = time_per_task * 1e6) %>%
  group_by(ndeps) %>%
  mutate(median_time = median(time_per_task)) %>%
  ungroup()

dpi = 300
h = 6
w = 15

# ---------------------------------------------------------------------

p = ggplot(df, aes(x=ndeps, y=time_per_task)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x="Number of dependencies per task",
    y="Creation and registration time per task (us)",
    title="bench6.register_deps: registration time vs number of dependencies")

ggsave(sprintf("%s.png", input_file), plot=p, width=w, height=h, dpi=dpi)

