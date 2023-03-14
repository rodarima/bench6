library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)

# Load the arguments (argv)
args = commandArgs(trailingOnly=TRUE)

input_file = args[1]

df = read_delim(input_file, delim=",", show_col_types = FALSE)

dpi = 96
h = 2
w = 6

# ---------------------------------------------------------------------

p = ggplot(df, aes(time_ms)) +
  geom_histogram(bins=80) +
  theme_bw() +
  labs(x = "Time (ms)", title="Nanos6: readywave time")
  # TODO: Add ntasks and taskwork to labels

ggsave(sprintf("%s.png", input_file), plot=p, width=w, height=h, dpi=dpi)
