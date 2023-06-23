library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)

# Load the arguments (argv)
args = commandArgs(trailingOnly=TRUE)

input_file = "data/readywave-instr.csv"

df = read_delim(input_file, delim=",", show_col_types = FALSE) %>%
  mutate(instr = as.factor(instr))

dpi = 150
h = 2
w = 7

# ---------------------------------------------------------------------

p = ggplot(df, aes(time_ms, fill=instr)) +
  geom_histogram(color="white", bins=50) +
  #theme_bw() +
  labs(x = "Time (ms)", title="bench6.readywave -r 100 -t 5000 -w 10")
  # TODO: Add ntasks and taskwork to labels

ggsave(sprintf("%s.png", input_file), plot=p, width=w, height=h, dpi=dpi)
