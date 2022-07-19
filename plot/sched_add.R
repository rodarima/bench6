library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)

# Load the arguments (argv)
args = commandArgs(trailingOnly=TRUE)

input_file = args[1]

df = read_delim(input_file, delim=",", comment="#") %>%
  mutate(run = as.factor(run)) %>%
  mutate(time_per_task_per_cpu = time_per_task_per_cpu * 1e9)

dpi = 150
h = 2
w = 6

# ---------------------------------------------------------------------

#p = ggplot(df, aes(x=run, y=time_per_task)) +
#  geom_point() +
#  theme_bw() +
#  labs(
#       x = "Number of run",
#       y="get_ready_task() time (ns)",
#       title="Nanos6: get ready task time")

p = ggplot(df, aes(x=time_per_task_per_cpu)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_continuous(breaks = breaks_pretty(10)) +
  labs(
       x="Duration per task per CPU (ns / task * CPU)",
       title="bench6.sched_add: time to unblock N tasks")

ggsave(sprintf("%s.png", input_file), plot=p, width=w, height=h, dpi=dpi)

