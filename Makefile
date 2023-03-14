CC=clang
CFLAGS=-O3 -fompss-2

BENCHMARKS=\
	sched_get \
	sched_add \
	register_deps \
	readywave

BIN=$(addprefix bench6.,$(BENCHMARKS)) 
DATA=$(addsuffix .csv, $(addprefix data/,$(BENCHMARKS)))
PLOT=$(DATA:=.png)

all: $(BIN) $(DATA) $(PLOT)

bench6.%: src/%.c src/common.c
	$(CC) $(CFLAGS) -o $@ $^

data/%.csv.png: data/%.csv plot/%.R
	Rscript plot/$(*F).R $<

data/%.csv: bench6.%
	mkdir -p data
	./$^ > $@
