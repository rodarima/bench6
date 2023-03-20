CC=clang
CFLAGS=-O3 -fompss-2 -pthread
LDFLAGS=-lm
#CFLAGS+=-g -fno-omit-frame-pointer -gdwarf-4
DESTDIR?=/usr

BENCHMARKS=\
	sched_get \
	sched_add \
	register_deps \
	readywave

BIN=$(addprefix bench6.,$(BENCHMARKS)) 
DATA=$(addsuffix .csv, $(addprefix data/,$(BENCHMARKS)))
PLOT=$(DATA:=.png)

all: $(BIN)

clean:
	rm -f $(BIN)

install: $(BIN)
	mkdir -p $(DESTDIR)/bin
	cp bench6* $(DESTDIR)/bin

plot: $(BIN) $(DATA) $(PLOT)

bench6.%: src/%.c src/common.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

data/%.csv.png: data/%.csv plot/%.R
	Rscript plot/$(*F).R $<

data/%.csv: bench6.%
	mkdir -p data
	./$^ > $@
