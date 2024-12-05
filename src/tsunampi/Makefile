CC=clang
MPICC=I_MPI_CC=$(CC) OMPI_CC=$(CC) mpicc
MPILINK=$(MPICC)

CFLAGS+=-O3

CFLAGS+=-fompss-2
LDFLAGS+=-fompss-2
LIBS+=-lm 

BIN=tsunampi.tampi \
    tsunampi.tagaspi

tampi_LIBS+=-lstdc++ -l:libtampi-c.a
tagaspi_LIBS+=-lGPI2 -ltagaspi

all: $(BIN) copy

tsunampi.%: main.c %.c
	$(MPICC) $(CFLAGS) $($*_CFLAGS) $(LDFLAGS) $($*_LDFLAGS) $^ -o $@ $(LIBS) $($*_LIBS)

copy: $(BIN)
	scp $^ mn3:tsunampi/
