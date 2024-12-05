BUILDDIR?=.

CC = clang
MPICC = mpicc
WRAPPERS = MPICH_CC=$(CC) I_MPI_CC=$(CC) OMPI_CC=$(CC)

ifdef I_MPI_ROOT
MPICC=mpiicc
endif

CPPFLAGS=-Isrc
CFLAGS=-fompss-2 -Ofast -march=native -ffast-math -std=c11 -D_POSIX_C_SOURCE=199309L

LIBS=-lrt -lm

OPENBLAS_CFLAGS=-I$(OPENBLAS_HOME)/include
OPENBLAS_LDFLAGS=-L$(OPENBLAS_HOME)/lib -L$(OPENBLAS_HOME)/lib64 -lopenblas

TAMPI_CPPFLAGS=-I$(TAMPI_HOME)/include -DTAMPI
TAMPI_LDFLAGS=-ltampi-c -L$(TAMPI_HOME)/lib -Wl,--enable-new-dtags -Wl,-rpath=$(TAMPI_HOME)/lib

BIN= \
	$(BUILDDIR)/01.matmul_ompss2_tampi.bin \
	$(BUILDDIR)/02.matmul_ompss2_itampi.bin

all: $(BUILDDIR) $(BIN)

$(BUILDDIR)/01.matmul_ompss2_tampi.bin: src/common/matmul.c src/mpi/main.c src/mpi/01.matmul_ompss2_tampi.c
	$(WRAPPERS) $(MPICC) $(CPPFLAGS) $(CFLAGS) $(OPENBLAS_CFLAGS) $(TAMPI_CPPFLAGS) $^ -o $@ $(TAMPI_LDFLAGS) $(OPENBLAS_LDFLAGS) $(LIBS)

$(BUILDDIR)/02.matmul_ompss2_itampi.bin: src/common/matmul.c src/mpi/main.c src/mpi/02.matmul_ompss2_itampi.c
	$(WRAPPERS) $(MPICC) $(CPPFLAGS) $(CFLAGS) $(OPENBLAS_CFLAGS) $(TAMPI_CPPFLAGS) $^ -o $@ $(TAMPI_LDFLAGS) $(OPENBLAS_LDFLAGS) $(LIBS)

$(BUILDDIR):
	mkdir $(BUILDDIR)

clean:
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.bin
