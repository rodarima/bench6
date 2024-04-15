# Compilers
CC=clang
MCC=clang
MPICC=mpicc
CC_WRAPPERS=I_MPI_CC=$(CC) MPICH_CC=$(CC) OMPI_CC=$(CC)
MCC_WRAPPERS=I_MPI_CC=$(MCC) MPICH_CC=$(MCC) OMPI_CC=$(MCC)

# Disable SIMD kernel by default
SIMD?=0

ifdef I_MPI_ROOT
MPICC=mpiicc
endif

# Preprocessor flags
CPPFLAGS=-Isrc

# Compiler flags
CFLAGS=-Ofast -march=native -ffast-math -std=gnu99
MCCFLAGS=-fompss-2 -Ofast -march=native -ffast-math 

# Kernel flags
KFLAGS=$(CFLAGS)

# Compile kernel with SIMD if requested
ifneq ($(SIMD),0)
CPPFLAGS+=-DSIMD
KFLAGS+=-fopenmp-simd
endif

# Linker flags
LDFLAGS=-lrt -lm

# TAMPI flags
TAMPI_CPPFLAGS=-I$(TAMPI_HOME)/include -DTAMPI
TAMPI_LDFLAGS=-L$(TAMPI_HOME)/lib -L$(TAMPI_HOME)/lib64 -l:libtampi-c.a -lpthread -lstdc++

# GASPI flags
GASPI_CPPFLAGS=-I$(GASPI_HOME)/include
GASPI_LDFLAGS=-L$(GASPI_HOME)/lib -L$(GASPI_HOME)/lib64 -l:libGPI2.a -lpthread -lrt -libverbs

# GASPI extension flags
GASPI_EXT_CPPFLAGS=-I$(GASPI_EXT_HOME)/include
GASPI_EXT_LDFLAGS=-L$(GASPI_EXT_HOME)/lib -L$(GASPI_EXT_HOME)/lib64 -l:libGPI2.a -lpthread -lrt -libverbs

# TAGASPI flags
TAGASPI_CPPFLAGS=-I$(TAGASPI_HOME)/include -DTAGASPI
TAGASPI_LDFLAGS=-L$(TAGASPI_HOME)/lib -L$(TAGASPI_HOME)/lib64 -l:libtagaspi.a -lnuma -lstdc++

# List of programs
BIN=01.heat_seq.bin \
	02.heat_ompss2.bin \
	03.heat_ompss2_residual.bin \
	04.heat_mpi.bin \
	05.heat_mpi_nbuffer.bin \
	06.heat_mpi_ompss2_forkjoin.bin \
	07.heat_mpi_ompss2_tasks.bin \
	15.heat_mpirma_ompss2_tasks.bin \
	18.heat_mpirma_nbuffer.bin

ifdef TAMPI_HOME
BIN+=08.heat_tampi_ompss2_tasks.bin
BIN+=09.heat_itampi_ompss2_tasks.bin
ifdef IFENCE
BIN+=16.heat_tampirma_ompss2_tasks.bin
BIN+=17.heat_tampirma_ompss2_tasks.bin
endif
endif

ifdef GASPI_HOME
BIN+=10.heat_gaspi.bin
BIN+=11.heat_gaspi_nbuffer.bin
BIN+=12.heat_gaspi_ompss2_forkjoin.bin
BIN+=13.heat_gaspi_ompss2_tasks.bin
endif

ifdef GASPI_EXT_HOME
ifdef TAGASPI_HOME
BIN+=14.heat_tagaspi_ompss2_tasks.bin
endif
endif

# Sources
SMP_SRC=src/common/misc.c src/smp/main.c kernel.o
MPI_SRC=src/common/misc.c src/mpi/utils.c src/mpi/main.c kernel.o
GASPI_SRC=src/common/misc.c src/gaspi/utils.c src/gaspi/main.c kernel.o

all: $(BIN)

01.heat_seq.bin: $(SMP_SRC) src/smp/01.solver_seq.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

02.heat_ompss2.bin: $(SMP_SRC) src/smp/02.solver_ompss2.c
	$(MCC) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

03.heat_ompss2_residual.bin: $(SMP_SRC) src/smp/03.solver_ompss2_residual.c
	$(MCC) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

04.heat_mpi.bin: $(MPI_SRC) src/mpi/04.solver_mpi.c
	$(CC_WRAPPERS) $(MPICC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

05.heat_mpi_nbuffer.bin: $(MPI_SRC) src/mpi/05.solver_mpi_nbuffer.c
	$(CC_WRAPPERS) $(MPICC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

06.heat_mpi_ompss2_forkjoin.bin: $(MPI_SRC) src/mpi/06.solver_mpi_ompss2_forkjoin.c
	$(MCC_WRAPPERS) $(MPICC) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

07.heat_mpi_ompss2_tasks.bin: $(MPI_SRC) src/mpi/07.solver_mpi_ompss2_tasks.c
	$(MCC_WRAPPERS) $(MPICC) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

08.heat_tampi_ompss2_tasks.bin: $(MPI_SRC) src/mpi/08.solver_tampi_ompss2_tasks.c
	$(MCC_WRAPPERS) $(MPICC) $(CPPFLAGS) $(TAMPI_CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(TAMPI_LDFLAGS)

09.heat_itampi_ompss2_tasks.bin: $(MPI_SRC) src/mpi/09.solver_itampi_ompss2_tasks.c
	$(MCC_WRAPPERS) $(MPICC) $(CPPFLAGS) $(TAMPI_CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(TAMPI_LDFLAGS)

10.heat_gaspi.bin: $(GASPI_SRC) src/gaspi/10.solver_gaspi.c
	$(CC_WRAPPERS) $(MPICC) $(CPPFLAGS) $(GASPI_CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(GASPI_LDFLAGS)

11.heat_gaspi_nbuffer.bin: $(GASPI_SRC) src/gaspi/11.solver_gaspi_nbuffer.c
	$(CC_WRAPPERS) $(MPICC) -DNBUFFER $(CPPFLAGS) $(GASPI_CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(GASPI_LDFLAGS)

12.heat_gaspi_ompss2_forkjoin.bin: $(GASPI_SRC) src/gaspi/12.solver_gaspi_ompss2_forkjoin.c
	$(MCC_WRAPPERS) $(MPICC) $(CPPFLAGS) $(GASPI_CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(GASPI_LDFLAGS)

13.heat_gaspi_ompss2_tasks.bin: $(GASPI_SRC) src/gaspi/13.solver_gaspi_ompss2_tasks.c
	$(MCC_WRAPPERS) $(MPICC) $(CPPFLAGS) $(GASPI_CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(GASPI_LDFLAGS)

14.heat_tagaspi_ompss2_tasks.bin: $(GASPI_SRC) src/gaspi/14.solver_tagaspi_ompss2_tasks.c
	$(MCC_WRAPPERS) $(MPICC) $(CPPFLAGS) $(TAGASPI_CPPFLAGS) $(GASPI_EXT_CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(TAGASPI_LDFLAGS) $(GASPI_EXT_LDFLAGS)

15.heat_mpirma_ompss2_tasks.bin: $(MPI_SRC) src/mpi/15.solver_mpirma_ompss2_tasks.c
	$(MCC_WRAPPERS) $(MPICC) -DMPIRMA $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

16.heat_tampirma_ompss2_tasks.bin: $(MPI_SRC) src/mpi/16.solver_tampirma_ompss2_tasks.c
	$(MCC_WRAPPERS) $(MPICC) -DMPIRMA $(CPPFLAGS) $(TAMPI_CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(TAMPI_LDFLAGS)

17.heat_tampirma_ompss2_tasks.bin: $(MPI_SRC) src/mpi/17.solver_tampirma_ompss2_tasks.c
	$(MCC_WRAPPERS) $(MPICC) -DMPIRMA $(CPPFLAGS) $(TAMPI_CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(TAMPI_LDFLAGS)

18.heat_mpirma_nbuffer.bin: $(MPI_SRC) src/mpi/18.solver_mpirma_nbuffer.c
	$(CC_WRAPPERS) $(MPICC) -DMPIRMA $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

kernel.o: src/common/kernel.c
	$(CC) $(CPPFLAGS) $(KFLAGS) -c -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o *.bin
