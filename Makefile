MPICC ?= mpicc
HIPCC ?= hipcc
FFTW_PREFIX ?= $(CURDIR)/.local/fftw-3.3.10
HIP_PATH ?= $(DTKROOT)/hip
HOST_FLAGS = -O3 -march=native -ffp-contract=off -Wall
HIP_FLAGS = -O3 -ffp-contract=off -Wall -DMPICH_SKIP_MPICXX --offload-arch=$(PFC_GPU_ARCH) -I$(DTKROOT)/include -I$(HIP_PATH)/include
FFTW_FLAGS = -I$(FFTW_PREFIX)/include -L$(FFTW_PREFIX)/lib -Wl,-rpath,$(FFTW_PREFIX)/lib
MPI_FLAGS = -I$(I_MPI_ROOT)/include -L$(I_MPI_ROOT)/lib/release -Wl,-rpath,$(I_MPI_ROOT)/lib/release

.PHONY: all
all: bin/polygrpfc bin/polygrpfc_cpu_ref

.build bin:
	mkdir -p $@

.build/pfc.o: src/pfc.c src/pfc_dcu.h | .build
	$(MPICC) -std=c11 $(HOST_FLAGS) -DPFC_DCU $(FFTW_FLAGS) -c $< -o $@

.build/pfc_dcu.o: src/pfc_dcu.hip.cpp src/pfc_dcu.h | .build
	$(HIPCC) -std=c++14 $(HIP_FLAGS) $(MPI_FLAGS) -c $< -o $@

bin/polygrpfc: .build/pfc.o .build/pfc_dcu.o | bin
	$(HIPCC) $(HIP_FLAGS) $^ $(FFTW_FLAGS) $(MPI_FLAGS) -L$(DTKROOT)/lib -Wl,-rpath,$(DTKROOT)/lib -lfftw3_mpi -lfftw3 -lhipfft -lmpi -lm -o $@

bin/polygrpfc_cpu_ref: .build/pfc_cpu_ref.c | bin
	$(MPICC) -std=c11 $(HOST_FLAGS) $(FFTW_FLAGS) $< -lfftw3_mpi -lfftw3 -lm -o $@

.build/pfc_cpu_ref.c: | .build
	git show c86efb51cf03086fb98897ab8abfe2704a62096e:src/pfc.c > $@
