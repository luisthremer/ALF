# make BACKEND=cuda|host     (nvfortran 26.9; `. ~/alf-gpu/env/nvhpc.sh` first)
# make BACKEND=onemkl [LAPACK=dispatch|shim]   (ifx + icpx + oneMKL; `. ~/alf-gpu/env/intel.sh` first, run with
#                                               ./run_onemkl.sh, see README)
BACKEND ?= cuda
LAPACK  ?= dispatch
FC      = nvfortran
BLAS    = -lblas_seq -llapack_seq
CXXOBJ  =
EXTRALIB =
ifeq ($(BACKEND),cuda)
  FFLAGS = -O2 -g -traceback -mp=gpu -gpu=cc80 -cudalib=cublas,cusolver -DBENCH_LA_CUDA
else ifeq ($(BACKEND),host)
  FFLAGS = -O2 -g -traceback -mp -DBENCH_LA_HOST
else ifeq ($(BACKEND),onemkl)
  # JIT (spir64) only: no ocloc on nividia, so no spir64_gen AOT.  -qmkl adds the oneMKL SYCL libraries
  # (libmkl_sycl_*) automatically when -fopenmp-targets is given (verified with ldd).
  FC     = ifx
  CXX    = icpx
  MKLINC = $(MKLROOT)/include
  OFFLOAD = -fiopenmp -fopenmp-targets=spir64
  FFLAGS = -O2 -g -traceback $(OFFLOAD) -qmkl -I$(MKLINC) -DBENCH_LA_ONEMKL
  BLAS   =
  ifeq ($(LAPACK),shim)
    FFLAGS += -DBENCH_LA_SHIM
    CXXOBJ  = $(BUILD)/bench_sycl_shim.o
    EXTRALIB = -lmkl_sycl_lapack -lsycl -lstdc++
  else ifneq ($(LAPACK),dispatch)
    $(error LAPACK must be dispatch or shim)
  endif
else
  $(error BACKEND must be cuda, host or onemkl)
endif
ifeq ($(BACKEND),onemkl)
  BUILD = build/onemkl_$(LAPACK)
else
  BUILD = build/$(BACKEND)
endif
SRC   = src/bench_la.F90 src/bench_kern.F90 src/bench_util.F90 src/bench_ref.F90 src/bench_sec.F90 src/bench_main.F90
OBJ   = $(patsubst src/%.F90,$(BUILD)/%.o,$(SRC)) $(CXXOBJ)
EXE   = $(BUILD)/bench

all: $(EXE)
$(BUILD):
	mkdir -p $(BUILD)
$(BUILD)/%.o: src/%.F90 | $(BUILD)
	$(FC) $(FFLAGS) -module $(BUILD) -c $< -o $@
$(EXE): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@ $(BLAS) $(EXTRALIB)
# C++ SYCL shim (icpx -fsycl).  It contains no device kernels (only oneMKL calls), so the object is a plain
# host object and can be linked by ifx.
$(BUILD)/bench_sycl_shim.o: src/bench_sycl_shim.cpp | $(BUILD)
	$(CXX) -fsycl -O2 -g -std=c++17 -I$(MKLINC) -c $< -o $@
# standalone shim check against host LAPACK (make BACKEND=onemkl shimtest)
SHIMTEST = build/shim_test
shimtest: $(SHIMTEST)
build/bench_sycl_shim.o: src/bench_sycl_shim.cpp
	mkdir -p build
	$(CXX) -fsycl -O2 -g -std=c++17 -I$(MKLROOT)/include -c $< -o $@
build/shim_test.o: src/shim_test.F90
	mkdir -p build
	ifx -O2 -g -traceback -qmkl -c $< -o $@
$(SHIMTEST): build/shim_test.o build/bench_sycl_shim.o
	ifx -O2 -g -qmkl $^ -o $@ -lmkl_sycl_lapack -lsycl -lstdc++
# module dependencies (compile order)
$(BUILD)/bench_kern.o: $(BUILD)/bench_la.o
$(BUILD)/bench_util.o: $(BUILD)/bench_la.o
$(BUILD)/bench_ref.o:  $(BUILD)/bench_la.o
$(BUILD)/bench_sec.o:  $(BUILD)/bench_la.o $(BUILD)/bench_kern.o $(BUILD)/bench_util.o $(BUILD)/bench_ref.o
$(BUILD)/bench_main.o: $(BUILD)/bench_sec.o
clean:
	rm -rf build
.PHONY: all clean shimtest
