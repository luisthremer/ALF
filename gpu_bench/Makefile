# make BACKEND=cuda|host   (nvfortran 26.9; `. ~/alf-gpu/env/nvhpc.sh` first)
BACKEND ?= cuda
FC      = nvfortran
BLAS    = -lblas_seq -llapack_seq
ifeq ($(BACKEND),cuda)
  FFLAGS = -O2 -g -traceback -mp=gpu -gpu=cc80 -cudalib=cublas,cusolver -DBENCH_LA_CUDA
else ifeq ($(BACKEND),host)
  FFLAGS = -O2 -g -traceback -mp -DBENCH_LA_HOST
else
  $(error BACKEND must be cuda or host)
endif
BUILD = build/$(BACKEND)
SRC   = src/bench_la.F90 src/bench_kern.F90 src/bench_util.F90 src/bench_ref.F90 src/bench_sec.F90 src/bench_main.F90
OBJ   = $(patsubst src/%.F90,$(BUILD)/%.o,$(SRC))
EXE   = $(BUILD)/bench

all: $(EXE)
$(BUILD):
	mkdir -p $(BUILD)
$(BUILD)/%.o: src/%.F90 | $(BUILD)
	$(FC) $(FFLAGS) -module $(BUILD) -c $< -o $@
$(EXE): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@ $(BLAS)
# module dependencies (compile order)
$(BUILD)/bench_kern.o: $(BUILD)/bench_la.o
$(BUILD)/bench_util.o: $(BUILD)/bench_la.o
$(BUILD)/bench_ref.o:  $(BUILD)/bench_la.o
$(BUILD)/bench_sec.o:  $(BUILD)/bench_la.o $(BUILD)/bench_kern.o $(BUILD)/bench_util.o $(BUILD)/bench_ref.o
$(BUILD)/bench_main.o: $(BUILD)/bench_sec.o
clean:
	rm -rf build
.PHONY: all clean
