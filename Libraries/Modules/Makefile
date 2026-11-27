.PHONY : lib clean

LIB=modules_90.a

SRCS = Mat_subroutines_mod.F90 \
	alf_hdf5_mod.F90 \
	check_omp_num_threads_mod.F90 \
	comm_mod.F90 \
	entanglement_mod.F90 \
	errors_mod.F90 \
	files_mod.F90 \
	fourier_mod.F90 \
	histograms_mod.F90 \
	histograms_v2_mod.F90 \
	iso8601_datetime_mod.F90 \
	lattices_interface_mod.F90 \
	lattices_v3_mod.F90 \
	linalg_mod.F90 \
	log_mesh_mod.F90 \
	matrix_mod.F90 \
	maxent_mod.F90 \
	maxent_stoch_mod.F90 \
	mpi_shared_mem_mod_v2.F90 \
	mymats_mod.F90 \
	natural_constants_mod.F90 \
	precdef_mod.F90 \
	random_wrap_mod.F90 \
	runtime_error_mod.F90

# GPU=intel (configure.sh): oneMKL device backend of linalg_mod
ifeq ($(strip $(ALF_GPU)),intel)
SRCS += linalg_onemkl_mod.F90
SHIMOBJ = linalg_sycl_shim.o
endif

# CUDA backend of linalg_mod (GPU=nvidia only)
ifeq ($(strip $(ALF_GPU)),nvidia)
SRCS += linalg_cuda_mod.F90
endif

# GPU=amd (configure.sh exports ALF_GPU=amd): the HIP device path of linalg_mod (Fortran layer, C++ hipBLAS/hipSOLVER shim,
# HIP kernels and QR); ROCM_PATH and ALF_AMD_ARCH come from configure.sh
ifeq ($(ALF_GPU),amd)
SRCS += linalg_hip_mod.F90 alf_hipk_mod.F90
ALF_AMD_ARCH ?= gfx90a
ROCM_INC ?= $(if $(ROCM_PATH),$(ROCM_PATH),/opt/rocm)/include
HIPCC ?= $(if $(ROCM_PATH),$(ROCM_PATH),/opt/rocm)/bin/hipcc
HIPOBJS = alf_la_hip_shim.o alf_hip_kernels.o alf_hip_qr.o
endif

# Dependency helper variables
DEPSFILE := .deps.mk
DEPSGEN := ../../Prog/gen_deps.py

lib: $(LIB) libfakeintel.so

# OBJS and MODS are auto-generated into $(DEPSFILE) by gen_deps.py
-include $(DEPSFILE)

$(LIB): $(DEPSFILE) $(OBJS) $(SHIMOBJ) $(HIPOBJS)
	ar -r $(LIB) $(OBJS) $(SHIMOBJ) $(HIPOBJS)

alf_la_hip_shim.o: alf_la_hip_shim.cpp
	g++ -O2 -std=c++17 -D__HIP_PLATFORM_AMD__ -I$(ROCM_INC) -c $< -o $@
alf_hip_kernels.o: alf_hip_kernels.hip.cpp
	$(HIPCC) -O2 --offload-arch=$(ALF_AMD_ARCH) -I$(ROCM_INC) -c $< -o $@
alf_hip_qr.o: alf_hip_qr.hip.cpp
	$(HIPCC) -O2 --offload-arch=$(ALF_AMD_ARCH) -I$(ROCM_INC) -c $< -o $@

lattices_x86-64.so: runtime_error_mod.F90 lattices_interface_mod.F90 matrix_mod.F90 lattices_v3_mod.F90
	(make clean)
	gfortran -O3 -march=x86-64 -c runtime_error_mod.F90
	gfortran -O3 -march=x86-64 -c matrix_mod.F90
	gfortran -O3 -march=x86-64 -c lattices_v3_mod.F90
	gfortran -O3 -march=x86-64 -shared -fPIC lattices_interface_mod.F90 -o lattices_x86-64.so \
			runtime_error_mod.o matrix_mod.o lattices_v3_mod.o

lattices_armv8.4-a.so: runtime_error_mod.F90 lattices_interface_mod.F90 matrix_mod.F90 lattices_v3_mod.F90
	(make clean)
	gfortran -O3 -march=armv8.4-a -c runtime_error_mod.F90
	gfortran -O3 -march=armv8.4-a -c matrix_mod.F90
	gfortran -O3 -march=armv8.4-a -c lattices_v3_mod.F90
	gfortran -O3 -march=armv8.4-a -shared -fPIC lattices_interface_mod.F90 -o lattices_armv8.4-a.so \
			runtime_error_mod.o matrix_mod.o lattices_v3_mod.o

%.o: %.F90
	$(ALF_FC) -c $(ALF_FLAGS_MODULES) $<

libfakeintel.so: fakeintel.c
	@for cc in gcc clang icc icx cc; do \
		if command -v $$cc > /dev/null 2>&1; then \
			$$cc -shared -fPIC -o libfakeintel.so fakeintel.c; \
			break; \
		fi; \
	done; \
	if [ ! -f libfakeintel.so ]; then \
		echo "No C compiler found, libfakeintel.so not created." 1>&2; \
	fi

linalg_sycl_shim.o: linalg_sycl_shim.cpp
	icpx -fsycl -O2 -std=c++17 $(ALF_SYCL_FLAGS) -I$(MKLROOT)/include -c $<

entanglement_mod.o: entanglement_mod.F90
	$(ALF_FC) -c $(shell echo $(ALF_FLAGS_MODULES) | sed 's| -pedantic||') $<

$(DEPSFILE): $(SRCS) $(DEPSGEN)
	python3 $(DEPSGEN) $(SRCS) > $(DEPSFILE)

clean: $(DEPSFILE)
	rm -f $(OBJS) $(MODS) $(LIB) $(DEPSFILE) libfakeintel.so linalg_sycl_shim.o alf_la_hip_shim.o alf_hip_kernels.o alf_hip_qr.o

#unused modules: histograms histograms_v2 log_mesh natural_constants
#unused, except for in Prog/FFA_Orginals/: precdef
