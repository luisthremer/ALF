# Tests of the delayed and submatrix sequential updates (lowrank_update_mod). Run `make all` before cmake.

# 36-lowrank-unit: kernel-level harness. It drives the real Upgrade2 / Op_Wrapup / Op_Wrapdo and the
# Lowrank_update_mod kernels against an exactly recomputed Green function, so it links every Prog object
# except main.o (Prog/*.o) plus the Hamiltonian submodules.
add_executable(36-lowrank-unit 36-lowrank-unit.F90)
file(GLOB LOWRANK_UNIT_PROG_OBJS ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/*.o)
list(FILTER LOWRANK_UNIT_PROG_OBJS EXCLUDE REGEX ".*/main\\.o$")
file(GLOB LOWRANK_UNIT_HAM_OBJS ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/Hamiltonians/*.o)
target_link_libraries(36-lowrank-unit ${LOWRANK_UNIT_PROG_OBJS} ${LOWRANK_UNIT_HAM_OBJS} ${LIB_LIST})
foreach(iflag ${FLAG_LIST})
  target_compile_options(36-lowrank-unit PRIVATE ${iflag})
endforeach(iflag)
target_include_directories(36-lowrank-unit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/)
target_include_directories(36-lowrank-unit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules/)
enable_testing()

# MPI builds must be started through the MPI launcher (a bare singleton start needs the daemon of the MPI
# installation, which is missing in some containers). As root Open MPI needs two extra flags.
set(LOWRANK_MPIRUN "")
set(LOWRANK_BUILD serial)
if("${FLAG_LIST}" MATCHES "-D(MPI|TEMPERING|PARALLEL_PARAMS)")
  find_program(LOWRANK_MPIEXEC NAMES mpiexec mpirun REQUIRED)
  set(LOWRANK_MPIRUN ${LOWRANK_MPIEXEC})
  execute_process(COMMAND ${LOWRANK_MPIEXEC} --version OUTPUT_VARIABLE LOWRANK_MPI_VERSION ERROR_QUIET)
  if("${LOWRANK_MPI_VERSION}" MATCHES "Open MPI|OpenRTE")
    list(APPEND LOWRANK_MPIRUN --allow-run-as-root --oversubscribe)
  endif()
  set(LOWRANK_BUILD mpi)
  if("${FLAG_LIST}" MATCHES "-DTEMPERING")
    set(LOWRANK_BUILD tempering)
  endif()
  if("${FLAG_LIST}" MATCHES "-DPARALLEL_PARAMS")
    set(LOWRANK_BUILD parallel_params)
  endif()
endif()
if(LOWRANK_MPIRUN)
  set(LOWRANK_UNIT_COMMAND ${LOWRANK_MPIRUN} -np 1 $<TARGET_FILE:36-lowrank-unit>)
else()
  set(LOWRANK_UNIT_COMMAND 36-lowrank-unit)
endif()
add_test(NAME 36-lowrank-unit COMMAND ${LOWRANK_UNIT_COMMAND})
# Lowrank_set_N_delay must terminate with pending accepts / above the allocated capacity (timed automatic N_delay).
# They pass only on the specific error message (not on any non-zero exit). Own working directories: Lowrank_init
# appends to ./info, which 36-lowrank-unit checks and removes.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/36-pending ${CMAKE_CURRENT_BINARY_DIR}/36-capacity)
add_test(NAME 36-lowrank-unit-pending COMMAND ${LOWRANK_UNIT_COMMAND} pending WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/36-pending)
set_tests_properties(36-lowrank-unit-pending PROPERTIES
  PASS_REGULAR_EXPRESSION "Lowrank_set_N_delay: 1 accepted vertices pending \\(G must be flushed first\\)")
add_test(NAME 36-lowrank-unit-capacity COMMAND ${LOWRANK_UNIT_COMMAND} capacity WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/36-capacity)
set_tests_properties(36-lowrank-unit-capacity PROPERTIES
  PASS_REGULAR_EXPRESSION "Lowrank_set_N_delay: N_delay = 3 outside 1 \\.\\. 2 \\(allocated capacity\\)")

# Delayed and submatrix against the standard update on fixed seeds, whole runs (compare_schemes.py): the
# final field configuration must be byte-identical, the observables agree to 1e-8, 2 bins per scenario.
# lowrank-compare-schemes: a representative subset (every model and the ALF features, plus the scenarios that
# need MPI ranks in the MPI builds); lowrank-compare-schemes-full: every scenario (about 30 minutes,
# only with -DALF_LOWRANK_FULL_TESTS=ON).
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  set(LOWRANK_COMPARE_MPI "")
  if(LOWRANK_MPIRUN)
    string(REPLACE ";" " " LOWRANK_MPIRUN_STRING "${LOWRANK_MPIRUN}")
    set(LOWRANK_COMPARE_MPI --build ${LOWRANK_BUILD} --mpiexec "${LOWRANK_MPIRUN_STRING}")
  endif()
  # builds that contain the test-only Hamiltonians (add_test_models.sh) also run the scenarios that need them
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/Hamiltonians/Hamiltonian_Hubbard_Mixed_smod.F90)
    list(APPEND LOWRANK_COMPARE_MPI --test-models)
  endif()
  set(LOWRANK_COMPARE_COMMON ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_schemes.py
      ${LOWRANK_COMPARE_MPI}
      --exe ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/ALF.out
      --start ${CMAKE_CURRENT_SOURCE_DIR}/../../Scripts_and_Parameters_files/Start)
  # The quick set takes one bin per run (and, with MPI, only the automatic N_delay): the slow compilers (Intel
  # with run-time checks) need an hour otherwise. The full set below keeps 2 bins and N_delay = 0 and 3.
  set(LOWRANK_QUICK_OPTIONS --quick --max-bins 1)
  if(LOWRANK_MPIRUN)
    list(APPEND LOWRANK_QUICK_OPTIONS --n-delay 0)
  endif()
  add_test(NAME lowrank-compare-schemes
           COMMAND ${LOWRANK_COMPARE_COMMON} ${LOWRANK_QUICK_OPTIONS} --work ${CMAKE_CURRENT_BINARY_DIR}/compare_schemes_runs)
  # Global_tau_moves with a sequential range that leaves vertices out (KNOWN_ISSUES.md, H5), all three schemes.
  # One parameter set: not in the tempering / PARALLEL_PARAMS builds (they need an even number of replicas).
  if(LOWRANK_BUILD STREQUAL "serial" OR LOWRANK_BUILD STREQUAL "mpi")
    set(GLOBAL_TAU_RANGE_MPI "")
    if(LOWRANK_MPIRUN)
      set(GLOBAL_TAU_RANGE_MPI --mpiexec "${LOWRANK_MPIRUN_STRING}")
    endif()
    add_test(NAME 52-global-tau-range
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/global_tau_range.py ${GLOBAL_TAU_RANGE_MPI}
                     --exe ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/ALF.out
                     --start ${CMAKE_CURRENT_SOURCE_DIR}/../../Scripts_and_Parameters_files/Start
                     --work ${CMAKE_CURRENT_BINARY_DIR}/global_tau_range_runs)
    set_tests_properties(52-global-tau-range PROPERTIES PASS_REGULAR_EXPRESSION "PASS: global tau range" TIMEOUT 600)
  endif()
  option(ALF_LOWRANK_FULL_TESTS "Run every delayed/submatrix comparison scenario (about 30 minutes)" OFF)
  if(ALF_LOWRANK_FULL_TESTS)
    add_test(NAME lowrank-compare-schemes-full
             COMMAND ${LOWRANK_COMPARE_COMMON} --max-bins 2 --work ${CMAKE_CURRENT_BINARY_DIR}/compare_schemes_runs_full)
    set_tests_properties(lowrank-compare-schemes-full PROPERTIES TIMEOUT 7200)
  endif()
endif()
