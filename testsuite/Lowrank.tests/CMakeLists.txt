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

