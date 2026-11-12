# End-to-end checks of fork features with the normal ALF executable (see run_checks.py). Included from
# ../CMakeLists.txt only for a serial build without HDF5 (plain-text observables, no mpirun).
set(PROG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog")
find_program(PYTHON_EXECUTABLE NAMES python3)
if(NOT PYTHON_EXECUTABLE)
  message(WARNING "Run.tests: python3 not found, end-to-end checks not registered")
  return()
endif()
enable_testing()
foreach(scenario schemes pinning local)
  add_test(NAME run-${scenario}
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_checks.py" "${PROG_DIR}/ALF.out" ${scenario})
  set_tests_properties(run-${scenario} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 120)
endforeach()
