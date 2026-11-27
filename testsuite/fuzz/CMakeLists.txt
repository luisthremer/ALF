# P2.8 input fuzzing (fuzz_inputs.py, KNOWN_ISSUES.md): each listed case must stop with a message or run cleanly.
# Serial builds only (the harness starts ALF.out without mpiexec). HDF5 builds run the HDF5 restart cases,
# the other builds the text-file cases.
enable_testing()
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND AND NOT "${FLAG_LIST}" MATCHES "-D(MPI|TEMPERING|PARALLEL_PARAMS)")
  set(FUZZ_COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fuzz_inputs.py
      ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/ALF.out)
  if("${FLAG_LIST}" MATCHES "-DHDF5")
    set(FUZZ_HDF5 TRUE)
  else()
    set(FUZZ_HDF5 FALSE)
  endif()
  # fuzz_test(<name> <hdf5: TRUE/FALSE> <case> ...)
  function(fuzz_test name hdf5)
    if(NOT hdf5 STREQUAL FUZZ_HDF5)
      return()
    endif()
    string(REPLACE ";" "," cases "${ARGN}")
    set(mode "")
    if(hdf5)
      set(mode --hdf5)
    endif()
    add_test(NAME ${name} COMMAND ${FUZZ_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${name} ${mode} --only ${cases})
    set_tests_properties(${name} PROPERTIES PASS_REGULAR_EXPRESSION "FAILED: 0" ENVIRONMENT "OMP_NUM_THREADS=1"
                         TIMEOUT 900)
  endfunction()
endif()
