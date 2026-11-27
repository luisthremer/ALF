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

  # inputs that stopped (or ran) correctly before the P2.8 fixes: keep them that way
  fuzz_test(54-fuzz-baseline FALSE ok_base ok_crlf ok_tabs par_missing_file par_empty_file par_no_VAR_QMC
            par_no_VAR_lattice par_no_VAR_Model_Generic par_no_VAR_Hubbard par_no_VAR_ham_name par_truncated_in_QMC
            par_no_slash_QMC par_unknown_key_QMC par_unknown_key_Hubbard par_L1_real par_L1_huge par_Mz_int
            par_ham_name_unknown par_dup_key par_Beta_nan par_Dtau_nan par_Dtau_1e-300 par_Theta_huge_proj par_Ltau_2
            par_N_SUN_neg par_N_FL_0 par_L1_neg seeds_missing seeds_empty seeds_garbage seeds_huge seeds_zero
            confin_empty confin_truncated confin_garbage confin_field_0 confin_field_3 confin_field_neg2
            cont_ok_restart cont_field_nan confin_and_h5 size_L1_1 size_chain_L2_1 size_Ltrot_1 size_Ltrot_2e9)
  fuzz_test(54-fuzz-baseline-hdf5 TRUE h5_ok_restart h5_empty h5_no_configuration)
endif()
