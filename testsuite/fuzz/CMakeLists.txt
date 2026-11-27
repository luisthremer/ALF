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
    # h5fc of the HDF5 the build links (configure.sh: -L<HDF5_DIR>/lib), for h5mut.f90
    set(FUZZ_H5FC "")
    foreach(item ${LIB_LIST})
      if(item MATCHES "^-L(.*)/lib$")
        if(EXISTS "${CMAKE_MATCH_1}/bin/h5fc")
          set(FUZZ_H5FC --h5fc "${CMAKE_MATCH_1}/bin/h5fc")
        endif()
      endif()
    endforeach()
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
      set(mode --hdf5 ${FUZZ_H5FC})
    endif()
    add_test(NAME ${name} COMMAND ${FUZZ_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${name} ${mode} --only ${cases})
    set_tests_properties(${name} PROPERTIES PASS_REGULAR_EXPRESSION "FAILED: 0" ENVIRONMENT "OMP_NUM_THREADS=1"
                         TIMEOUT 900)
  endfunction()

  # inputs that stopped (or ran) correctly before the P2.8 fixes: keep them that way
  # (not par_Dtau_1e-300: nint(Beta/Dtau) traps in -ffpe-trap=invalid builds, KNOWN_ISSUES.md Z9)
  fuzz_test(54-fuzz-baseline FALSE ok_base ok_crlf ok_tabs par_missing_file par_empty_file par_no_VAR_QMC
            par_no_VAR_lattice par_no_VAR_Model_Generic par_no_VAR_Hubbard par_no_VAR_ham_name par_truncated_in_QMC
            par_no_slash_QMC par_unknown_key_QMC par_unknown_key_Hubbard par_L1_real par_L1_huge par_Mz_int
            par_ham_name_unknown par_dup_key par_Beta_nan par_Dtau_nan par_Theta_huge_proj par_Ltau_2
            par_N_SUN_neg par_N_FL_0 par_L1_neg seeds_missing seeds_empty seeds_garbage seeds_huge seeds_zero
            confin_empty confin_truncated confin_garbage confin_field_0 confin_field_3 confin_field_neg2
            cont_ok_restart cont_field_nan confin_and_h5 size_L1_1 size_chain_L2_1 size_Ltrot_1 size_Ltrot_2e9)
  fuzz_test(54-fuzz-baseline-hdf5 TRUE h5_ok_restart h5_empty h5_no_configuration)
  # NaN/Inf in a real model parameter (generated readers, parse_ham_mod.py)
  fuzz_test(55-fuzz-model-params FALSE par_U_nan par_T_nan par_PhiX_nan par_U2_nan par_chem_nan ok_base)
  # NaN/Inf in a real parameter of &VAR_QMC
  fuzz_test(56-fuzz-qmc-params FALSE par_CPU_MAX_nan par_Amplitude_inf_cont par_dtLangevin_nan par_MaxForce_inf
            ok_base)
  # no bins or no measured time slice
  fuzz_test(57-fuzz-no-measurements FALSE par_NBin_0 par_NBin_neg par_LOBS_neg par_LOBS_beyond ok_base)
  # NaN field values (confin_0.h5)
  fuzz_test(58-fuzz-field-nan TRUE h5_field_nan h5_ok_restart)
  # confin_0.h5: unreadable file, missing datasets, datasets of another shape
  fuzz_test(59-fuzz-h5-restart TRUE h5_empty h5_truncated h5_garbage h5_from_L6 h5_from_L2 h5_no_configuration
            h5_no_seed h5_short_seed h5_config_int h5_ok_restart)
  # text confin_0 with more values than fields (another lattice size or Beta)
  fuzz_test(60-fuzz-confin-size FALSE confin_from_L6 confin_extra_lines confin_from_L2 ok_base)
  # data.h5 that HDF5 cannot open
  fuzz_test(61-fuzz-data-h5 TRUE datah5_truncated datah5_garbage h5_ok_restart)
endif()
