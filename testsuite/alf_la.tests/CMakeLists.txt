cmake_minimum_required(VERSION 3.5...4.2)
Project(alf_la.tests C Fortran)

# Unit tests of Libraries/Modules/alf_la_mod.F90 (ticket 202). One program per routine family; every test exits
# nonzero on a violated check. HOST backend: bit-identical to the direct BLAS/LAPACK call. CUDA backend (the flags in
# ALF_FLAGS_PROG contain -DALF_LA_CUDA): tolerances of the ticket. t_timing is a benchmark, not registered.
add_library(la_test_util STATIC la_test_util.F90)
set(alfla_tests t_checks t_gemm t_geru t_trsm t_lu t_qr t_qrpp t_geqp3 t_order t_fail t_info t_determinism)
foreach(iter ${alfla_tests} t_timing)
  add_executable(${iter} ${iter}.F90)
endforeach(iter)

foreach(iter la_test_util ${alfla_tests} t_timing)
  foreach(iflag ${FLAG_LIST})
    target_compile_options(${iter} PRIVATE ${iflag})
  endforeach(iflag)
  target_include_directories(${iter} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules/)
endforeach(iter)

# The check helpers must see NaN/Inf: compile them without finite-math assumptions (the ALF test flags contain
# -ffast-math / -fp-model=fast=2 / -fast). These options come after the flags above on the command line.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(LA_IEEE_OPT -fno-finite-math-only)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(LA_IEEE_OPT -fp-model=precise)
else()                       # NVHPC / PGI
  set(LA_IEEE_OPT -Kieee)
endif()
set_source_files_properties(la_test_util.F90 PROPERTIES COMPILE_OPTIONS "${LA_IEEE_OPT}")

foreach(iter ${alfla_tests} t_timing)
  target_link_libraries(${iter} la_test_util ${LIB_LIST})
endforeach(iter)

enable_testing()
foreach(iter ${alfla_tests})
  if(NOT iter STREQUAL "t_info")
    add_test(${iter} ${iter})
  endif()
endforeach(iter)
add_test(t_fail_exit t_fail)
add_test(t_info_recycle t_info recycle)
add_test(t_info_sync t_info sync)
set(alfla_all ${alfla_tests} t_fail_exit t_info_recycle t_info_sync)
list(REMOVE_ITEM alfla_all t_info)
foreach(iter ${alfla_all})
  set_tests_properties(${iter} PROPERTIES TIMEOUT 900)
  if("-DALF_LA_CUDA" IN_LIST FLAG_LIST OR "-DALF_LA_ONEMKL" IN_LIST FLAG_LIST)
    # fail loudly instead of silently running the target regions on the host
    set_property(TEST ${iter} APPEND PROPERTY ENVIRONMENT OMP_TARGET_OFFLOAD=MANDATORY)
  endif()
endforeach(iter)
# t_fail: must print the alf_la termination message (PASS_REGULAR_EXPRESSION) ...
set_tests_properties(t_fail PROPERTIES PASS_REGULAR_EXPRESSION "alf_la: [^\n]*(failed|status)")
# ... and, as a separate test, exit nonzero
set_tests_properties(t_fail_exit PROPERTIES WILL_FAIL TRUE)
set_tests_properties(t_info_recycle t_info_sync PROPERTIES PASS_REGULAR_EXPRESSION "alf_la: [^\n]*geqrf failed")
