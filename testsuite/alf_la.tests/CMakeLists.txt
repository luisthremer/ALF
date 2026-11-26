cmake_minimum_required(VERSION 3.5...4.2)
Project(alf_la.tests C Fortran)

# Unit tests of Libraries/Modules/alf_la_mod.F90 (HOST backend). One program per routine family; every test exits
# nonzero on a violated check. HOST backend: bit-identical to the direct BLAS/LAPACK call.
add_library(la_test_util STATIC la_test_util.F90)
set(alfla_tests t_checks t_gemm t_geru t_trsm t_lu t_qr t_qrpp t_geqp3 t_fail t_determinism t_lwork)
foreach(iter ${alfla_tests})
  add_executable(${iter} ${iter}.F90)
endforeach(iter)

foreach(iter la_test_util ${alfla_tests})
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

foreach(iter ${alfla_tests})
  target_link_libraries(${iter} la_test_util ${LIB_LIST})
endforeach(iter)

enable_testing()
foreach(iter ${alfla_tests})
  add_test(${iter} ${iter})
  set_tests_properties(${iter} PROPERTIES TIMEOUT 900)
endforeach(iter)
# t_fail must print the alf_la termination message (PASS_REGULAR_EXPRESSION) ...
set_tests_properties(t_fail PROPERTIES PASS_REGULAR_EXPRESSION "alf_la: [^\n]*(failed|status)")
# ... and, as a separate test, exit nonzero
add_test(t_fail_exit t_fail)
set_tests_properties(t_fail_exit PROPERTIES WILL_FAIL TRUE TIMEOUT 900)
