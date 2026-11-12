cmake_minimum_required(VERSION 3.5...4.2)
Project(Analysis.tests C Fortran)

add_executable(1-maxent-kernel-overflow 1-maxent-kernel-overflow.F90)

set(analysistests 1-maxent-kernel-overflow)

foreach(iter ${analysistests})
  foreach(iflag ${FLAG_LIST})
    target_compile_options(${iter} PRIVATE ${iflag})
  endforeach(iflag)
  target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Analysis/maxent_wrapper_mod.o)
  target_link_libraries(${iter} ${LIB_LIST})
  target_include_directories(${iter} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Analysis/)
  target_include_directories(${iter} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules/)
endforeach(iter)

enable_testing()
add_test(1-maxent-kernel-overflow 1-maxent-kernel-overflow)

# End-to-end test of the Gaussian bootstrap in Max_SAC (needs a built Analysis/Max_SAC.out).
# See README_bootstrap_test.md.
find_program(PYTHON3_EXE python3)
if(PYTHON3_EXE)
  add_test(NAME 2-maxsac-bootstrap
           COMMAND ${PYTHON3_EXE} ${CMAKE_CURRENT_SOURCE_DIR}/2-maxsac-bootstrap.py
                   ${CMAKE_CURRENT_SOURCE_DIR}/../../Analysis/Max_SAC.out)
endif()
