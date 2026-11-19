# Tests of the AI feature LOWRANK (master_ai; written with AI assistance, not part of the ALF baseline).
# Added by testsuite/CMakeLists.txt only for builds with -DALF_LOWRANK. Run `make lib program` before cmake.

# 36-lowrank-unit: kernel-level harness. It drives the real Upgrade2 / Op_Wrapup / Op_Wrapdo and the
# Lowrank_update_mod kernels against an exactly recomputed Green function, so it links every Prog object
# except main.o (Prog/*.o and Prog/Lowrank/*.o) plus the Hamiltonian submodules.
add_executable(36-lowrank-unit 36-lowrank-unit.F90)
file(GLOB LOWRANK_UNIT_PROG_OBJS ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/*.o ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/Lowrank/*.o)
list(FILTER LOWRANK_UNIT_PROG_OBJS EXCLUDE REGEX ".*/main\\.o$")
file(GLOB LOWRANK_UNIT_HAM_OBJS ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/Hamiltonians/*.o)
target_link_libraries(36-lowrank-unit ${LOWRANK_UNIT_PROG_OBJS} ${LOWRANK_UNIT_HAM_OBJS} ${LIB_LIST})
foreach(iflag ${FLAG_LIST})
  target_compile_options(36-lowrank-unit PRIVATE ${iflag})
endforeach(iflag)
target_include_directories(36-lowrank-unit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/)
target_include_directories(36-lowrank-unit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules/)
enable_testing()
add_test(NAME 36-lowrank-unit COMMAND 36-lowrank-unit)

# Delayed and submatrix against the standard update on fixed seeds, whole runs (all serial scenarios of
# test_branch_parameters.json plus Global_tau, MALA, continuous-field and projector runs; 2 bins each).
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  add_test(NAME lowrank-compare-schemes
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_schemes.py
                   --exe ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/ALF.out
                   --start ${CMAKE_CURRENT_SOURCE_DIR}/../../Scripts_and_Parameters_files/Start
                   --work ${CMAKE_CURRENT_BINARY_DIR}/compare_schemes_runs --max-bins 2)
endif()
