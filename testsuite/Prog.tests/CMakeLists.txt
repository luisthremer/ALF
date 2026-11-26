cmake_minimum_required(VERSION 3.5...4.2)
Project(matmod.tests C Fortran)

# add_executable(1-copy-sel-rows 1-copy-sel-rows.F90)
# add_executable(2-copy-sel-columns 2-copy-sel-columns.F90)
# add_executable(3-opmult 3-opmult.F90)
# add_executable(4-opmultct 4-opmultct.F90)
# add_executable(6-opexpmult-1 6-opexpmult-1.F90)
# add_executable(6-opexpmult-2 6-opexpmult-2.F90)
# add_executable(7-opexpmultct-1 7-opexpmultct-1.F90)
# add_executable(7-opexpmultct-2 7-opexpmultct-2.F90)
#add_executable(8-scalematrix-1 8-scalematrix-1.F90)
#add_executable(8-scalematrix-2 8-scalematrix-2.F90)
add_executable(9-Op-Phase 9-Op-Phase.F90)
add_executable(10-Op-mmultL 10-Op-mmultL.F90)
add_executable(11-Op-mmultR 11-Op-mmultR.F90)
add_executable(13-Op-Wrapup 13-Op-Wrapup.F90)
add_executable(14-Op-Wrapdo 14-Op-Wrapdo.F90)
add_executable(15-cgr 15-cgr.F90)
add_executable(16-get-blocks 16-get-blocks.F90)
add_executable(17-solve-extended-system 17-solve-extended-system.F90)
add_executable(18-ul-update-matrices 18-ul-update-matrices.F90)
add_executable(19-ur-update-matrices 19-ur-update-matrices.F90)
add_executable(20-qdrp 20-qdrp.F90)
add_executable(21-Op-make 21-Op-make.F90)
add_executable(22-sl-mat-mults 21-sl-mat-mults.F90)
add_executable(23-cgrp 23-cgrp.F90)
add_executable(24-udv 24-udv.F90)
add_executable(25-assign-UDV-state 25-assign-UDV-state.F90)
add_executable(26-Test-Polymorphic-Fortran 26-Test-Polymorphic-Fortran.F90)
add_executable(27-Op-mmultL-m1.F90 27-Op-mmultL-m1.F90)
add_executable(28-Op-validate 28-Op-validate.F90)
add_executable(29-Op-validate-hermitian-fail 29-Op-validate-hermitian-fail.F90)
add_executable(30-Op-validate-type-fail 30-Op-validate-type-fail.F90)
add_executable(31-Op-validate-projector-fail 31-Op-validate-projector-fail.F90)
add_executable(32-Op-validate-1x1-fail 32-Op-validate-1x1-fail.F90)
add_executable(33-Op-equal 33-Op-equal.F90)
add_executable(34-Op-V-symm-pass 34-Op-V-symm-pass.F90)
add_executable(35-Op-V-symm-fail 35-Op-V-symm-fail.F90)
# set(operatorparts 1-copy-sel-rows 2-copy-sel-columns 3-opmult 4-opmultct 5-FillExpOps 6-opexpmult-1 6-opexpmult-2 7-opexpmultct-1 7-opexpmultct-2 9-Op-Phase 10-Op-mmultL 11-Op-mmultR 13-Op-Wrapup 14-Op-Wrapdo 21-Op-make 22-sl-mat-mults)
set(operatorparts 9-Op-Phase 10-Op-mmultL 11-Op-mmultR 13-Op-Wrapup 14-Op-Wrapdo 21-Op-make 22-sl-mat-mults 26-Test-Polymorphic-Fortran 27-Op-mmultL-m1.F90 28-Op-validate 29-Op-validate-hermitian-fail 30-Op-validate-type-fail 31-Op-validate-projector-fail 32-Op-validate-1x1-fail 33-Op-equal 34-Op-V-symm-pass 35-Op-V-symm-fail)
set(cgr2_2tests 16-get-blocks 17-solve-extended-system)
#set(cgr2_1tests 8-scalematrix-1 8-scalematrix-2)
set(cgr1tests 15-cgr 20-qdrp 23-cgrp)
set(udvstatetests 18-ul-update-matrices 19-ur-update-matrices 24-udv 25-assign-UDV-state)
set(progtests ${operatorparts} ${cgr1tests} ${cgr2_1tests} ${cgr2_2tests} ${udvstatetests})


foreach(iter ${operatorparts})
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/Operator_mod.o ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/Fields_mod.o)
endforeach(iter)

foreach(iter ${cgr2_2tests})
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/control_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/cgr2_2_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/UDV_WRAP_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/QDRP_decompose_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/udv_state_mod.o)
endforeach(iter)

#foreach(iter ${cgr2_1tests})
#target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/cgr2_1.o)
#target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/UDV_WRAP.o)
#target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/cgr1.o)
#endforeach(iter)

foreach(iter ${cgr1tests})
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/cgr1_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/control_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/UDV_WRAP_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/QDRP_decompose_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/udv_state_mod.o)
endforeach(iter)

foreach(iter ${udvstatetests})
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/control_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/UDV_WRAP_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/udv_state_mod.o)
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/QDRP_decompose_mod.o)
endforeach(iter)


foreach(iter ${progtests})
  foreach(iflag ${FLAG_LIST})
    target_compile_options(${iter} PRIVATE ${iflag})
  endforeach(iflag)
  target_include_directories(${iter} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/)
  target_include_directories(${iter} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules/)
  target_link_libraries(${iter} ${LIB_LIST})
endforeach(iter)


target_include_directories(26-Test-Polymorphic-Fortran PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/)
target_link_libraries(26-Test-Polymorphic-Fortran ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/OpTTypes_mod.o ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/DynamicMatrixArray_mod.o ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/ContainerElementBase_mod.o)
foreach(iflag ${FLAG_LIST})
  target_compile_options(26-Test-Polymorphic-Fortran PRIVATE ${iflag})
endforeach(iflag)

enable_testing()
# add_test(1-copy-sel-rows 1-copy-sel-rows)
# add_test(2-copy-sel-columns 2-copy-sel-columns)
# add_test(3-opmult 3-opmult)
# add_test(4-opmultct 4-opmultct)
# add_test(6-opexpmult-1 6-opexpmult-1)
# add_test(6-opexpmult-2 6-opexpmult-2)
# add_test(7-opexpmultct-1 7-opexpmultct-1)
# add_test(7-opexpmultct-2 7-opexpmultct-2)
#add_test(8-scalematrix-1 8-scalematrix-1)
#add_test(8-scalematrix-2 8-scalematrix-2)
add_test(9-Op-Phase 9-Op-Phase)
add_test(10-Op-mmultL 10-Op-mmultL)
add_test(11-Op-mmultR 11-Op-mmultR)
add_test(13-Op-Wrapup 13-Op-Wrapup)
#set_tests_properties(13-Op-Wrapup PROPERTIES WILL_FAIL TRUE)
add_test(14-Op-Wrapdo 14-Op-Wrapdo)
add_test(15-cgr 15-cgr)
add_test(16-get-blocks 16-get-blocks)
add_test(17-solve-extended-system 17-solve-extended-system)
add_test(18-ul-update-matrices 18-ul-update-matrices)
add_test(19-ur-update-matrices 19-ur-update-matrices)
add_test(20-qdrp 20-qdrp)
add_test(21-Op-make 21-Op-make)
add_test(22-sl-mat-mults 22-sl-mat-mults)
add_test(23-cgrp 23-cgrp)
add_test(24-udv 24-udv)
add_test(25-assign-UDV-state 25-assign-UDV-state)
add_test(26-Test-Polymorphic-Fortran 26-Test-Polymorphic-Fortran)
add_test(27-Op-mmultL-m1.F90 27-Op-mmultL-m1.F90)
add_test(28-Op-validate 28-Op-validate)
add_test(29-Op-validate-hermitian-fail 29-Op-validate-hermitian-fail)
set_tests_properties(29-Op-validate-hermitian-fail PROPERTIES WILL_FAIL TRUE)
add_test(30-Op-validate-type-fail 30-Op-validate-type-fail)
set_tests_properties(30-Op-validate-type-fail PROPERTIES WILL_FAIL TRUE)
add_test(31-Op-validate-projector-fail 31-Op-validate-projector-fail)
set_tests_properties(31-Op-validate-projector-fail PROPERTIES WILL_FAIL TRUE)
add_test(32-Op-validate-1x1-fail 32-Op-validate-1x1-fail)
set_tests_properties(32-Op-validate-1x1-fail PROPERTIES WILL_FAIL TRUE)
add_test(33-Op-equal 33-Op-equal)
add_test(34-Op-V-symm-pass 34-Op-V-symm-pass)
add_test(35-Op-V-symm-fail 35-Op-V-symm-fail)
set_tests_properties(35-Op-V-symm-fail PROPERTIES WILL_FAIL TRUE)

# Hamiltonian smoke test: Alloc_Ham + Ham_set + Validate_Ham_Variables with
# default parameters, one ctest per model in Prog/Hamiltonians.list (and
# Prog/Hamiltonians.list.d/*), so newly registered models are covered
# automatically. The executable links all Prog objects except main.o (object
# list from Prog/.deps.mk, as in MALA.tests). Each model runs in its own
# directory because Ham_set reads ./parameters and writes ./info. With
# TEMPERING the parameters are read from Temp_<n>/, so the test is only
# built without it.
if("${FLAG_LIST}" MATCHES "-DTEMPERING")
  message(STATUS "Hamiltonian smoke tests skipped for TEMPERING builds")
else()
  if("${FLAG_LIST}" MATCHES "-DMPI")
    find_package(MPI REQUIRED COMPONENTS Fortran)
  endif()
  set(PROG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog")
  file(STRINGS "${PROG_DIR}/Hamiltonians.list" ham_lines)
  file(GLOB ham_list_d "${PROG_DIR}/Hamiltonians.list.d/*")
  foreach(list_file ${ham_list_d})
    file(STRINGS "${list_file}" more_lines)
    list(APPEND ham_lines ${more_lines})
  endforeach()
  set(ham_names)
  foreach(line ${ham_lines})
    string(STRIP "${line}" line)
    if(line AND NOT line MATCHES "^#")
      string(REGEX MATCH "^[^ \t]+" ham_name "${line}")
      list(APPEND ham_names ${ham_name})
    endif()
  endforeach()
  file(READ "${PROG_DIR}/.deps.mk" dependencies)
  string(FIND "${dependencies}" "\n\n" objects_end)
  string(SUBSTRING "${dependencies}" 0 ${objects_end} object_list)
  string(REGEX MATCHALL "[A-Za-z0-9_./+-]+\\.o" objects "${object_list}")
  list(REMOVE_ITEM objects main.o)
  set(prog_objects)
  foreach(object ${objects})
    list(APPEND prog_objects "${PROG_DIR}/${object}")
  endforeach()
  add_executable(36-Ham-smoke 36-Ham-smoke.F90)
  target_compile_options(36-Ham-smoke PRIVATE ${FLAG_LIST})
  target_include_directories(36-Ham-smoke PRIVATE "${PROG_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(36-Ham-smoke ${prog_objects} ${LIB_LIST})
  foreach(ham_name ${ham_names})
    set(run_dir "${CMAKE_CURRENT_BINARY_DIR}/36-Ham-smoke.${ham_name}")
    file(MAKE_DIRECTORY "${run_dir}")
    if("${FLAG_LIST}" MATCHES "-DMPI")
      add_test(NAME 36-Ham-smoke-${ham_name}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:36-Ham-smoke> ${MPIEXEC_POSTFLAGS} ${ham_name}
        WORKING_DIRECTORY "${run_dir}")
    else()
      add_test(NAME 36-Ham-smoke-${ham_name} COMMAND 36-Ham-smoke ${ham_name}
        WORKING_DIRECTORY "${run_dir}")
    endif()
    set_tests_properties(36-Ham-smoke-${ham_name} PROPERTIES
      ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 60)
  endforeach()

  # Unit test for qmc_walker_mod (walker layout, allocation, deallocation).
  add_executable(40-walker 40-walker.F90)
  target_compile_options(40-walker PRIVATE ${FLAG_LIST})
  target_include_directories(40-walker PRIVATE "${PROG_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(40-walker ${prog_objects} ${LIB_LIST})
  if("${FLAG_LIST}" MATCHES "-DMPI")
    add_test(NAME 40-walker
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:40-walker> ${MPIEXEC_POSTFLAGS})
  else()
    add_test(NAME 40-walker COMMAND 40-walker)
  endif()
  # Control_Check_Phase must abort on a NaN phase (KNOWN_ISSUES.md, determinant phase underflow).
  add_executable(41-check-phase-fail 41-check-phase-fail.F90)
  target_compile_options(41-check-phase-fail PRIVATE ${FLAG_LIST})
  target_include_directories(41-check-phase-fail PRIVATE "${PROG_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(41-check-phase-fail ${prog_objects} ${LIB_LIST})
  if("${FLAG_LIST}" MATCHES "-DMPI")
    add_test(NAME 41-check-phase-fail
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:41-check-phase-fail> ${MPIEXEC_POSTFLAGS})
  else()
    add_test(NAME 41-check-phase-fail COMMAND 41-check-phase-fail)
  endif()
  set_tests_properties(41-check-phase-fail PROPERTIES
    PASS_REGULAR_EXPRESSION "phase\\| /= 1: *\\( *NaN")

  # Start-up checks outside QMC_runtime_var (KNOWN_ISSUES.md H4, M1, M3): each case must stop with its message.
  add_executable(43-prog-checks-fail 43-prog-checks-fail.F90)
  target_compile_options(43-prog-checks-fail PRIVATE ${FLAG_LIST})
  target_include_directories(43-prog-checks-fail PRIVATE "${PROG_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(43-prog-checks-fail ${prog_objects} ${LIB_LIST})
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/43-prog-checks-fail.run")
  foreach(case_msg "field0|does not fit to type" "ltrot|Dtau must be > 0" "nwrap|Nwrap must be >= 1")
    string(REPLACE "|" ";" case_msg_list "${case_msg}")
    list(GET case_msg_list 0 case_name)
    list(GET case_msg_list 1 case_regex)
    if("${FLAG_LIST}" MATCHES "-DMPI")
      add_test(NAME 43-prog-checks-fail-${case_name}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:43-prog-checks-fail> ${MPIEXEC_POSTFLAGS} ${case_name}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/43-prog-checks-fail.run")
    else()
      add_test(NAME 43-prog-checks-fail-${case_name} COMMAND 43-prog-checks-fail ${case_name}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/43-prog-checks-fail.run")
    endif()
    set_tests_properties(43-prog-checks-fail-${case_name} PROPERTIES
      PASS_REGULAR_EXPRESSION "${case_regex}" TIMEOUT 60)
  endforeach()

  # Pot_scal weights of the Hubbard bilayer lattices (KNOWN_ISSUES.md H7)
  add_executable(44-Hubbard-pot-weights 44-Hubbard-pot-weights.F90)
  target_compile_options(44-Hubbard-pot-weights PRIVATE ${FLAG_LIST})
  target_include_directories(44-Hubbard-pot-weights PRIVATE "${PROG_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(44-Hubbard-pot-weights ${prog_objects} ${LIB_LIST})
  foreach(lattice Bilayer_square Bilayer_honeycomb)
    set(run_dir "${CMAKE_CURRENT_BINARY_DIR}/44-Hubbard-pot-weights.${lattice}")
    file(MAKE_DIRECTORY "${run_dir}")
    if("${FLAG_LIST}" MATCHES "-DMPI")
      add_test(NAME 44-Hubbard-pot-weights-${lattice}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:44-Hubbard-pot-weights> ${MPIEXEC_POSTFLAGS} ${lattice}
        WORKING_DIRECTORY "${run_dir}")
    else()
      add_test(NAME 44-Hubbard-pot-weights-${lattice} COMMAND 44-Hubbard-pot-weights ${lattice}
        WORKING_DIRECTORY "${run_dir}")
    endif()
    set_tests_properties(44-Hubbard-pot-weights-${lattice} PROPERTIES
      ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 60)
  endforeach()
endif()

# Unit test for QMC_runtime_var (defaults, set/get round trip, lock).
# In MPI builds it is started with mpiexec, like 36-Ham-smoke.
if("${FLAG_LIST}" MATCHES "-DMPI")
  find_package(MPI REQUIRED COMPONENTS Fortran)
endif()
add_executable(37-QMC-runtime-var 37-QMC-runtime-var.F90)
target_compile_options(37-QMC-runtime-var PRIVATE ${FLAG_LIST})
target_include_directories(37-QMC-runtime-var PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
target_link_libraries(37-QMC-runtime-var
  "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/QMC_runtime_var_mod.o" ${LIB_LIST})
if("${FLAG_LIST}" MATCHES "-DMPI")
  add_test(NAME 37-QMC-runtime-var
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:37-QMC-runtime-var> ${MPIEXEC_POSTFLAGS})
else()
  add_test(NAME 37-QMC-runtime-var COMMAND 37-QMC-runtime-var)
endif()

# A setter after the lock must terminate the program.
add_executable(38-QMC-runtime-lock-fail 38-QMC-runtime-lock-fail.F90)
target_compile_options(38-QMC-runtime-lock-fail PRIVATE ${FLAG_LIST})
target_include_directories(38-QMC-runtime-lock-fail PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
target_link_libraries(38-QMC-runtime-lock-fail
  "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/QMC_runtime_var_mod.o" ${LIB_LIST})
if("${FLAG_LIST}" MATCHES "-DMPI")
  add_test(NAME 38-QMC-runtime-lock-fail
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:38-QMC-runtime-lock-fail> ${MPIEXEC_POSTFLAGS})
else()
  add_test(NAME 38-QMC-runtime-lock-fail COMMAND 38-QMC-runtime-lock-fail)
endif()
set_tests_properties(38-QMC-runtime-lock-fail PROPERTIES WILL_FAIL TRUE TIMEOUT 60)

# Start-up checks in QMC_runtime_var (KNOWN_ISSUES.md H3, M2, M4): each case must stop with its message.
add_executable(42-QMC-runtime-checks-fail 42-QMC-runtime-checks-fail.F90)
target_compile_options(42-QMC-runtime-checks-fail PRIVATE ${FLAG_LIST})
target_include_directories(42-QMC-runtime-checks-fail PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
target_link_libraries(42-QMC-runtime-checks-fail
  "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/QMC_runtime_var_mod.o" ${LIB_LIST})
foreach(case_msg "nsweep|NSweep must be >= 1" "lobs|LOBS_ST > LOBS_EN" "projector|does not support Ltau = 1"
    "ntseq_start|Nt_sequential_start" "ntseq_end|Nt_sequential_end")
  string(REPLACE "|" ";" case_msg_list "${case_msg}")
  list(GET case_msg_list 0 case_name)
  list(GET case_msg_list 1 case_regex)
  if("${FLAG_LIST}" MATCHES "-DMPI")
    add_test(NAME 42-QMC-runtime-checks-fail-${case_name}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:42-QMC-runtime-checks-fail> ${MPIEXEC_POSTFLAGS} ${case_name})
  else()
    add_test(NAME 42-QMC-runtime-checks-fail-${case_name} COMMAND 42-QMC-runtime-checks-fail ${case_name})
  endif()
  set_tests_properties(42-QMC-runtime-checks-fail-${case_name} PROPERTIES
    PASS_REGULAR_EXPRESSION "${case_regex}" TIMEOUT 60)
endforeach()

# Unit test for the communication layer comm_mod (Libraries/Modules).
add_executable(39-comm-mod 39-comm-mod.F90)
target_compile_options(39-comm-mod PRIVATE ${FLAG_LIST})
target_include_directories(39-comm-mod PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
target_link_libraries(39-comm-mod ${LIB_LIST})
if("${FLAG_LIST}" MATCHES "-DMPI")
  add_test(NAME 39-comm-mod
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:39-comm-mod> ${MPIEXEC_POSTFLAGS})
else()
  add_test(NAME 39-comm-mod COMMAND 39-comm-mod)
endif()

# ALF_UNSTABLE_THRESHOLD (control_mod): a deviation of 0.5 passes with the default (10) and with 1.0,
# aborts with 0.25; a non-positive or unparsable value ends the program with exit status 2 (alf_exit).
add_executable(45-unstable-threshold 45-unstable-threshold.F90)
target_compile_options(45-unstable-threshold PRIVATE ${FLAG_LIST})
target_include_directories(45-unstable-threshold PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
target_link_libraries(45-unstable-threshold "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/control_mod.o" ${LIB_LIST})
add_test(45-unstable-threshold-default 45-unstable-threshold)
set_tests_properties(45-unstable-threshold-default PROPERTIES ENVIRONMENT "ALF_UNSTABLE_THRESHOLD="
  PASS_REGULAR_EXPRESSION "PASS: deviation")
add_test(45-unstable-threshold-1.0 45-unstable-threshold)
set_tests_properties(45-unstable-threshold-1.0 PROPERTIES ENVIRONMENT "ALF_UNSTABLE_THRESHOLD=1.0"
  PASS_REGULAR_EXPRESSION "PASS: deviation")
if(NOT "${FLAG_LIST}" MATCHES "-DTEMPERING")
  # (with TEMPERING the abort branch first writes Temp_<igroup>/info, a directory a unit test does not have)
  add_test(45-unstable-threshold-0.25 45-unstable-threshold)
  set_tests_properties(45-unstable-threshold-0.25 PROPERTIES ENVIRONMENT "ALF_UNSTABLE_THRESHOLD=0.25"
    PASS_REGULAR_EXPRESSION "exceeding the threshold" FAIL_REGULAR_EXPRESSION "PASS: deviation" TIMEOUT 60)
endif()
set(EXPECT_EXIT ${CMAKE_CURRENT_SOURCE_DIR}/expect_exit.cmake)
foreach(bad abc -1 0)
  add_test(NAME 45-unstable-threshold-bad-${bad}
    COMMAND ${CMAKE_COMMAND} -DCMD=$<TARGET_FILE:45-unstable-threshold> -DEXPECTED=2
      -DENVVAR=ALF_UNSTABLE_THRESHOLD=${bad} -P ${EXPECT_EXIT})
  set_tests_properties(45-unstable-threshold-bad-${bad} PROPERTIES
    PASS_REGULAR_EXPRESSION "not a positive number" TIMEOUT 60)
endforeach()

# alf_exit (runtime_error_mod): exit status 0 / 1 / 2, and 1 for any other code.
add_executable(46-alf-exit 46-alf-exit.F90)
target_compile_options(46-alf-exit PRIVATE ${FLAG_LIST})
target_include_directories(46-alf-exit PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
target_link_libraries(46-alf-exit ${LIB_LIST})
foreach(pair "0;0" "1;1" "2;2" "7;1")
  list(GET pair 0 code)
  list(GET pair 1 expected)
  add_test(NAME 46-alf-exit-${code}
    COMMAND ${CMAKE_COMMAND} -DCMD=$<TARGET_FILE:46-alf-exit>\;${code} -DEXPECTED=${expected} -P ${EXPECT_EXIT})
  set_tests_properties(46-alf-exit-${code} PROPERTIES TIMEOUT 60)
endforeach()

# Trace_updates_mod: the trace line format of the ALF_TRACE_UPDATES debug build.
add_executable(47-trace-updates 47-trace-updates.F90)
target_compile_options(47-trace-updates PRIVATE ${FLAG_LIST})
target_include_directories(47-trace-updates PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog")
target_link_libraries(47-trace-updates "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/Trace_updates_mod.o" ${LIB_LIST})
add_test(NAME 47-trace-updates COMMAND 47-trace-updates WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(47-trace-updates PROPERTIES PASS_REGULAR_EXPRESSION "PASS: trace" TIMEOUT 60)

# Timers_mod with -DALF_TIMERS (the program build leaves it off): stack timers, counters, timers.txt.
# The module source is compiled into the test with its own module directory, so the Prog build's
# timers_mod.mod (without the macro) is not used.
add_executable(48-timers 48-timers.F90 "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/Timers_mod.F90")
target_compile_options(48-timers PRIVATE ${FLAG_LIST} -DALF_TIMERS)
set_target_properties(48-timers PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_48")
target_include_directories(48-timers PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
target_link_libraries(48-timers ${LIB_LIST})
if("${FLAG_LIST}" MATCHES "-DMPI")
  add_test(NAME 48-timers
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:48-timers> ${MPIEXEC_POSTFLAGS} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
else()
  add_test(NAME 48-timers COMMAND 48-timers WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
set_tests_properties(48-timers PROPERTIES PASS_REGULAR_EXPRESSION "PASS: timers" TIMEOUT 60)

# Prog/Makefile (KNOWN_ISSUES.md, MK): a changed .mod file of Libraries/Modules must rebuild the Prog objects.
# make -n -W pretends that runtime_error_mod.mod was just rewritten and prints the commands it would run.
find_program(ALF_MAKE NAMES make gmake)
if(ALF_MAKE)
  add_test(NAME 50-make-lib-deps
    COMMAND ${ALF_MAKE} -n -W ../Libraries/Modules/runtime_error_mod.mod ALF.out
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog)
  set_tests_properties(50-make-lib-deps PROPERTIES PASS_REGULAR_EXPRESSION "-o [A-Za-z0-9_/]+\\.o " TIMEOUT 60)
endif()

# Precision Green warning (control_mod, KNOWN_ISSUES.md M6): a deviation of 0.5 is reported, one of 1e-6 is not.
if(NOT "${FLAG_LIST}" MATCHES "-DTEMPERING")
  # (with TEMPERING Control_Print writes Temp_<igroup>/info, a directory a unit test does not have)
  add_executable(51-precision-green-warning 51-precision-green-warning.F90)
  target_compile_options(51-precision-green-warning PRIVATE ${FLAG_LIST})
  target_include_directories(51-precision-green-warning PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(51-precision-green-warning "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/control_mod.o" ${LIB_LIST})
  foreach(case large small)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/51-${case})
    if("${FLAG_LIST}" MATCHES "-DMPI")
      add_test(NAME 51-precision-green-warning-${case}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:51-precision-green-warning> ${MPIEXEC_POSTFLAGS} ${case}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/51-${case})
    else()
      add_test(NAME 51-precision-green-warning-${case} COMMAND 51-precision-green-warning ${case}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/51-${case})
    endif()
  endforeach()
  set_tests_properties(51-precision-green-warning-large PROPERTIES
    PASS_REGULAR_EXPRESSION "Warning: Precision Green Max" TIMEOUT 60)
  set_tests_properties(51-precision-green-warning-small PROPERTIES
    PASS_REGULAR_EXPRESSION "PASS: summary written" FAIL_REGULAR_EXPRESSION "Warning: Precision Green" TIMEOUT 60)
endif()
