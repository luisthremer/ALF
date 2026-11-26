# Tests of the GPU code path (P2.4 step 1b, 1c, 1d): gpu_state_mod (device mirror, gpu_supported, env guard, walker pack/unpack),
# gpu_wrap_mod (device hopping and Op_V kernels, self-test), gpu_stab_mod (device stabilisation, self-test, phase test), gpu_obs_mod (53), gpu_tau_mod (54).
# Registered only when the build has -DALF_GPU (configure.sh GPU=host); `make all` must have run before cmake.
# One executable, one ctest per mode (see the header of 50-gpu-state.F90); each test has its own working directory
# because Ham_set and the mirror write ./info. Skipped for TEMPERING builds (Ham_set reads Temp_<n>/ there), like
# 36-Ham-smoke.
enable_testing()
if("${FLAG_LIST}" MATCHES "-DTEMPERING")
  message(STATUS "GPU tests skipped for TEMPERING builds")
else()
  if("${FLAG_LIST}" MATCHES "-DMPI")
    find_package(MPI REQUIRED COMPONENTS Fortran)
  endif()
  set(PROG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog")
  # object list from Prog/.deps.mk (all Prog objects, gpu_state_mod.o included, except main.o)
  file(READ "${PROG_DIR}/.deps.mk" dependencies)
  string(FIND "${dependencies}" "\n\n" objects_end)
  string(SUBSTRING "${dependencies}" 0 ${objects_end} object_list)
  string(REGEX MATCHALL "[A-Za-z0-9_./+-]+\\.o" objects "${object_list}")
  list(REMOVE_ITEM objects main.o)
  set(prog_objects)
  foreach(object ${objects})
    list(APPEND prog_objects "${PROG_DIR}/${object}")
  endforeach()
  add_executable(50-gpu-state 50-gpu-state.F90)
  target_compile_options(50-gpu-state PRIVATE ${FLAG_LIST})
  target_include_directories(50-gpu-state PRIVATE "${PROG_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(50-gpu-state ${prog_objects} ${LIB_LIST})

  add_executable(51-gpu-wrap 51-gpu-wrap.F90)
  target_compile_options(51-gpu-wrap PRIVATE ${FLAG_LIST})
  target_include_directories(51-gpu-wrap PRIVATE "${PROG_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(51-gpu-wrap ${prog_objects} ${LIB_LIST})

  add_executable(52-gpu-stab 52-gpu-stab.F90)
  target_compile_options(52-gpu-stab PRIVATE ${FLAG_LIST})
  target_include_directories(52-gpu-stab PRIVATE "${PROG_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(52-gpu-stab ${prog_objects} ${LIB_LIST})

  add_executable(53-gpu-obs 53-gpu-obs.F90)
  target_compile_options(53-gpu-obs PRIVATE ${FLAG_LIST})
  target_include_directories(53-gpu-obs PRIVATE "${PROG_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(53-gpu-obs ${prog_objects} ${LIB_LIST})

  add_executable(54-gpu-tau 54-gpu-tau.F90)
  target_compile_options(54-gpu-tau PRIVATE ${FLAG_LIST})
  target_include_directories(54-gpu-tau PRIVATE "${PROG_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
  target_link_libraries(54-gpu-tau ${prog_objects} ${LIB_LIST})

  # gpu_state_test(<name> <env list as NAME=value;...> <regex or "" > <mode and arguments...>)
  function(gpu_state_test name env regex)
    set(run_dir "${CMAKE_CURRENT_BINARY_DIR}/${name}")
    file(MAKE_DIRECTORY "${run_dir}")
    if("${FLAG_LIST}" MATCHES "-DMPI")
      add_test(NAME ${name}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:50-gpu-state> ${MPIEXEC_POSTFLAGS} ${ARGN}
        WORKING_DIRECTORY "${run_dir}")
    else()
      add_test(NAME ${name} COMMAND 50-gpu-state ${ARGN} WORKING_DIRECTORY "${run_dir}")
    endif()
    set(test_env "OMP_NUM_THREADS=1")
    if(NOT "${env}" STREQUAL "")
      list(APPEND test_env ${env})
    endif()
    set_tests_properties(${name} PROPERTIES ENVIRONMENT "${test_env}" TIMEOUT 300)
    if(NOT "${regex}" STREQUAL "")
      set_tests_properties(${name} PROPERTIES PASS_REGULAR_EXPRESSION "${regex}")
    endif()
  endfunction()

  # A: table self-test (flat tables against Hop_mod_* and Op_*) on Hubbard 4x4, checkerboard and dense hopping
  gpu_state_test(gpu-state-tables "ALF_GPU_SELFTEST=1" "GPU self-test PASSED" tables)
  gpu_state_test(gpu-state-tables-dense "ALF_GPU_SELFTEST=1" "GPU self-test PASSED" tables dense)
  # fault controls: a perturbed or NaN table must make the self-test fail
  foreach(fault fault fault_v nan nan_v)
    gpu_state_test(gpu-state-fault-${fault} "ALF_GPU_SELFTEST=${fault}" "GPU self-test FAILED" tables)
  endforeach()
  # B: bitwise round trip walker -> mirror -> walker (GR, udvst, udvl, udvr; D or L as the build has it)
  gpu_state_test(gpu-state-roundtrip "" "" roundtrip)
  # C: refusals, one per runtime-detectable unsupported option
  gpu_state_test(gpu-state-refusals "" "" refusals)
  # D: ALF_GPU_DISABLE, env guard
  gpu_state_test(gpu-state-disable "ALF_GPU_DISABLE=1" "" inactive "disabled by ALF_GPU_DISABLE=1")
  gpu_state_test(gpu-state-selftest-inactive "ALF_GPU_DISABLE=1;ALF_GPU_SELFTEST=1"
    "ALF_GPU_SELFTEST is set but the GPU path is inactive" plain)
  gpu_state_test(gpu-state-env-refused "ALF_GPU_NOUPDATE=1" "GPU debug env refused: ALF_GPU_NOUPDATE" plain)
  gpu_state_test(gpu-state-fault-value "ALF_GPU_DEBUG=1;ALF_GPU_FAULT=bogus" "GPU fault value invalid: ALF_GPU_FAULT=bogus" plain)
  # ALF_GPU_CGR_LU (CGR phase source) is a debug variable with the values 0, 1, 2
  gpu_state_test(gpu-state-cgr-lu-refused "ALF_GPU_CGR_LU=0" "GPU debug env refused: ALF_GPU_CGR_LU" plain)
  gpu_state_test(gpu-state-cgr-lu-value "ALF_GPU_DEBUG=1;ALF_GPU_CGR_LU=3" "GPU env value invalid: ALF_GPU_CGR_LU=3" plain)
  gpu_state_test(gpu-state-envreader "ALF_GPU_DISABLE=1;ALF_GPU_FOO=1" "" envreader)

  # F: device wrap kernels (gpu_wrap_mod self-test run by Setup_gpu with ALF_GPU_SELFTEST=2), 51-gpu-wrap
  function(gpu_wrap_test name env regex)
    set(run_dir "${CMAKE_CURRENT_BINARY_DIR}/${name}")
    file(MAKE_DIRECTORY "${run_dir}")
    if("${FLAG_LIST}" MATCHES "-DMPI")
      add_test(NAME ${name}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:51-gpu-wrap> ${MPIEXEC_POSTFLAGS} ${ARGN}
        WORKING_DIRECTORY "${run_dir}")
    else()
      add_test(NAME ${name} COMMAND 51-gpu-wrap ${ARGN} WORKING_DIRECTORY "${run_dir}")
    endif()
    set_tests_properties(${name} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;${env}" TIMEOUT 600
      PASS_REGULAR_EXPRESSION "${regex}")
  endfunction()
  gpu_wrap_test(gpu-wrap-selftest "ALF_GPU_SELFTEST=2" "GPU wrap self-test PASSED" selftest)
  gpu_wrap_test(gpu-wrap-selftest-dense "ALF_GPU_SELFTEST=2" "GPU wrap self-test PASSED" selftest dense)
  # the opt-in dense hopping product Hd (ALF_GPU_HOPD=1, KNOWN_ISSUES G2) stays tested
  gpu_wrap_test(gpu-wrap-selftest-hopd "ALF_GPU_SELFTEST=2;ALF_GPU_HOPD=1" "GPU wrap self-test PASSED" selftest)
  # fault controls: a perturbed or NaN device copy must make the wrap self-test fail. (ALF_GPU_FAULT=hopc needs complex
  # hopping elements and =order a non-palindromic element order; the Hubbard checkerboard has neither.)
  foreach(fault hop vexp nan)
    gpu_wrap_test(gpu-wrap-fault-${fault} "ALF_GPU_SELFTEST=2;ALF_GPU_FAULT=${fault}" "GPU self-test FAILED" selftest)
  endforeach()

  # G: device stabilisation (gpu_stab_mod self-test run by the test executable with ALF_GPU_SELFTEST=3), 52-gpu-stab
  function(gpu_stab_test name env regex)
    set(run_dir "${CMAKE_CURRENT_BINARY_DIR}/${name}")
    file(MAKE_DIRECTORY "${run_dir}")
    if("${FLAG_LIST}" MATCHES "-DMPI")
      add_test(NAME ${name}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:52-gpu-stab> ${MPIEXEC_POSTFLAGS} ${ARGN}
        WORKING_DIRECTORY "${run_dir}")
    else()
      add_test(NAME ${name} COMMAND 52-gpu-stab ${ARGN} WORKING_DIRECTORY "${run_dir}")
    endif()
    set_tests_properties(${name} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;${env}" TIMEOUT 900
      PASS_REGULAR_EXPRESSION "${regex}")
  endfunction()
  gpu_stab_test(gpu-stab-selftest "ALF_GPU_SELFTEST=3" "GPU stab self-test PASSED" selftest)
  gpu_stab_test(gpu-stab-selftest-nwrap7 "ALF_GPU_SELFTEST=3" "GPU stab self-test PASSED" selftest 7)
  # fault controls: a wrong pivot permutation, a flipped CGR phase, a NaN scale or a flipped tracked det(U) phase must make the
  # self-test fail (the self-test runs in check mode: tracked phase used and compared with the LU phase, whatever the default is)
  foreach(fault perm phase nan det)
    gpu_stab_test(gpu-stab-fault-${fault} "ALF_GPU_SELFTEST=3;ALF_GPU_FAULT_STAB=${fault}" "GPU self-test FAILED" selftest)
  endforeach()
  # determinant phase of the LU factors at large N (raw product underflows), selection bitwise vs LU_PHASE_SAFE
  gpu_stab_test(gpu-stab-phase "" "gpu-stab test passed: phase" phase)
  # tracked det(U) (speedup S1): the reflector formula against the LU phase (N = 4 .. 1000), unitary, tau = 0, complex last reflector,
  # negative controls; n = 3200 (raw product underflows); an untracked slot stops CGR in the modes that use the tracked phase.
  # Every test names its mode (ALF_GPU_CGR_LU), so the results do not depend on the default of the build.
  gpu_stab_test(gpu-stab-refdet "" "gpu-stab test passed: refdet" refdet)
  gpu_stab_test(gpu-stab-refdet-3200 "" "gpu-stab test passed: refdet3200" refdet3200)
  gpu_stab_test(gpu-stab-untracked-tracked "ALF_GPU_DEBUG=1;ALF_GPU_CGR_LU=0" "untracked U in CGR" untracked)
  gpu_stab_test(gpu-stab-untracked-check "ALF_GPU_DEBUG=1;ALF_GPU_CGR_LU=2" "untracked U in CGR" untracked)
  gpu_stab_test(gpu-stab-untracked-lu "ALF_GPU_DEBUG=1;ALF_GPU_CGR_LU=1" "gpu-stab test passed: untracked" untracked)

  # H: device observables (gpu_obs_mod), 53-gpu-obs: random G through host Obser/OBSERT and the device, every observable to 1e-12
  function(gpu_obs_test name env regex)
    set(run_dir "${CMAKE_CURRENT_BINARY_DIR}/${name}")
    file(MAKE_DIRECTORY "${run_dir}")
    if("${FLAG_LIST}" MATCHES "-DMPI")
      add_test(NAME ${name}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:53-gpu-obs> ${MPIEXEC_POSTFLAGS} ${ARGN}
        WORKING_DIRECTORY "${run_dir}")
    else()
      add_test(NAME ${name} COMMAND 53-gpu-obs ${ARGN} WORKING_DIRECTORY "${run_dir}")
    endif()
    set_tests_properties(${name} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;${env}" TIMEOUT 900
      PASS_REGULAR_EXPRESSION "${regex}")
  endfunction()
  foreach(obs_case sq1 sq1s sq2 sq2s sq2h bil1 bil2)
    gpu_obs_test(gpu-obs-${obs_case} "" "gpu-obs test passed: eq ${obs_case}" eq ${obs_case})
  endforeach()
  # negative controls: device G + 1e-6, weights U on every orbital (both must exceed the tolerance by 1e3), and the fault
  # ALF_GPU_FAULT=obs (device weight of U x (1 + 1e-3)) must make the comparison fail
  gpu_obs_test(gpu-obs-control-g "" "gpu-obs test passed: control 1" control 1)
  gpu_obs_test(gpu-obs-control-weights "" "gpu-obs test passed: control 2" control 2)
  gpu_obs_test(gpu-obs-fault-obs "ALF_GPU_DEBUG=1;ALF_GPU_FAULT=obs" "gpu-obs test FAILED" eq sq2s)

  # I: device tau pass (gpu_tau_mod self-test run by the test executable with ALF_GPU_SELFTEST=4), 54-gpu-tau
  function(gpu_tau_test name env regex)
    set(run_dir "${CMAKE_CURRENT_BINARY_DIR}/${name}")
    file(MAKE_DIRECTORY "${run_dir}")
    if("${FLAG_LIST}" MATCHES "-DMPI")
      add_test(NAME ${name}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:54-gpu-tau> ${MPIEXEC_POSTFLAGS} ${ARGN}
        WORKING_DIRECTORY "${run_dir}")
    else()
      add_test(NAME ${name} COMMAND 54-gpu-tau ${ARGN} WORKING_DIRECTORY "${run_dir}")
    endif()
    set_tests_properties(${name} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;${env}" TIMEOUT 900
      PASS_REGULAR_EXPRESSION "${regex}")
  endfunction()
  gpu_tau_test(gpu-tau-selftest-mz "ALF_GPU_SELFTEST=4" "GPU tau self-test PASSED" selftest 5 mz)
  gpu_tau_test(gpu-tau-selftest-su2 "ALF_GPU_SELFTEST=4" "GPU tau self-test PASSED" selftest 5 su2)
  gpu_tau_test(gpu-tau-selftest-nwrap7 "ALF_GPU_SELFTEST=4" "GPU tau self-test PASSED" selftest 7 mz)
  # fault controls: a wrong permutation, propagation, precision reduction or block must make the self-test fail
  foreach(fault perm prop prec blocks)
    gpu_tau_test(gpu-tau-fault-${fault} "ALF_GPU_SELFTEST=4;ALF_GPU_FAULT_TAU=${fault}" "GPU self-test FAILED" selftest 5 mz)
  endforeach()

  # GPU=amd only: device path of linalg_mod against the host path (own HIP QR and hipSOLVER), mixed presence must stop
  if("${FLAG_LIST}" MATCHES "-DALF_LA_HIP")
    add_executable(55-linalg-hip 55-linalg-hip.F90)
    target_compile_options(55-linalg-hip PRIVATE ${FLAG_LIST})
    target_include_directories(55-linalg-hip PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
    target_link_libraries(55-linalg-hip ${LIB_LIST})
    add_test(NAME linalg-hip-hipqr COMMAND 55-linalg-hip compare)
    add_test(NAME linalg-hip-hipsolver COMMAND 55-linalg-hip compare)
    set_tests_properties(linalg-hip-hipsolver PROPERTIES ENVIRONMENT "ALF_HIP_QR=0")
    set_tests_properties(linalg-hip-hipqr linalg-hip-hipsolver PROPERTIES PASS_REGULAR_EXPRESSION "linalg HIP test PASSED")
    add_test(NAME linalg-hip-mixed COMMAND 55-linalg-hip mixed)
    set_tests_properties(linalg-hip-mixed PROPERTIES WILL_FAIL TRUE)
  endif()

  # E: whole ALF.out runs, Hubbard 4x4 (see gpu_e2e.py for the cases): GPU path active against the same binary with
  # ALF_GPU_DISABLE=1 (and, if ALF_GPU_REF_EXE names the ALF.out of a GPU=none build, against that): confout byte-identical,
  # observables equal to 1e-11; refusal of the timed N_delay; negative control ALF_GPU_FAULT=upd.
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_FOUND)
    set(gpu_e2e_args --exe ${PROG_DIR}/ALF.out --start ${CMAKE_CURRENT_SOURCE_DIR}/../../Scripts_and_Parameters_files/Start
        --work ${CMAKE_CURRENT_BINARY_DIR}/e2e_runs)
    if(DEFINED ENV{ALF_GPU_REF_EXE})
      list(APPEND gpu_e2e_args --ref-exe $ENV{ALF_GPU_REF_EXE})
    endif()
    if("${FLAG_LIST}" MATCHES "-DMPI")
      string(REPLACE ";" " " mpirun_string "${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};1;${MPIEXEC_PREFLAGS}")
      list(APPEND gpu_e2e_args --mpiexec "${mpirun_string}")
    endif()
    add_test(NAME gpu-state-e2e COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gpu_e2e.py ${gpu_e2e_args})
    set_tests_properties(gpu-state-e2e PROPERTIES TIMEOUT 1200)
  endif()
  # CUDA backend of linalg_mod (GPU=nvidia, needs the GPU): device wrappers vs host (55) and the stream-order canary (56)
  if("${FLAG_LIST}" MATCHES "-DALF_LA_CUDA" AND NOT "${FLAG_LIST}" MATCHES "-DMPI")
    foreach(t 55-cuda-linalg 56-cuda-order)
      add_executable(${t} ${t}.F90)
      target_compile_options(${t} PRIVATE ${FLAG_LIST})
      target_include_directories(${t} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
      target_link_libraries(${t} ${LIB_LIST})
      add_test(NAME ${t} COMMAND ${t})
    endforeach()
    add_test(NAME 55-cuda-linalg-mixed COMMAND 55-cuda-linalg mixed)
    set_tests_properties(55-cuda-linalg-mixed PROPERTIES WILL_FAIL TRUE)
  endif()
endif()
