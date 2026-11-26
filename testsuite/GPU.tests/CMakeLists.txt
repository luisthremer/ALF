# Tests of the GPU code path (P2.4 step 1b, 1c): gpu_state_mod (device mirror, gpu_supported, env guard, walker pack/unpack),
# gpu_wrap_mod (device hopping and Op_V kernels, self-test).
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
  # fault controls: a perturbed or NaN device copy must make the wrap self-test fail. (ALF_GPU_FAULT=hopc needs complex
  # hopping elements and =order a non-palindromic element order; the Hubbard checkerboard has neither.)
  foreach(fault hop vexp nan)
    gpu_wrap_test(gpu-wrap-fault-${fault} "ALF_GPU_SELFTEST=2;ALF_GPU_FAULT=${fault}" "GPU self-test FAILED" selftest)
  endforeach()

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
endif()
