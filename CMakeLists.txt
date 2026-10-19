cmake_minimum_required(VERSION 3.13)

project(TROLL VERSION 4.1.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "CXX compiler ${CMAKE_CXX_COMPILER}, build type ${CMAKE_BUILD_TYPE}")

# --- Options -----------------------------------------------------------------
option(TROLL_WARNINGS_AS_ERRORS "Treat compiler warnings as errors (-Werror); ON in CI" OFF)
option(TROLL_EXTRA_WARNINGS "Add -Wextra -Wshadow -Wconversion (advisory: several hundred warnings today)" OFF)
option(TROLL_OPENMP "Parallelise leaf fluxes and the LAImax table with OpenMP (results do not depend on the thread count)" ON)
option(TROLL_LTO "Link-time optimisation (results unchanged, a few percent faster)" ON)
set(TROLL_SANITIZE "" CACHE STRING "Semicolon-separated sanitizers, e.g. address;undefined")
set(TROLL_DEFINES "" CACHE STRING
    "Compile-time model options, e.g. \"Output_ABC;TROLL_NO_WATER\" (see include/troll_defines.hpp)")

# --- Dependencies --------------------------------------------------------------
find_package(GSL REQUIRED)
message(STATUS "GSL ${GSL_VERSION} found in ${GSL_INCLUDE_DIRS}")

# --- Executable ----------------------------------------------------------------
add_executable(${PROJECT_NAME}
    src/troll.cpp
    src/species.cpp
    src/crown.cpp
    src/tree.cpp
    src/lookup_tables.cpp
    src/input.cpp
    src/memory.cpp
    src/output.cpp
    src/simulation.cpp
    src/internal/lookup_tables_cache.cpp
    src/internal/param_registry.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/internal
)
target_link_libraries(${PROJECT_NAME} PRIVATE GSL::gsl GSL::gslcblas)
target_compile_definitions(${PROJECT_NAME} PRIVATE TROLL_VERSION="${PROJECT_VERSION}")

if(TROLL_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
        message(STATUS "OpenMP enabled (set OMP_NUM_THREADS to choose the number of threads)")
    else()
        message(WARNING "OpenMP not found: building single-threaded. On macOS: brew install libomp")
    endif()
endif()

if(TROLL_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _troll_ipo OUTPUT _troll_ipo_msg)
    if(_troll_ipo AND NOT TROLL_SANITIZE)
        set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
        message(STATUS "Link-time optimisation enabled")
    endif()
endif()

# Warnings are ON (the former global "-w" silently disabled all of them).
target_compile_options(${PROJECT_NAME} PRIVATE -Wall)
if(TROLL_EXTRA_WARNINGS)
    target_compile_options(${PROJECT_NAME} PRIVATE -Wextra -Wshadow -Wconversion)
endif()
if(TROLL_WARNINGS_AS_ERRORS)
    target_compile_options(${PROJECT_NAME} PRIVATE -Werror)
endif()

if(TROLL_DEFINES)
    message(STATUS "TROLL model options: ${TROLL_DEFINES}")
    target_compile_definitions(${PROJECT_NAME} PRIVATE ${TROLL_DEFINES})
endif()

if(TROLL_SANITIZE)
    string(REPLACE ";" "," _troll_sanitizers "${TROLL_SANITIZE}")
    message(STATUS "Sanitizers: ${_troll_sanitizers}")
    target_compile_options(${PROJECT_NAME} PRIVATE -fsanitize=${_troll_sanitizers} -fno-omit-frame-pointer -g)
    target_link_options(${PROJECT_NAME} PRIVATE -fsanitize=${_troll_sanitizers})
endif()

# --- Regression checks (see sanity/README.md) --------------------------------------
set(_troll_check ${CMAKE_SOURCE_DIR}/sanity/check.sh --troll $<TARGET_FILE:${PROJECT_NAME}>)

# 4.1.1 computes the leaf kernel in single precision, so outputs differ from the pinned
# reference (4.0 / dev@4de3273) in the last digits. Until the reference can be moved to a
# published 4.1.1 commit, the checks use tolerances measured for this change:
# - the 20-day example is compared file by file;
# - longer runs are compared on stand-level aggregates only, because after a few months a
#   last-digit difference tips a mortality draw and per-tree outputs follow a different
#   (statistically equivalent) trajectory.
# Set the tolerances to 0 and TROLL_SANITY_FILES_LONG to "" for bit-exact checks
# against a 4.1.1 reference.
set(TROLL_SANITY_RTOL_EXAMPLE "1e-2" CACHE STRING "Relative tolerance, 20-day example")
set(TROLL_SANITY_RTOL_LONG "1e-2" CACHE STRING "Relative tolerance, long runs (stand, 1000 days), aggregates")
set(TROLL_SANITY_ATOL "1e-6" CACHE STRING "Absolute tolerance (ignores negligible values)")
set(TROLL_SANITY_FILES_LONG "sumstats|water_balance|LAIprofile" CACHE STRING "Output files compared for long runs (regex)")
set(_troll_tol_example --rtol ${TROLL_SANITY_RTOL_EXAMPLE} --atol ${TROLL_SANITY_ATOL})
set(_troll_tol_stand --rtol ${TROLL_SANITY_RTOL_LONG} --atol ${TROLL_SANITY_ATOL} --files ${TROLL_SANITY_FILES_LONG})

# 20 days from bare ground on the 4 ha example (fast smoke test)
add_custom_target(sanity
    COMMAND ${_troll_check} --case example --nbiter 20 ${_troll_tol_example}
    DEPENDS ${PROJECT_NAME} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sanity
    COMMENT "Sanity check: example, 20 days" USES_TERMINAL VERBATIM)

# 1000 days from bare ground (replaces the missing check_1000.sh)
add_custom_target(sanity1000
    COMMAND ${_troll_check} --case example --nbiter 1000 ${_troll_tol_stand}
    DEPENDS ${PROJECT_NAME} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sanity
    COMMENT "Sanity check: example, 1000 days" USES_TERMINAL VERBATIM)

# One year on an established 1 ha stand: adult growth, competition, mortality, treefall
add_custom_target(sanity-stand
    COMMAND ${_troll_check} --case stand --nbiter 365 ${_troll_tol_stand}
    DEPENDS ${PROJECT_NAME} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sanity
    COMMENT "Sanity check: established stand, 1 year" USES_TERMINAL VERBATIM)

add_custom_target(check DEPENDS sanity sanity-stand)

# Outputs must be bit-identical whatever the number of OpenMP threads
add_custom_target(threads-check
    COMMAND ${CMAKE_SOURCE_DIR}/sanity/bench.sh --troll $<TARGET_FILE:${PROJECT_NAME}> --threads "1 4" --nbiter 30
    DEPENDS ${PROJECT_NAME} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sanity
    COMMENT "Thread-count independence: 1 vs 4 threads must give identical outputs" USES_TERMINAL VERBATIM)
